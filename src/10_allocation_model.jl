# =============================================================================
# 10_allocation_model.jl
# Household-to-residence allocation via MILP / LP (Section 3)
#
# Implements the allocation model (P2) from the paper:
#   - Compatibility-based class grouping of families and residences
#   - MILP formulation with composite objective (quality + volume − penalty)
#   - Solved as LP relaxation (Corollary 1: integrality is guaranteed)
#   - Uses HiGHS with the primal simplex method
# =============================================================================

using CSV, DataFrames, StatsBase, JuMP, HiGHS, Dates

include(joinpath(@__DIR__, "..", "config.jl"))
include(joinpath(@__DIR__, "utils.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Census gender marginals — spatial targets M_ℓ (Section 5.1)
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_gender_marginals(df_sectors, df_aggregates)
        -> (Dict{(String,Int),Int}, Dict{(String,Int),Int})

Pre-compute the census target counts M_ℓ for male-headed and female-headed
households, indexed by (census_tract, bathroom_category).

Returns two dictionaries: `male_targets` and `female_targets`.
"""
function load_gender_marginals(df_sectors::DataFrame, df_aggregates::DataFrame)
    male_targets = Dict{Tuple{String, Int}, Int}()
    female_targets = Dict{Tuple{String, Int}, Int}()

    for row in eachrow(df_sectors)
        row.N_RPO == 0 && continue
        agg_row = df_aggregates[df_aggregates.setor .== row.COD_SETOR_NOVO, :]
        isempty(agg_row) && continue
        marginal = agg_row[1, :]

        for nb in 1:7
            col_offset = 295 + (nb - 1) * 2
            col_male = Symbol("V$(lpad(col_offset, 5, '0'))")
            col_female = Symbol("V$(lpad(col_offset + 1, 5, '0'))")
            male_targets[(row.COD_SETOR, nb)] =
                max(0, Int(to_number(marginal[col_male])))
            female_targets[(row.COD_SETOR, nb)] =
                max(0, Int(to_number(marginal[col_female])))
        end
    end

    return male_targets, female_targets
end

# ─────────────────────────────────────────────────────────────────────────────
# Solution unpacking — Map class-level decisions to individual IDs
# ─────────────────────────────────────────────────────────────────────────────

"""
    unpack_solution(pairs, y, Kf, Kh, families_by_class, residences_by_class)
        -> DataFrame

Translate the optimization solution `y` (number of allocations per class pair)
into a DataFrame of individual (ID_Familia, ID_R) pairings.
"""
function unpack_solution(pairs, y, Kf, Kh, families_by_class, residences_by_class)
    ids_fam = Int[]
    ids_res = String[]

    # Pointer-based unpacking (avoids popfirst! mutation issues)
    fam_ptr = Dict(k => 1 for k in keys(families_by_class))
    res_ptr = Dict(k => 1 for k in keys(residences_by_class))

    # Process pairs in decreasing allocation order for consistency
    sorted_pairs = sort(pairs, by=p -> JuMP.value(y[p]), rev=true)

    for p in sorted_pairs
        qty = round(Int, JuMP.value(y[p]))
        qty <= 0 && continue

        fk = Kf[p[1]]
        rk = Kh[p[2]]
        fam_vec = families_by_class[fk]
        res_vec = residences_by_class[rk]

        (isempty(fam_vec) || isempty(res_vec)) && continue

        for _ in 1:qty
            fi = min(fam_ptr[fk], length(fam_vec))
            ri = min(res_ptr[rk], length(res_vec))
            push!(ids_fam, fam_vec[fi])
            push!(ids_res, res_vec[ri])
            fam_ptr[fk] += 1
            res_ptr[rk] += 1
        end
    end

    return DataFrame(ID_Familia=ids_fam, ID_R=ids_res)
end

# ─────────────────────────────────────────────────────────────────────────────
# Core optimization — Model (P2)
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_allocation(df_residences, df_families, male_targets, female_targets)
        -> DataFrame

Solve the household-to-residence allocation model (P2) for a single
municipality. By Corollary 1, the LP relaxation yields integral solutions,
so the model is solved as a continuous LP using the primal simplex method.

Returns a DataFrame of (ID_Familia, ID_R) pairings.
"""
function solve_allocation(df_residences::AbstractDataFrame,
                          df_families::DataFrame,
                          male_targets::Dict, female_targets::Dict)
    # ── Compute tract-level target ratio τ_s ──────────────────────────────
    mean_nm = isempty(df_families) ? 0.0 : mean(df_families.N_M)

    mask_valid = (df_residences.ESPECIE .== 1) .&
                 (df_residences.TIPO_ESPECIE .!= 104) .&
                 (df_residences.N_B .∈ Ref(1:7))
    valid_caps = effective_capacity.(df_residences.N_B[mask_valid])
    mean_bt = isempty(valid_caps) ? 1.0 : mean(valid_caps)
    τ_local = mean_bt > 0 ? (mean_nm / mean_bt) : 2.0

    # ── Group residences into classes (subtype, ethnicity, tract, NB) ─────
    residences_by_class = Dict{Tuple{Any,Any,Any,Int}, Vector{String}}()
    for row in eachrow(df_residences)
        if row.ESPECIE == 1 && row.TIPO_ESPECIE != 104 && 1 <= row.N_B <= 7
            k = (row.TIPO_ESPECIE, row.PR_R, row.SETOR, Int(row.N_B))
            push!(get!(residences_by_class, k, String[]), string(row.ID))
        end
    end

    # ── Group families into classes (subtype, ethnicity, gender, NM) ──────
    families_by_class = Dict{Tuple{Any,Any,Int,Int}, Vector{Int}}()
    for row in eachrow(df_families)
        k = (row.TIPO_ESPECIE, row.PR_R, gender_code(row.PR_G), Int(row.N_M))
        push!(get!(families_by_class, k, Int[]), row.ID_Familia)
    end

    Kf = collect(keys(families_by_class))
    Kh = collect(keys(residences_by_class))

    if isempty(Kf) || isempty(Kh)
        return DataFrame(ID_Familia=Int64[], ID_R=String[])
    end

    # ── Admissible pairs: matching subtype and head ethnicity (Eq. 1) ─────
    pairs = [(i, j) for i in eachindex(Kf) for j in eachindex(Kh)
             if Kf[i][1] == Kh[j][1] && Kf[i][2] == Kh[j][2]]

    if isempty(pairs)
        return DataFrame(ID_Familia=Int64[], ID_R=String[])
    end

    sup = Dict(k => length(v) for (k, v) in families_by_class)
    cap = Dict(k => length(v) for (k, v) in residences_by_class)

    # Unique (tract, NB) combinations for spatial target constraints
    tract_nb_combos = unique([(Kh[j][3], Kh[j][4]) for j in eachindex(Kh)])

    # ── Build and solve LP (P2 relaxation) ────────────────────────────────
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    set_attribute(model, "solver", "simplex")  # Primal simplex (Section 5.2)

    # Decision variables: y^r_f >= 0 (continuous by Corollary 1)
    @variable(model, y[p in pairs] >= 0)
    # Slack variables: δ_ℓ >= 0 (P2e)
    @variable(model, δ[tract_nb_combos, 1:2] >= 0)

    # Constraint (P2b): family class capacity
    for i in eachindex(Kf)
        ps = [p for p in pairs if p[1] == i]
        isempty(ps) || @constraint(model, sum(y[p] for p in ps) <= sup[Kf[i]])
    end

    # Constraint (P2c): residence class capacity
    for j in eachindex(Kh)
        ps = [p for p in pairs if p[2] == j]
        isempty(ps) || @constraint(model, sum(y[p] for p in ps) <= cap[Kh[j]])
    end

    # Constraint (P2d): spatial target adherence with slack
    for (s, nb) in tract_nb_combos, g in 1:2
        ps = [p for p in pairs
              if Kf[p[1]][3] == g && Kh[p[2]][3] == s && Kh[p[2]][4] == nb]
        isempty(ps) && continue
        target_dict = (g == 1) ? male_targets : female_targets
        M_ℓ = get(target_dict, (s, nb), 0)
        @constraint(model, sum(y[p] for p in ps) <= M_ℓ + δ[(s, nb), g])
    end

    # Objective (P2a): maximize weighted allocation − penalty for slack
    α = 1.0
    β = 1.0
    @objective(model, Max,
        sum((α + allocation_quality(Kf[p[1]][4], Kh[p[2]][4], τ_local)) * y[p]
            for p in pairs)
        - β * sum(δ[(s, nb), g]
            for (s, nb) in tract_nb_combos for g in 1:2))

    optimize!(model)

    return unpack_solution(pairs, y, Kf, Kh, families_by_class,
                           residences_by_class)
end

# ─────────────────────────────────────────────────────────────────────────────
# Pipeline entry point — Process all 295 municipalities
# ─────────────────────────────────────────────────────────────────────────────

function run_allocation()
    println("Loading datasets...")
    DFF = CSV.read(joinpath(PATH_STAGE3, "dataset_familias_2025.csv"), DataFrame)
    DFR = CSV.read(joinpath(PATH_STAGE4, "residencias_2025.csv"), DataFrame)
    df_sectors = CSV.read(joinpath(PATH_STAGE2, "correspondencia_setores.csv"),
                          DataFrame)
    df_agg = CSV.read(joinpath(PATH_RAW_AGGR,
                "Agregados_por_setores_caracteristicas_domicilio2_SC.csv"),
                      DataFrame)

    ensure_mun_column!(DFR)
    ensure_mun_column!(DFF)

    println("Pre-computing census gender marginals...")
    male_targets, female_targets = load_gender_marginals(df_sectors, df_agg)

    # Initialize residence ID column
    DFF[!, :ID_R] = Vector{Union{Missing, String}}(missing, nrow(DFF))

    # Group by municipality, process smallest first
    mun_groups = groupby(DFR, :mun)
    total_muns = length(mun_groups)

    mun_sizes = [(idx=i, mun=g.mun[1], n=count(DFF.mun .== g.mun[1]))
                 for (i, g) in enumerate(mun_groups)]
    sort!(mun_sizes, by=x -> x.n)

    results = Vector{DataFrame}(undef, total_muns)
    partials_dir = joinpath(PATH_RESULTS, "partials")
    mkpath(partials_dir)

    t_start = now()

    for info in mun_sizes
        mun_code = info.mun
        partial_file = joinpath(partials_dir, "muni_$(mun_code).csv")

        # Skip if already computed
        if isfile(partial_file)
            results[info.idx] = CSV.read(partial_file, DataFrame)
            println("Municipality $mun_code loaded from cache.")
            continue
        end

        df_mun = mun_groups[info.idx]
        fam_mun = DFF[DFF.mun .== mun_code, :]

        if nrow(df_mun) > 0 && nrow(fam_mun) > 0
            df_res = solve_allocation(df_mun, fam_mun,
                                      male_targets, female_targets)
            CSV.write(partial_file, df_res)
            results[info.idx] = df_res

            pct = info.n > 0 ? round(nrow(df_res) / info.n * 100, digits=2) : 0.0
            println("Municipality $mun_code: $(nrow(df_res)) allocations ($pct%)")
        else
            empty_df = DataFrame(ID_Familia=Int64[], ID_R=String[])
            CSV.write(partial_file, empty_df)
            results[info.idx] = empty_df
        end
    end

    # Consolidate results
    println("Consolidating results...")
    valid = [results[i] for i in 1:total_muns if isassigned(results, i)]
    all_pairs = reduce(vcat, valid;
                       init=DataFrame(ID_Familia=Int64[], ID_R=String[]))

    if nrow(all_pairs) > 0
        lookup = Dict(r.ID_Familia => r.ID_R for r in eachrow(all_pairs))
        for row in eachrow(DFF)
            haskey(lookup, row.ID_Familia) && (row.ID_R = lookup[row.ID_Familia])
        end
    end

    # Clean up and save
    :mun in propertynames(DFR) && select!(DFR, Not(:mun))
    :mun in propertynames(DFF) && select!(DFF, Not(:mun))

    output_file = joinpath(PATH_RESULTS, "dataset_familias_com_ID_R.csv")
    CSV.write(output_file, DFF)

    elapsed = canonicalize(Dates.CompoundPeriod(now() - t_start))
    println("Allocation complete. Total time: $elapsed")
    println("Output: $output_file")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_allocation()
end
