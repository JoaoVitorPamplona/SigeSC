# =============================================================================
# benchmark_lp_vs_milp.jl
# Computational benchmark: LP relaxation vs. MILP formulation (Table 1)
#
# For the 10 largest municipalities in Santa Catarina, this script:
#   1. Solves the allocation model (P2) as an LP (simplex)
#   2. Solves the same model as a MILP (branch-and-bound)
#   3. Verifies that the LP solution is naturally integral (Theorem 2)
#   4. Reports solve times and speedup factors
#
# All instances use HiGHS as the solver.
# =============================================================================

using CSV, DataFrames, StatsBase, JuMP, HiGHS, Dates

include(joinpath(@__DIR__, "..", "config.jl"))
include(joinpath(@__DIR__, "..", "src", "utils.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

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

# ─────────────────────────────────────────────────────────────────────────────
# Pre-compute census gender marginals
# ─────────────────────────────────────────────────────────────────────────────

male_targets = Dict{Tuple{String, Int}, Int}()
female_targets = Dict{Tuple{String, Int}, Int}()

for row in eachrow(df_sectors)
    row.N_RPO == 0 && continue
    agg_row = df_agg[df_agg.setor .== row.COD_SETOR_NOVO, :]
    isempty(agg_row) && continue
    marginal = agg_row[1, :]
    for nb in 1:7
        offset = 295 + (nb - 1) * 2
        male_targets[(row.COD_SETOR, nb)] =
            max(0, Int(to_number(marginal[Symbol("V$(lpad(offset, 5, '0'))")])))
        female_targets[(row.COD_SETOR, nb)] =
            max(0, Int(to_number(marginal[Symbol("V$(lpad(offset+1, 5, '0'))")])))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Benchmark function: solve both LP and MILP for one municipality
# ─────────────────────────────────────────────────────────────────────────────

function benchmark_municipality(df_mun, fam_mun)
    # Compute τ_local
    mean_nm = isempty(fam_mun) ? 0.0 : mean(fam_mun.N_M)
    mask_bt = (df_mun.ESPECIE .== 1) .& (df_mun.TIPO_ESPECIE .!= 104) .&
              (df_mun.N_B .∈ Ref(1:7))
    valid_caps = effective_capacity.(df_mun.N_B[mask_bt])
    mean_bt = isempty(valid_caps) ? 1.0 : mean(valid_caps)
    τ_local = mean_bt > 0 ? (mean_nm / mean_bt) : 2.0

    # Class grouping
    res_classes = Dict{Tuple{Any,Any,Any,Int}, Vector{String}}()
    for r in eachrow(df_mun)
        if r.ESPECIE == 1 && r.TIPO_ESPECIE != 104 && 1 <= r.N_B <= 7
            push!(get!(res_classes, (r.TIPO_ESPECIE, r.PR_R, r.SETOR, Int(r.N_B)),
                       String[]), string(r.ID))
        end
    end

    fam_classes = Dict{Tuple{Any,Any,Int,Int}, Vector{Int}}()
    for r in eachrow(fam_mun)
        push!(get!(fam_classes, (r.TIPO_ESPECIE, r.PR_R, gender_code(r.PR_G),
                                Int(r.N_M)), Int[]), r.ID_Familia)
    end

    Kf = collect(keys(fam_classes))
    Kh = collect(keys(res_classes))
    pairs = [(i, j) for i in eachindex(Kf) for j in eachindex(Kh)
             if Kf[i][1] == Kh[j][1] && Kf[i][2] == Kh[j][2]]
    combos = unique([(Kh[j][3], Kh[j][4]) for j in eachindex(Kh)])

    isempty(pairs) && return missing

    sup = Dict(k => length(v) for (k, v) in fam_classes)
    cap = Dict(k => length(v) for (k, v) in res_classes)
    α, β = 1.0, 1.0

    # ── Helper: build constraints and objective ───────────────────────────
    function build_model!(model, y_var, slack_var)
        for i in eachindex(Kf)
            ps = [p for p in pairs if p[1] == i]
            isempty(ps) || @constraint(model, sum(y_var[p] for p in ps) <= sup[Kf[i]])
        end
        for j in eachindex(Kh)
            ps = [p for p in pairs if p[2] == j]
            isempty(ps) || @constraint(model, sum(y_var[p] for p in ps) <= cap[Kh[j]])
        end
        for (s, nb) in combos, g in 1:2
            ps = [p for p in pairs
                  if Kf[p[1]][3] == g && Kh[p[2]][3] == s && Kh[p[2]][4] == nb]
            isempty(ps) && continue
            tgt = get(g == 1 ? male_targets : female_targets, (s, nb), 0)
            @constraint(model, sum(y_var[p] for p in ps) <= tgt + slack_var[(s,nb),g])
        end
        @objective(model, Max,
            sum((α + allocation_quality(Kf[p[1]][4], Kh[p[2]][4], τ_local)) *
                y_var[p] for p in pairs)
            - β * sum(slack_var[(s,nb),g] for (s,nb) in combos for g in 1:2))
    end

    # ── 1. LP relaxation ──────────────────────────────────────────────────
    m_lp = Model(HiGHS.Optimizer)
    set_silent(m_lp)
    set_attribute(m_lp, "solver", "simplex")
    @variable(m_lp, y_lp[p in pairs] >= 0)
    @variable(m_lp, f_lp[combos, 1:2] >= 0)
    build_model!(m_lp, y_lp, f_lp)

    t_lp = @elapsed optimize!(m_lp)
    lp_integral = all(v -> isapprox(v, round(v), atol=1e-5), value.(y_lp))

    # ── 2. MILP formulation ───────────────────────────────────────────────
    m_milp = Model(HiGHS.Optimizer)
    set_silent(m_milp)
    set_attribute(m_milp, "time_limit", 120.0)
    @variable(m_milp, y_m[p in pairs] >= 0, Int)
    @variable(m_milp, f_m[combos, 1:2] >= 0)
    build_model!(m_milp, y_m, f_m)

    t_milp = @elapsed optimize!(m_milp)
    status_milp = String(Symbol(termination_status(m_milp)))
    gap = try round(relative_gap(m_milp) * 100, digits=4) catch; NaN end

    return (vars=length(pairs), t_lp=round(t_lp, digits=4),
            t_milp=round(t_milp, digits=4), lp_integral=lp_integral,
            status=status_milp, gap_pct=gap)
end

# ─────────────────────────────────────────────────────────────────────────────
# Run benchmark on the 10 largest municipalities
# ─────────────────────────────────────────────────────────────────────────────

mun_groups = groupby(DFR, :mun)
mun_list = [(mun=g.mun[1], n=count(DFF.mun .== g.mun[1]), idx=i)
            for (i, g) in enumerate(mun_groups)]
sort!(mun_list, by=x -> x.n)

# Select 10 largest
test_cities = reverse(mun_list[end-9:end])

results = DataFrame(
    Municipality=Int64[], Families=Int64[], Variables=Int64[],
    LP_Time_s=Float64[], MILP_Time_s=Float64[],
    LP_Integral=Bool[], MILP_Status=String[], MILP_Gap_Pct=Float64[]
)

println("\nBenchmark: LP relaxation vs. MILP — 10 largest municipalities")
println("=" ^ 75)

for city in test_cities
    df_mun = mun_groups[city.idx]
    fam_mun = DFF[DFF.mun .== city.mun, :]

    res = benchmark_municipality(df_mun, fam_mun)
    if !ismissing(res)
        push!(results, (city.mun, city.n, res.vars, res.t_lp, res.t_milp,
                        res.lp_integral, res.status, res.gap_pct))
        speedup = res.t_milp / max(res.t_lp, 0.001)
        println("  $(city.mun) | Families: $(city.n) | Vars: $(res.vars) | " *
                "LP: $(res.t_lp)s | MILP: $(res.t_milp)s | " *
                "Speedup: $(round(speedup, digits=2))x | " *
                "LP integral: $(res.lp_integral)")
    end
end

# Save results
output_file = joinpath(PATH_BENCHMARKS, "benchmark_lp_vs_milp.csv")
CSV.write(output_file, results)
println("\nResults saved to $output_file")

# Summary statistics
println("\n" * "=" ^ 75)
println("SUMMARY")
println("  Average speedup: $(round(mean(results.MILP_Time_s ./ max.(results.LP_Time_s, 0.001)), digits=2))x")
println("  Median speedup:  $(round(median(results.MILP_Time_s ./ max.(results.LP_Time_s, 0.001)), digits=2))x")
println("  All LP solutions integral: $(all(results.LP_Integral))")
