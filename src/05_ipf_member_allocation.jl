# =============================================================================
# 05_ipf_member_allocation.jl
# Sequential household member allocation via IPF and integer rounding (MILP)
#
# Implements Sections 2.2–2.3 and Algorithm 1 of the paper:
#   - Seed matrix construction from 2010 census microdata
#   - Column marginal adjustment (Algorithm 1)
#   - Iterative Proportional Fitting (Section 2.2.3)
#   - Integer rounding via MILP (P1) solved with HiGHS (Section 2.3)
#   - Member-to-household assignment (adjourning orders k = 2, …, 15)
# =============================================================================

using CSV, DataFrames, StatsBase, Random, LinearAlgebra, JuMP, HiGHS

include(joinpath(@__DIR__, "..", "config.jl"))
include(joinpath(@__DIR__, "utils.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# IPF: Iterative Proportional Fitting (Section 2.2.3)
# ─────────────────────────────────────────────────────────────────────────────

"""
    ipf(seed, R, C; max_iter=1000, tol=1e-3) -> Matrix{Float64}

Standard IPF (Deming & Stephan, 1940). Starting from `seed`, alternately
scale rows to match marginals `R` and columns to match marginals `C`.

Convergence is declared when max(|row_sums − R|, |col_sums − C|) < `tol`.
"""
function ipf(seed::Matrix{Float64}, R::Vector, C::Vector;
             max_iter::Int=1000, tol::Float64=1e-3)
    X = copy(seed)
    total_target = sum(R)
    s = sum(X)
    if s > 0
        X .*= (total_target / s)
    end

    for iter in 1:max_iter
        # Row scaling
        rs = sum(X, dims=2)
        X .= X .* (R ./ [s > 0 ? s : 1.0 for s in rs])

        # Column scaling
        cs = sum(X, dims=1)
        X .= X .* (C' ./ [s > 0 ? s : 1.0 for s in cs])

        # Convergence check (maximum absolute error)
        err_R = maximum(abs.(sum(X, dims=2) .- R))
        err_C = maximum(abs.(sum(X, dims=1) .- C'))
        if err_R < tol && err_C < tol
            return X
        end
    end

    @warn "IPF did not converge. Final error: $(maximum(abs.(sum(X, dims=2) .- R)))"
    return X
end

# ─────────────────────────────────────────────────────────────────────────────
# Integer Rounding via MILP (P1) — Section 2.3
# ─────────────────────────────────────────────────────────────────────────────

"""
    integer_round(X_star, row_marginals, col_marginals) -> Matrix{Int}

Solve the MILP (P1) that finds an integer matrix Y minimizing the total
absolute deviation from the continuous IPF solution `X_star`, subject to
exact row and column marginal constraints.

Uses HiGHS as the solver.
"""
function integer_round(X_star::Matrix{Float64},
                       row_marginals::Vector{Int},
                       col_marginals::Vector{Int})
    N = size(X_star, 1)

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    set_attribute(model, "threads", 1)
    set_attribute(model, "random_seed", 42)

    @variable(model, Y[1:N, 1:N] >= 0, Int)

    # Row marginal constraints (P1b)
    @constraint(model, [i=1:N], sum(Y[i, j] for j in 1:N) == row_marginals[i])
    # Column marginal constraints (P1c)
    @constraint(model, [j=1:N], sum(Y[i, j] for i in 1:N) == col_marginals[j])

    # Absolute deviation linearization (P1d)
    @variable(model, t[1:N, 1:N] >= 0)
    for i in 1:N, j in 1:N
        @constraint(model, Y[i, j] - X_star[i, j] <= t[i, j])
        @constraint(model, X_star[i, j] - Y[i, j] <= t[i, j])
    end

    # Minimize total deviation (P1a)
    @objective(model, Min, sum(t))

    optimize!(model)

    return round.(Int, value.(Y))
end

# ─────────────────────────────────────────────────────────────────────────────
# Column Marginal Adjustment — Algorithm 1
# ─────────────────────────────────────────────────────────────────────────────

"""
    adjust_column_marginals(seed, row_marginals, stock) -> Vector{Int}

Compute adjusted column marginals following Algorithm 1 of the paper.

Given the seed matrix `seed`, desired row marginals `row_marginals`, and
available population stock per type `stock`, this function:
1. Computes an expansion factor to scale seed column sums to the total demand.
2. Caps each column marginal by the available stock.
3. Redistributes any residual cyclically to types with remaining stock.

Returns a vector C satisfying sum(C) == sum(row_marginals).
"""
function adjust_column_marginals(seed::Matrix{Float64},
                                  row_marginals::Vector{Int},
                                  stock::Vector{Int})
    N = size(seed, 1)

    # Step 1: Seed column sums (V_j)
    V = vec(sum(seed, dims=1))

    # Step 2: Expansion factor
    total_demand = sum(row_marginals)
    total_seed = sum(seed)
    expansion_factor = total_demand / total_seed

    # Step 3: Expanded demand per type
    demand = V .* expansion_factor

    # Step 4: Initial allocation capped by stock
    C = Int.(floor.(min.(demand, Float64.(stock))))

    # Step 5: Redistribute residual
    residual = total_demand - sum(C)
    if residual > 0
        remaining_stock = stock .- C
        sorted_idx = sortperm(remaining_stock, rev=true)
        for ℓ in 1:residual
            idx = sorted_idx[mod1(ℓ, N)]
            if remaining_stock[idx] > 0
                C[idx] += 1
                remaining_stock[idx] -= 1
            end
        end
    end

    return C
end

# ─────────────────────────────────────────────────────────────────────────────
# Member Assignment — Translates integer matrix Y* into household pairings
# ─────────────────────────────────────────────────────────────────────────────

"""
    assign_members!(ord, Y, group_2022, DFP, type_labels)

For each cell Y[j, i] > 0, select Y[j, i] heads of type j and Y[j, i]
unassigned individuals of type i, and assign each individual to the
corresponding head's household with order `ord`.

Modifies `DFP` in place (columns :Ordem and :ID_Familia).
"""
function assign_members!(ord::Int, Y::Matrix{Int},
                         group_2022::DataFrame, DFP::DataFrame,
                         type_labels::Vector{String})
    rng = MersenneTwister(123 + ord)
    N = length(type_labels)

    # Build ID-to-row index for fast lookup
    id_to_row = Dict(id => i for (i, id) in enumerate(DFP.ID))

    # Initialize pools by demographic type
    head_pool = Dict{String, Vector{Tuple{Any, Any}}}()
    candidate_pool = Dict{String, Vector{Any}}()
    for t in type_labels
        head_pool[t] = Tuple{Any, Any}[]
        candidate_pool[t] = Any[]
    end

    # Populate pools from the municipal subset
    for row in eachrow(group_2022)
        if row.Ordem == 1 && row.N_M >= ord
            push!(head_pool[row.Tipo], (row.ID, row.ID_Familia))
        elseif row.Ordem == 0
            push!(candidate_pool[row.Tipo], row.ID)
        end
    end

    # Shuffle pools for randomized tie-breaking
    for t in type_labels
        shuffle!(rng, head_pool[t])
        shuffle!(rng, candidate_pool[t])
    end

    # Main assignment loop following the cardinality matrix Y
    for j in 1:N
        head_type = type_labels[j]
        available_heads = head_pool[head_type]

        if isempty(available_heads) || sum(Y[j, :]) <= 0
            continue
        end

        head_ptr = 1

        for i in 1:N
            member_type = type_labels[i]
            qty = Y[j, i]
            qty <= 0 && continue

            n_members = length(candidate_pool[member_type])
            n_heads = length(available_heads) - head_ptr + 1
            actual = min(qty, n_members, n_heads)
            actual <= 0 && continue

            # Select heads
            head_slice = available_heads[head_ptr:(head_ptr + actual - 1)]
            family_ids = [x[2] for x in head_slice]

            # Select and remove candidates from pool
            member_ids = splice!(candidate_pool[member_type], 1:actual)

            # Update DFP
            rows_to_edit = [id_to_row[id] for id in member_ids]
            DFP.Ordem[rows_to_edit] .= ord
            DFP.ID_Familia[rows_to_edit] .= family_ids

            head_ptr += actual
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Load and prepare 2010 Census microdata
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_microdata_2010() -> DataFrame

Load and prepare the 2010 Census microdata sample (D). Selects relevant
columns, recodes gender and ethnicity, and adds the age bracket variable.
"""
function load_microdata_2010()
    DF10 = CSV.read(joinpath(PATH_RAW_MICRO, "pop_filtrado.csv"), DataFrame)

    cols = [:code_muni, :code_weighting, :V0300, :V0010, :V0502, :V0504,
            :V0601, :V6036, :V0606]
    DF10 = DF10[:, cols]
    rename!(DF10, [:mun, :area_pod, :cod_res, :peso, :Relacao_PR, :Ordem,
                   :genero, :idade, :raça])

    replace!(DF10.genero, "Feminino" => "F", "Masculino" => "M")

    race_map = Dict("Branca" => 1, "Preta" => 2, "Amarela" => 3,
                    "Parda" => 4, "Indígena" => 5)
    DF10.raça = get.(Ref(race_map), DF10.raça, 0)

    return DF10
end

# ─────────────────────────────────────────────────────────────────────────────
# Main entry point: allocate members of order k
# ─────────────────────────────────────────────────────────────────────────────

"""
    allocate_order(ord::Int)

Run the full IPF + integer rounding + assignment pipeline for household
member order `ord` (k = 2, …, 15). Reads the person dataset with orders
up to `ord - 1` already assigned, and writes the updated dataset.

For orders k ≥ 7, uses the extended family dataset that includes households
with more than 6 members (see `expand_large_households`).
"""
function allocate_order(ord::Int)
    # Load family dataset (extended version for k >= 7)
    family_file = ord >= 7 ? "dataset_familia_extended.csv" : "dataset_familia.csv"
    DFF = CSV.read(joinpath(PATH_STAGE3, family_file), DataFrame)

    # Load person dataset with orders up to (ord - 1) assigned
    prev = ord - 1
    DFP = CSV.read(joinpath(PATH_STAGE3, "dataset_pessoas_ordem$(prev).csv"), DataFrame)
    DFP.faixa = age_bracket.(DFP.idade)

    # Load and prepare 2010 microdata
    DF10 = load_microdata_2010()
    DF10.idade = ifelse.(DF10.idade .>= 100, 100, DF10.idade)
    DF10.faixa = age_bracket.(DF10.idade)

    # Handle municipalities split between 2010 and 2022
    mun_10 = unique(DF10.mun)
    mun_22 = unique(DFF.mun)
    common = filter(x -> x in mun_22, mun_10)

    # Balneário Rincão (from Içara) and Pescaria Brava (from Laguna)
    split_pairs = [4207007 4209409]  # [parent_mun new_mun]
    split_new = setdiff(mun_22, mun_10)

    # Replicate parent municipality data for newly created municipalities
    for i in 1:size(split_pairs, 1)
        parent = split_pairs[i, 1]
        child = split_new[i]  # assumes same ordering
        rows_copy = deepcopy(DF10[DF10.mun .== parent, :])
        rows_copy.mun .= child
        append!(DF10, rows_copy)
    end

    groups_10 = groupby(DF10, :mun)

    # Process each municipality
    for n_group in 1:295
        group = groups_10[n_group]
        mun_code = group.mun[1]
        println("Order $ord — municipality $n_group ($mun_code)")

        # Build demographic type labels for the 2022 subset
        subset_22 = DFP[DFP.mun .== mun_code, :]
        subset_22.Tipo = string.(subset_22.genero, "_", subset_22.raça, "_",
                                 subset_22.faixa)

        # Attach household size (N_M) to heads for filtering
        subset_22 = leftjoin(subset_22, select(DFF, :ID_Familia, :N_M),
                             on=:ID_Familia)

        type_labels = sort(unique(subset_22.Tipo))
        type_index = Dict(type_labels .=> 1:length(type_labels))
        subset_22.seqtipo = [type_index[t] for t in subset_22.Tipo]

        # Check if there are households requiring order k
        mask = (DFF.mun .== mun_code) .& (DFF.N_M .>= ord)
        if sum(mask) > 0
            Y = build_cardinality_matrix(ord, group, mun_code, DFF,
                                          subset_22, type_labels)
            assign_members!(ord, Y, subset_22, DFP, type_labels)
        end
    end

    # Remove temporary column and save
    select!(DFP, Not(:faixa))
    CSV.write(joinpath(PATH_STAGE3, "dataset_pessoas_ordem$(ord).csv"), DFP)
end

"""
    build_cardinality_matrix(ord, group_10, mun_code, DFF, group_22, type_labels)
        -> Matrix{Int}

Construct the seed matrix from 2010 microdata, compute row/column marginals,
run IPF, and solve the integer rounding problem (P1).

Returns the integer cardinality matrix Y* of size N × N, where N is the
number of demographic types.
"""
function build_cardinality_matrix(ord::Int, group_10::AbstractDataFrame,
                                   mun_code::Int, DFF::DataFrame,
                                   group_22::DataFrame,
                                   type_labels::Vector{String})
    N = length(type_labels)
    type_index = Dict(type_labels .=> 1:N)

    # Build demographic type for 2010 data
    group_10 = copy(DataFrame(group_10))
    group_10.Tipo = string.(group_10.genero, "_", group_10.raça, "_",
                            group_10.faixa)

    # Keep only types that exist in 2022
    allowed = Set(type_labels)
    group_10 = filter(row -> row.Tipo in allowed, group_10)
    group_10.seqtipo = [type_index[t] for t in group_10.Tipo]

    # Heads (order 1) and members of target order in 2010 sample
    df_heads = filter(row -> row.Ordem == 1, group_10)
    df_members = filter(row -> row.Ordem == ord, group_10)

    # Join heads with members by household to build seed
    heads_join = select(df_heads, :cod_res, :seqtipo => :seqtype_head, :peso)
    members_join = select(df_members, :cod_res, :seqtipo => :seqtype_member)
    pairs = innerjoin(heads_join, members_join, on=:cod_res)

    observed = combine(groupby(pairs, [:seqtype_head, :seqtype_member]),
                       :peso => sum => :weight)

    # Build seed matrix with smoothing
    SMOOTHING = 0.001
    seed = fill(SMOOTHING, N, N)
    for row in eachrow(observed)
        seed[Int(row.seqtype_head), Int(row.seqtype_member)] = row.weight
    end

    # Zero out rows for heads below minimum age (bracket 1 = ages 0–13)
    underage_idx = findall(k -> endswith(k, "_1"), type_labels)
    seed[underage_idx, :] .= 0.0

    # Row marginals: heads in 2022 with N_M >= ord
    heads_22 = filter(row -> row.Ordem == 1 && row.N_M >= ord, group_22)
    row_counts = combine(groupby(heads_22, :seqtipo), nrow => :count)
    skeleton = DataFrame(seqtipo=1:N)
    row_marg_df = leftjoin(skeleton, row_counts, on=:seqtipo)
    sort!(row_marg_df, :seqtipo)
    row_marginals = coalesce.(row_marg_df.count, 0)

    # Column marginals: unassigned individuals in 2022
    unassigned_22 = filter(row -> row.Ordem == 0, group_22)
    col_counts = combine(groupby(unassigned_22, :seqtipo), nrow => :count)
    col_marg_df = leftjoin(skeleton, col_counts, on=:seqtipo)
    sort!(col_marg_df, :seqtipo)
    stock = coalesce.(col_marg_df.count, 0)

    # Adjust column marginals via Algorithm 1
    col_marginals = adjust_column_marginals(Float64.(seed), row_marginals, stock)

    # Run IPF
    X_star = ipf(Float64.(seed), row_marginals, col_marginals)

    # Integer rounding (P1)
    Y = integer_round(X_star, row_marginals, col_marginals)

    return Y
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
# To run for all orders:
#   for k in 2:6; allocate_order(k); end
#   expand_large_households()   # see Section 4.2.3
#   for k in 7:15; allocate_order(k); end
