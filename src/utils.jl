# =============================================================================
# utils.jl — Shared utility functions for the SisgeSC pipeline
# =============================================================================

"""
    to_number(x) -> Float64

Attempt to parse `x` as Float64. Returns 0.0 for missing, corrupted, or
non-numeric inputs.
"""
to_number(x) = try parse(Float64, string(x)) catch; 0.0 end

"""
    to_int(x) -> Int

Attempt to convert `x` to Int. Returns 0 for non-integer inputs.
"""
to_int(x) = x isa Int64 ? x : 0

"""
    gender_code(pr_g) -> Int

Map gender identifiers to numeric codes: 1 = Male, 2 = Female.
Accepts String ("M"/"F"), Int (1/2), or String-encoded Int ("1"/"2").
"""
gender_code(pr_g) = (pr_g == "M" || pr_g == 1 || pr_g == "1") ? 1 : 2

"""
    ensure_mun_column!(df::DataFrame) -> DataFrame

Derive a 7-digit IBGE municipality code from the census tract code (SETOR)
and add it as a `:mun` column if not already present.
"""
function ensure_mun_column!(df::DataFrame)
    if !hasproperty(df, :mun)
        hasproperty(df, :SETOR) ||
            error("DataFrame has neither :mun nor :SETOR to derive municipality code.")
        df[!, :mun] = [parse(Int64, first(string(s), 7)) for s in df.SETOR]
    end
    return df
end

"""
    age_bracket(age::Int) -> Int

Map an exact age to one of 9 demographic brackets used in the IPF stage.
Returns the bracket index (1–9).

| Index | Range     |
|-------|-----------|
| 1     | 0–13      |
| 2     | 14–17     |
| 3     | 18–28     |
| 4     | 29–39     |
| 5     | 40–49     |
| 6     | 50–59     |
| 7     | 60–69     |
| 8     | 70–79     |
| 9     | 80–100    |
"""
function age_bracket(age::Int)
    brackets = [(0,13),(14,17),(18,28),(29,39),(40,49),(50,59),(60,69),(70,79),(80,100)]
    for (i, (lo, hi)) in enumerate(brackets)
        lo <= age <= hi && return i
    end
    return length(brackets)  # fallback for age > 100
end

"""
    normalize_fill!(v::AbstractVector) -> Vector

Normalize `v` to sum to 1. Replace any zero entries with a small positive
value (second-smallest nonzero / length) to avoid excluding combinations
a priori, then renormalize.
"""
function normalize_fill!(v)
    using LinearAlgebra: normalize
    v .= normalize(v, 1)
    zeros_idx = findall(x -> x == 0, v)
    if !isempty(zeros_idx)
        v[zeros_idx] .= sort(unique(v))[min(2, end)] / length(v)
    end
    return normalize(v, 1)
end

"""
    effective_capacity(nb::Int) -> Float64

Map a census bathroom category `nb` ∈ {1,…,7} to its effective sanitary
capacity value b̃(nb), as defined in Section 3 of the paper.

Categories 1–4 map to their numeric value; 5 = shared (0.5),
6 = rudimentary (0.25), 7 = absent (0.1).
"""
function effective_capacity(nb::Int)
    nb <= 4 && return Float64(nb)
    nb == 5 && return 0.5
    nb == 6 && return 0.25
    return 0.1  # nb == 7
end

"""
    allocation_quality(nm::Int, nb::Int, τ::Float64) -> Float64

Compute the quality score q(NM, b̃) for assigning a household of size `nm`
to a residence with bathroom category `nb`, given tract-level target ratio `τ`.

    q = 1 / (1 + (NM / b̃ − τ)²)

Returns a value in (0, 1], maximized when NM/b̃ = τ.
"""
function allocation_quality(nm::Int, nb::Int, τ::Float64)
    b_tilde = effective_capacity(nb)
    ratio = nm / b_tilde
    return 1.0 / (1.0 + (ratio - τ)^2)
end
