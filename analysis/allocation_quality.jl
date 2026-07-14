# =============================================================================
# allocation_quality.jl
# Post-optimization analysis: allocation quality and constraint violations
#
# Generates Figures 4 and 5 of the paper:
#   - Figure 4: Distribution of mean allocation quality q per municipality
#   - Figure 5: Distribution of mean demographic constraint violation δ/M
# =============================================================================

using CSV, DataFrames, StatsBase, CairoMakie, Printf

include(joinpath(@__DIR__, "..", "config.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

results_dir = PATH_RESULTS

q_files = filter(f -> startswith(basename(f), "q_valores_"),
                 readdir(results_dir, join=true))
delta_files = filter(f -> startswith(basename(f), "delta_"),
                     readdir(results_dir, join=true))

df_q = reduce(vcat, [CSV.read(f, DataFrame) for f in q_files])
df_delta = reduce(vcat, [CSV.read(f, DataFrame) for f in delta_files])

# ─────────────────────────────────────────────────────────────────────────────
# Municipality size classification
# ─────────────────────────────────────────────────────────────────────────────

"""
    size_class(n) -> Int

Map number of families to a size category for grouped boxplots.
1 = All, 2 = < 5k, 3 = 5k–20k, 4 = 20k–50k, 5 = >= 50k
"""
function size_class(n)
    n < 5000  && return 2
    n < 20000 && return 3
    n < 50000 && return 4
    return 5
end

labels = ["All", "< 5k", "5k-20k", "20k-50k", ">= 50k"]

# ─────────────────────────────────────────────────────────────────────────────
# Allocation quality q (Figure 4)
# ─────────────────────────────────────────────────────────────────────────────

# Weighted mean q per municipality
q_by_mun = combine(
    groupby(df_q, :mun),
    [:q, :alocados] => ((q, a) -> sum(q .* a) / sum(a)) => :q_mean,
    :alocados => sum => :total_allocated
)

# Build data for boxplot (each municipality appears in "All" + its size class)
df_plot_q = DataFrame(q_mean=Float64[], class=Int[])
for row in eachrow(q_by_mun)
    push!(df_plot_q, (q_mean=row.q_mean, class=1))  # "All"
    push!(df_plot_q, (q_mean=row.q_mean, class=size_class(row.total_allocated)))
end

fig1 = Figure(size=(700, 400))
ax1 = Axis(fig1[1, 1],
    xlabel="Number of families per municipality",
    ylabel="Mean allocation quality (q) per municipality",
    title="Distribution of municipal allocation quality by municipality size",
    xticks=(1:5, labels))

boxplot!(ax1, df_plot_q.class, df_plot_q.q_mean,
         color=:steelblue, whiskerwidth=0.5, width=0.6)
hlines!(ax1, [median(q_by_mun.q_mean)], color=(:red, 0.4),
        linestyle=:dash, linewidth=1)

save(joinpath(results_dir, "fig4_boxplot_q.pdf"), fig1)
println("Saved: fig4_boxplot_q.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# Constraint violation δ/M (Figure 5)
# ─────────────────────────────────────────────────────────────────────────────

# Relative violation as percentage
df_delta[!, :violation_pct] = ifelse.(df_delta.M .> 0,
    100.0 .* df_delta.delta ./ df_delta.M, 0.0)

# Mean violation per municipality
delta_by_mun = combine(
    groupby(df_delta, :mun),
    :violation_pct => mean => :mean_violation,
    :delta => (d -> count(d .> 1e-6)) => :n_active,
    :delta => length => :n_total
)

# Join with allocation data for size classification
mun_df = innerjoin(q_by_mun, delta_by_mun, on=:mun)
mun_df[!, :class] = size_class.(mun_df.total_allocated)

# Build data for boxplot
df_plot_delta = DataFrame(violation=Float64[], class=Int[])
for row in eachrow(mun_df)
    push!(df_plot_delta, (violation=row.mean_violation, class=1))
    push!(df_plot_delta, (violation=row.mean_violation, class=row.class))
end

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    ylabel="Mean constraint violation (δ/M × 100%)",
    xlabel="Municipalities grouped by number of families",
    xticks=(1:5, labels),
    title="Demographic constraint violation per municipality")

boxplot!(ax2, df_plot_delta.class, df_plot_delta.violation,
         color=:coral, whiskerwidth=0.5, width=0.6)

# Add sample size annotations
for i in 1:5
    n = count(df_plot_delta.class .== i)
    text!(ax2, i, minimum(df_plot_delta.violation) - 0.5,
          text="n=$n", align=(:center, :top), fontsize=10, color=:gray40)
end

save(joinpath(results_dir, "fig5_boxplot_delta.pdf"), fig2)
println("Saved: fig5_boxplot_delta.pdf")

# ─────────────────────────────────────────────────────────────────────────────
# Scatter plot: q vs δ/M (supplementary)
# ─────────────────────────────────────────────────────────────────────────────

fig3 = Figure(size=(600, 500))
ax3 = Axis(fig3[1, 1],
    xlabel="Mean constraint violation (δ/M × 100%)",
    ylabel="Mean allocation quality (q)",
    title="Allocation quality vs. demographic constraint violation")

scatter!(ax3, mun_df.mean_violation, mun_df.q_mean,
         color=(:steelblue, 0.5), markersize=7)

# Linear trend line
x = mun_df.mean_violation
y = mun_df.q_mean
β1 = cov(x, y) / var(x)
β0 = mean(y) - β1 * mean(x)
xr = range(minimum(x), maximum(x), length=100)
lines!(ax3, collect(xr), β0 .+ β1 .* collect(xr),
       color=:red, linewidth=1.5, linestyle=:dash)

r = cor(x, y)
text!(ax3, maximum(x) * 0.6, maximum(y) * 0.97,
      text=@sprintf("r = %.3f", r), fontsize=12, color=:red)

save(joinpath(results_dir, "fig_scatter_q_vs_delta.pdf"), fig3)
println("Saved: fig_scatter_q_vs_delta.pdf")

println("\nAll figures saved to $results_dir")
