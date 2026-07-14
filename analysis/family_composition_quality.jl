# =============================================================================
# family_composition_quality.jl
# Validation of intra-household demographic correlations (Section 4.3)
#
# Generates Figures 1, 2, and 3 of the paper:
#   - Figure 1: Age correlations (scatter) and age differences (histograms)
#               between head and members of orders k = 2, 3, 4
#   - Figure 2: Distribution of gender configurations across municipalities
#   - Figure 3: Row-normalized heatmaps of ethnic composition
#
# =============================================================================

using CSV, DataFrames, StatsBase, Random, Plots, StatsPlots, Statistics
using Plots.PlotMeasures

include(joinpath(@__DIR__, "..", "config.jl"))
include(joinpath(@__DIR__, "..", "src", "utils.jl"))

# ─────────────────────────────────────────────────────────────────────────────
# Main validation and figure generation
# ─────────────────────────────────────────────────────────────────────────────

"""
    validate_family_composition(DFP, ord; output_dir=PATH_RESULTS)

For a given household member order `ord`, generate:
- Per-municipality validation metrics (CSV)
- Age scatter plot (head vs. member)
- Age difference histogram
- Ethnicity heatmap (for ord <= 4)
- Gender configuration boxplot (for ord == 2)

`DFP` is the synthetic person dataset with columns :mun, :ID_Familia,
:Ordem, :genero, :raça, :idade.
"""
function validate_family_composition(DFP::DataFrame, ord::Int;
                                      output_dir::String=PATH_RESULTS)
    mkpath(output_dir)

    # ── Build head–member pairs ───────────────────────────────────────────
    df_heads = DFP[DFP.Ordem .== 1,
                   [:mun, :ID_Familia, :genero, :raça, :idade]]
    df_members = DFP[DFP.Ordem .== ord,
                     [:mun, :ID_Familia, :genero, :raça, :idade]]

    df_pairs = innerjoin(df_heads, df_members,
                         on=[:mun, :ID_Familia], makeunique=true)
    # After join: :idade = head age, :idade_1 = member age
    #             :genero/:genero_1, :raça/:raça_1

    println("Order $ord: $(nrow(df_pairs)) head–member pairs")

    # ── Per-municipality metrics ──────────────────────────────────────────
    metrics = combine(groupby(df_pairs, :mun)) do sub
        (Households       = nrow(sub),
         Age_Correlation  = round(cor(sub.idade, sub.idade_1), digits=3),
         Mean_Age_Diff    = round(mean(sub.idade .- sub.idade_1), digits=2),
         Same_Gender_Rate = round(mean(sub.genero .== sub.genero_1), digits=4),
         Same_Ethnic_Rate = round(mean(sub.raça .== sub.raça_1), digits=3))
    end
    CSV.write(joinpath(output_dir, "validation_order$(ord).csv"), metrics)

    # ── Sample for scatter and histogram (cap at 100k for performance) ────
    n_sample = min(100_000, nrow(df_pairs))
    sample_idx = sample(1:nrow(df_pairs), n_sample, replace=false)
    df_sample = df_pairs[sample_idx, :]

    # ── Figure 1a/c/e: Age difference histogram ──────────────────────────
    df_sample.age_diff = df_sample.idade .- df_sample.idade_1

    histogram(df_sample.age_diff,
        title="Age Difference (Head - Member $ord)",
        xlabel="Years of Difference",
        ylabel="Frequency",
        legend=false,
        xticks=-100:10:100)
    savefig(joinpath(output_dir, "fig1_histogram_age_diff_order$(ord).pdf"))

    # ── Figure 1b/d/f: Age scatter plot ──────────────────────────────────
    scatter(df_sample.idade, df_sample.idade_1,
        alpha=0.15, markersize=2, legend=false,
        title="",
        xlabel="Age of head",
        ylabel="Age of Member $ord",
        xticks=0:10:100, dpi=300)

    # Reference lines
    plot!([0, 100], [0, 100], color=:red, linewidth=2, linestyle=:dash)  # y = x (same age)
    plot!([25, 100], [0, 75], color=:blue, linewidth=1)   # head is parent (~25 yr older)
    plot!([0, 75], [25, 100], color=:green, linewidth=1)  # head is child (~25 yr younger)
    savefig(joinpath(output_dir, "fig1_scatter_age_order$(ord).png"))

    # ── Figure 3: Ethnicity heatmap (row-normalized, orders 2–4) ─────────
    if ord <= 4
        ethnic_codes = sort(intersect(unique(df_pairs.raça), 1:5))
        ethnic_labels = Dict(1 => "White", 2 => "Black", 3 => "Asian",
                             4 => "Mixed", 5 => "Indigenous")
        labels_vec = [ethnic_labels[c] for c in ethnic_codes]

        if !isempty(ethnic_codes)
            # Confusion matrix: rows = head ethnicity, cols = member ethnicity
            mat = [nrow(df_pairs[(df_pairs.raça .== r1) .&
                                  (df_pairs.raça_1 .== r2), :])
                   for r1 in ethnic_codes, r2 in ethnic_codes]

            # Row-normalize (each row sums to 1)
            row_sums = sum(mat, dims=2)
            row_sums[row_sums .== 0] .= 1.0
            proportions = mat ./ row_sums

            h = heatmap(labels_vec, labels_vec, proportions',
                title=" ",
                xlabel="Race of Head",
                ylabel="Race of Member $ord",
                color=:viridis,
                clim=(0, 1),
                aspect_ratio=:equal,
                grid=true, framestyle=:box,
                size=(700, 700),
                left_margin=15mm, bottom_margin=15mm)
            savefig(h, joinpath(output_dir,
                                "fig3_heatmap_ethnicity_order$(ord).png"))
        end
    end

    # ── Figure 2: Gender configuration boxplot (order 2 only) ─────────────
    if ord == 2
        # Compute gender pairing proportions per municipality
        df_pairs.pair_type = Symbol.(df_pairs.genero, "_", df_pairs.genero_1)

        counts = combine(groupby(df_pairs, [:mun, :pair_type]), nrow => :n)
        props = transform(groupby(counts, :mun),
                          :n => (x -> x ./ sum(x)) => :prop)

        df_wide = unstack(props, :mun, :pair_type, :prop)
        df_wide = coalesce(df_wide, 0.0)

        # Ensure all four columns exist
        for col in [:F_M, :F_F, :M_F, :M_M]
            if !(string(col) in names(df_wide))
                df_wide[!, col] = zeros(Float64, nrow(df_wide))
            else
                df_wide[!, col] = Float64.(coalesce.(df_wide[:, col], 0.0))
            end
        end

        p = boxplot(
            ["Female-Male" "Female-Female" "Male-Female" "Male-Male"],
            [df_wide.F_M, df_wide.F_F, df_wide.M_F, df_wide.M_M],
            title="Distribution of Gender Configurations by Municipality",
            ylabel="Proportional Frequency",
            xlabel="Configuration (Head - Member 2)",
            legend=false,
            size=(1200, 800),
            titlefontsize=16, guidefontsize=14, tickfontsize=12,
            palette=:Paired, lw=1.5, fillalpha=0.9,
            left_margin=15mm, bottom_margin=10mm,
            right_margin=5mm, top_margin=5mm)

        savefig(joinpath(output_dir, "fig2_boxplot_gender_order2.pdf"))
    end

    println("Order $ord: figures saved to $output_dir")
end

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────
# Usage:
#   DFP = CSV.read(joinpath(PATH_STAGE3, "dataset_pessoas_com_Ordem_final.csv"), DataFrame)
#   for ord in [2, 3, 4]
#       validate_family_composition(DFP, ord)
#   end
