"""Draw grayscale residual and influence panels from private in-memory fits."""

import numpy as np


MODEL_LABELS = {"uds_only": "UDS-only", "mri_only": "MRI-only", "integrated": "UDS + MRI"}


def style(axis, xlabel, ylabel):
    from .phase1_figures import GRAYSCALE
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color(GRAYSCALE.outline)
    axis.tick_params(colors=GRAYSCALE.muted, length=3)
    axis.set(xlabel=xlabel, ylabel=ylabel)
    axis.grid(color=GRAYSCALE.grid, linewidth=.6)
    axis.set_axisbelow(True)


def scatter(axis, x, y):
    from .phase1_figures import GRAYSCALE
    axis.scatter(x, y, s=9, color=GRAYSCALE.accent, alpha=.45, edgecolors="none", rasterized=False)


def bin_line(axis, bins, key):
    from .phase1_figures import GRAYSCALE
    axis.plot([v["x_mean"] for v in bins], [v[key] for v in bins], "o-", color=GRAYSCALE.ink,
              linewidth=1.35, markersize=3, label="Bin means")


def residual_figure(report, context):
    from .phase1_figures import GRAYSCALE, plt
    figure, axes = plt.subplots(3, 3, figsize=(12.2, 10.1), layout="constrained")
    figure.suptitle("Residual diagnostics for the three registered models\n"
                    "n = 842 per model | Black lines: descriptive bin means; Q-Q: quartile reference",
                    fontsize=14, fontweight="bold", color=GRAYSCALE.ink)
    for row, name in enumerate(("uds_only", "mri_only", "integrated")):
        fit, a, model = context["fits"][name], context["arrays"][name], report["models"][name]
        axis = axes[row, 0]
        scatter(axis, fit.fitted_values, fit.residuals)
        axis.axhline(0, color=GRAYSCALE.outline, linestyle=":", linewidth=1)
        bin_line(axis, model["fitted_value_bins"], "residual_mean")
        style(axis, "Fitted score", "Residual (correct responses)")
        axis.set_title(MODEL_LABELS[name] + ": residual vs fitted", loc="left", fontsize=10, fontweight="bold")
        axis = axes[row, 1]
        scatter(axis, a["qq_theoretical"], a["qq_observed"])
        ends = np.asarray([min(a["qq_theoretical"]), max(a["qq_theoretical"])])
        axis.plot(ends, a["qq_line_intercept"] + a["qq_line_slope"] * ends, color=GRAYSCALE.ink, linewidth=1.3)
        style(axis, "Standard-normal quantile", "Internally studentized residual")
        axis.set_title(MODEL_LABELS[name] + ": normal Q-Q", loc="left", fontsize=10, fontweight="bold")
        axis = axes[row, 2]
        scatter(axis, fit.fitted_values, np.sqrt(np.abs(a["studentized"])))
        bin_line(axis, model["fitted_value_bins"], "mean_sqrt_abs_studentized")
        style(axis, "Fitted score", "sqrt(abs(studentized residual))")
        axis.set_title(MODEL_LABELS[name] + ": scale-location", loc="left", fontsize=10, fontweight="bold")
    return figure


def predictor_figure(report, context):
    from .phase1_figures import GRAYSCALE, plt
    from matplotlib.ticker import MaxNLocator
    figure, axes = plt.subplots(2, 2, figsize=(10.7, 7.2), layout="constrained")
    figure.suptitle("Integrated model: residuals by registered predictor\n"
                    "n = 842 | All observations retained; bin means are descriptive", fontsize=14, fontweight="bold", color=GRAYSCALE.ink)
    fit, design = context["fits"]["integrated"], context["designs"]["integrated"]
    for axis, (term, title, xlabel, reference, scale) in zip(axes.flat, (
        ("thickness_0_1mm", "Left inferior parietal thickness", "Thickness (mm)", 2.344, .1),
        ("age_c", "Age", "Age (years)", 68, 1),
        ("education_c", "Education", "Education (years)", 16, 1),
    )):
        j = design.columns.index(term)
        raw = [row[j] * scale + reference for row in design.matrix]
        scatter(axis, raw, fit.residuals)
        axis.axhline(0, color=GRAYSCALE.outline, linestyle=":", linewidth=1)
        bin_line(axis, report["integrated_predictor_bins"][term], "residual_mean")
        style(axis, xlabel, "Residual (correct responses)")
        axis.set_title(title, loc="left", fontsize=11, fontweight="bold")
        if scale == 1:
            axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    axis = axes[1, 1]
    j = design.columns.index("female")
    groups = [[e for row, e in zip(design.matrix, fit.residuals) if row[j] == code] for code in (0, 1)]
    axis.boxplot(groups, positions=[0, 1], patch_artist=True, widths=.45,
                         boxprops={"facecolor": GRAYSCALE.accent_light, "edgecolor": GRAYSCALE.outline},
                         medianprops={"color": GRAYSCALE.ink, "linewidth": 1.5},
                         whiskerprops={"color": GRAYSCALE.outline}, capprops={"color": GRAYSCALE.outline},
                         flierprops={"marker": "o", "markersize": 3, "markerfacecolor": "none", "markeredgecolor": GRAYSCALE.outline})
    axis.set_xticks([0, 1], [f"Recorded male\nn = {len(groups[0])}", f"Recorded female\nn = {len(groups[1])}"])
    axis.axhline(0, color=GRAYSCALE.outline, linestyle=":", linewidth=1)
    style(axis, "Recorded sex (male reference)", "Residual (correct responses)")
    axis.set_title("Residual spread by recorded sex", loc="left", fontsize=11, fontweight="bold")
    return figure


def influence_figure(report, context):
    from .phase1_figures import GRAYSCALE, plt
    from matplotlib.ticker import MaxNLocator
    figure, axes = plt.subplots(1, 3, figsize=(12.4, 4.15), layout="constrained")
    fit, a = context["fits"]["integrated"], context["arrays"]["integrated"]
    count = report["models"]["integrated"]["review_flags"]["union_count"]
    figure.suptitle(f"Integrated model: influence review\nFull primary sample n = 842 | {count} flagged observations; no automatic deletion",
                    fontsize=14, fontweight="bold", color=GRAYSCALE.ink)
    axis = axes[0]
    scatter(axis, fit.leverages, a["studentized"])
    for threshold in (-a["limits"]["studentized"], a["limits"]["studentized"]):
        axis.axhline(threshold, color=GRAYSCALE.ink, linestyle="--", linewidth=1)
    axis.axvline(a["limits"]["leverage"], color=GRAYSCALE.outline, linestyle=":", linewidth=1.2)
    style(axis, "Leverage", "Internally studentized residual")
    axis.set_title("Residual / leverage flags", loc="left", fontsize=10, fontweight="bold")
    axis = axes[1]
    scatter(axis, fit.leverages, a["cooks_distance"])
    axis.axhline(a["limits"]["cooks_distance"], color=GRAYSCALE.ink, linestyle="--", linewidth=1)
    axis.axvline(a["limits"]["leverage"], color=GRAYSCALE.outline, linestyle=":", linewidth=1.2)
    style(axis, "Leverage", "Cook's distance")
    axis.set_title("Cook / leverage flags", loc="left", fontsize=10, fontweight="bold")
    axis = axes[2]
    estimates = [q["estimate"] for q in context["omissions"]]
    if estimates:
        axis.hist(estimates, bins=min(12, max(3, int(np.sqrt(len(estimates))))), color=GRAYSCALE.accent_light,
                  edgecolor=GRAYSCALE.outline, linewidth=.8)
    else:
        axis.text(.5, .5, "No flagged observations", transform=axis.transAxes, ha="center")
    full_slope = report["integrated_influence_followup"]["primary_full_sample"]["thickness"]["estimate"]
    axis.axvline(full_slope, color=GRAYSCALE.ink, linestyle="--", linewidth=1.3, label="Full-sample estimate")
    axis.legend(frameon=False, fontsize=8, loc="upper left")
    style(axis, "Thickness coefficient per 0.1 mm", "Diagnostic omission fits")
    axis.set_title(f"One flagged omission per fit (n = 841)", loc="left", fontsize=10, fontweight="bold")
    axis.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
    return figure


def build_figures(report, context):
    return residual_figure(report, context), predictor_figure(report, context), influence_figure(report, context)
