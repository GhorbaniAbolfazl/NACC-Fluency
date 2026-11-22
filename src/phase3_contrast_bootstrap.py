"""Bootstrap the registered score contrast by resampling participants within MRI groups."""

import json
from hashlib import sha256
from math import isfinite, sqrt
from pathlib import Path
from platform import python_version
from random import Random
from statistics import NormalDist, fmean, pvariance, stdev

from .build_dataset import write_output
from .phase1_figures import GRAYSCALE, save_pdf
from .phase3_mean_intervals import PLAN_FILE, SAMPLE_NAME
from .phase3_mri_contrast import (
    ESTIMAND,
    GROUP_ORDER,
    REPORT_FILE as CONTRAST_FILE,
    build_report as build_contrast_report,
    load_mri_groups,
)
from .phase3_sampling import percentile_interval

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

REPORT_FILE = Path("reports/phase3/score_difference_bootstrap.json")
FIGURE_DIR = Path("figures/phase3")
FIGURE_NAME = "03_score_difference_bootstrap.pdf"


def resample_difference(groups: dict[str, list[dict]], column: str, replicates: int, seed: int) -> list[float]:
    """Keep each group's size and draw complete rows using one local random stream."""
    higher, lower = (groups[name] for name in GROUP_ORDER)
    if min(len(higher), len(lower), replicates) < 2:
        raise ValueError("A contrast bootstrap requires two participants per group and at least two replicates.")
    if not all(isfinite(row[column]) for members in (higher, lower) for row in members):
        raise ValueError("Every bootstrap participant must have a finite score.")
    generator = Random(seed)
    differences = []
    for _ in range(replicates):
        higher_draw = generator.choices(higher, k=len(higher))
        lower_draw = generator.choices(lower, k=len(lower))
        differences.append(fmean(row[column] for row in higher_draw) - fmean(row[column] for row in lower_draw))
    return differences


def build_report() -> dict:
    plan, groups, flow = load_mri_groups()
    standard = json.loads(CONTRAST_FILE.read_text(encoding="utf-8"))
    if standard != build_contrast_report():
        raise ValueError("The saved Welch contrast differs from the current data or plan. Review and rebuild Step 4.")
    settings = plan["bootstrap"]["mri_contrast"]
    if settings["sample"] != SAMPLE_NAME or settings["statistic"] != ESTIMAND:
        raise ValueError("The registered contrast bootstrap differs from the planned sample or statistic.")
    replicates = plan["bootstrap"]["replicates"]
    score_column = standard["score_column"]
    differences = resample_difference(groups, score_column, replicates, settings["seed"])
    observed = standard["estimate"]
    bootstrap_mean = fmean(differences)
    bootstrap_se = stdev(differences)
    interval = percentile_interval(differences, plan["confidence_level"])
    welch_interval = observed["confidence_interval"]
    exact_empirical_se = sqrt(sum(
        pvariance([row[score_column] for row in groups[name]]) / len(groups[name]) for name in GROUP_ORDER
    ))
    return {
        "stage": "Phase 3, Step 5",
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "analysis_plan_sha256": sha256(PLAN_FILE.read_bytes()).hexdigest(),
        "welch_contrast_sha256": sha256(CONTRAST_FILE.read_bytes()).hexdigest(),
        "sample": SAMPLE_NAME,
        "sample_size": sum(map(len, groups.values())),
        "estimand": plan["estimands"][ESTIMAND],
        "score_column": score_column,
        "mri_grouping": standard["mri_grouping"],
        "group_order": list(GROUP_ORDER),
        "group_sizes": {name: len(groups[name]) for name in GROUP_ORDER},
        "bootstrap": {
            "replicates": replicates,
            "seed": settings["seed"],
            "generator": "Python random.Random, Mersenne Twister, using choices with equal weights",
            "python_version": python_version(),
            "draw_rule": settings["draw_rule"],
            "draw_order": "Within each replicate, draw the at-or-above-cutoff group first, then the below-cutoff group, using one random stream.",
            "row_order": "Original eligible CSV order within each group, identified by the participant-table hash",
            "unit": "One complete participant row, preserving that person's UDS and MRI values together",
            "cutoff": plan["mri_grouping"]["bootstrap_rule"],
            "standard_error": "Sample standard deviation of the bootstrap mean differences, denominator B - 1",
            "exact_empirical_standard_error": "sqrt(v_higher / n_higher + v_lower / n_lower), where each empirical variance v uses denominator n",
            "percentile_interval": plan["bootstrap"]["confidence_interval"],
            "interval_scope": plan["testing_policy"]["interval_scope"],
        },
        "estimate": {
            "unit": plan["estimands"][ESTIMAND]["unit"],
            "observed_mean_difference": observed["mean_difference"],
            "welch_standard_error": observed["standard_error"],
            "exact_empirical_bootstrap_standard_error": exact_empirical_se,
            "bootstrap_mean_difference": bootstrap_mean,
            "estimated_bootstrap_bias": bootstrap_mean - observed["mean_difference"],
            "bootstrap_standard_error": bootstrap_se,
            "bootstrap_to_welch_se_ratio": bootstrap_se / observed["standard_error"],
            "bootstrap_mean_monte_carlo_standard_error": bootstrap_se / sqrt(replicates),
            "percentile_interval": interval,
            "welch_interval": welch_interval,
            "interval_comparison": {
                "lower_percentile_minus_welch": interval["lower"] - welch_interval["lower"],
                "upper_percentile_minus_welch": interval["upper"] - welch_interval["upper"],
                "percentile_to_welch_width_ratio": interval["interval_width"] / welch_interval["interval_width"],
                "welch_includes_zero": welch_interval["lower"] <= 0 <= welch_interval["upper"],
                "percentile_includes_zero": interval["lower"] <= 0 <= interval["upper"],
            },
        },
        "multiplicity": standard["multiplicity"],
        "interpretation": {
            "effect": "Higher MRI group's mean score minus lower MRI group's mean score, in correct animal names in 60 seconds",
            "standard_error": "Estimated spread of the mean difference over repeated independent samples with the same two group sizes",
            "mean_offset": "The empirical-bootstrap expectation of a difference of arithmetic means equals the observed difference. The reported small offset reflects finite-replicate Monte Carlo variation.",
            "comparison": settings["comparison"],
            "confidence_interval": plan["interpretation"]["confidence_interval"],
            "generalization": plan["sampling_model"]["scope"],
            "scope": "Unadjusted exploratory comparison, conditional on the fixed cutoff and observed group sizes; age, sex and education are addressed in Phase 4.",
            "limits": "Bootstrap approximations use the observed within-group empirical distributions. Agreement with Welch does not establish population coverage or practical equivalence.",
            "inference": "This step estimates uncertainty in the contrast. Its resampling distribution is centered near the observed effect and is not a null distribution for calculating a new p-value.",
        },
        "figure": str(FIGURE_DIR / FIGURE_NAME),
        "bootstrap_mean_differences": differences,
    }


def contrast_figure(report: dict):
    estimate = report["estimate"]
    values = report["bootstrap_mean_differences"]
    observed = estimate["observed_mean_difference"]
    welch_se = estimate["welch_standard_error"]
    figure, (density_axis, interval_axis) = plt.subplots(
        2, 1, figsize=(7.8, 6.1), sharex=True, layout="constrained",
        gridspec_kw={"height_ratios": [3.2, 1.1]},
    )
    figure.suptitle("Bootstrap uncertainty in the score difference", fontsize=16, fontweight="bold", color=GRAYSCALE.ink)
    sizes = report["group_sizes"]
    density_axis.set_title(
        f"Normal cognition   |   n = {sizes[GROUP_ORDER[0]]} + {sizes[GROUP_ORDER[1]]}   |   {report['bootstrap']['replicates']:,} within-group resamples",
        loc="left", fontsize=10, color=GRAYSCALE.muted, pad=10,
    )
    for axis in (density_axis, interval_axis):
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["left", "bottom"]].set_color(GRAYSCALE.outline)
        axis.tick_params(length=3, colors=GRAYSCALE.muted)
        axis.set_axisbelow(True)
    density_axis.grid(axis="y", color=GRAYSCALE.grid, linewidth=0.8)
    interval_axis.grid(axis="x", color=GRAYSCALE.grid, linewidth=0.8)
    density_axis.hist(values, bins=40, density=True, color=GRAYSCALE.accent_light, edgecolor=GRAYSCALE.outline, linewidth=0.55, label="Bootstrap differences")
    normal = NormalDist(observed, welch_se)
    lower = min(min(values), observed - 4 * welch_se)
    upper = max(max(values), observed + 4 * welch_se)
    x = [lower + (upper - lower) * i / 400 for i in range(401)]
    density_axis.plot(x, [normal.pdf(v) for v in x], color=GRAYSCALE.ink, linewidth=1.4, linestyle="--", label="Normal reference with Welch SE")
    density_axis.axvline(observed, ymax=0.73, color=GRAYSCALE.ink, linewidth=1.15, label="Observed difference")
    density_axis.axvline(0, ymax=0.73, color=GRAYSCALE.muted, linewidth=1.1, linestyle=":", label="Zero difference")
    density_axis.set(ylabel="Density (per correct response)", ylim=(0, None))
    density_axis.set_ylim(top=density_axis.get_ylim()[1] * 1.35)
    density_axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    density_axis.legend(loc="upper left", frameon=False, fontsize=8)
    annotation = f"Welch SE = {welch_se:.5f}\nBootstrap SE = {estimate['bootstrap_standard_error']:.5f}"
    density_axis.text(0.98, 0.97, annotation, transform=density_axis.transAxes, ha="right", va="top", fontsize=9, color=GRAYSCALE.ink)
    for position, key in ((1, "welch_interval"), (0, "percentile_interval")):
        interval = estimate[key]
        interval_axis.hlines(position, interval["lower"], interval["upper"], color=GRAYSCALE.ink, linewidth=1.8)
        interval_axis.vlines([interval["lower"], interval["upper"]], position - 0.10, position + 0.10, color=GRAYSCALE.ink, linewidth=1.3)
        interval_axis.plot(observed, position, marker="o", color=GRAYSCALE.ink, markersize=4)
    interval_axis.axvline(0, color=GRAYSCALE.muted, linestyle=":", linewidth=1.1)
    interval_axis.set_yticks([1, 0], ["Welch 95% CI", "Bootstrap 95% CI"])
    interval_axis.set(ylim=(-0.55, 1.55), xlabel="Mean score difference (higher - lower, correct responses)")
    interval_axis.xaxis.set_major_locator(MaxNLocator(nbins=7))
    interval_axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    return figure


def main() -> None:
    report = build_report()
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": GRAYSCALE.ink,
        "axes.labelpad": 9,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
    }):
        figure = contrast_figure(report)
        save_pdf(figure, FIGURE_NAME, FIGURE_DIR)
        plt.close(figure)
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()
