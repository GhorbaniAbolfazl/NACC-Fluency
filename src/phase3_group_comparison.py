"""Compare MRI thickness across UDS groups and complete the planned two-test family."""

import json
from collections import Counter
from hashlib import sha256
from importlib.metadata import version
from math import fsum, isfinite, sqrt
from pathlib import Path
from statistics import NormalDist, fmean, pvariance, variance

from scipy.stats import f

from .build_dataset import LEFT_THICKNESS, STATUS, write_output
from .data_audit import CATEGORY_LABELS
from .descriptive_statistics import numeric_summary
from .phase1_figures import GRAYSCALE, GROUP_LABELS, make_figure, save_pdf
from .phase2_samples import load_phase2_samples
from .phase3_mean_intervals import PLAN_FILE, load_registered_plan, mean_interval
from .phase3_mri_contrast import REPORT_FILE as CONTRAST_FILE, build_report as build_contrast_report

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

REPORT_FILE = Path("reports/phase3/thickness_by_status.json")
FAMILY_FILE = Path("reports/phase3/hypothesis_tests.json")
FIGURE_DIR = Path("figures/phase3")
MEANS_FIGURE = "04_thickness_group_means.pdf"
DIAGNOSTICS_FIGURE = "05_thickness_group_qq.pdf"
SAMPLE_NAME = "status_and_mri"
ESTIMAND = "thickness_by_status"


def encode_report(report: dict) -> bytes:
    return (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def load_status_groups() -> tuple[dict, dict[int, list[dict]], dict]:
    plan = load_registered_plan()
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    codes = plan["estimands"][ESTIMAND]["group_codes"]
    observed = Counter(CATEGORY_LABELS[STATUS][row[STATUS]] for row in rows)
    expected = plan["samples"][SAMPLE_NAME]
    if (
        len(rows) != expected["included"]
        or any(observed[label] != count for label, count in expected["included_by_cognitive_status"].items())
        or set(row[STATUS] for row in rows) != set(codes)
        or flow["participants_csv_sha256"] != plan["inputs"]["participants_csv_sha256"]
    ):
        raise ValueError("The cognitive-status comparison sample differs from the registered Phase 3 plan.")
    return plan, {code: [row for row in rows if row[STATUS] == code] for code in codes}, flow


def welch_anova(samples: list[list[float]]) -> dict:
    """Use Welch's weighted F statistic and its approximate denominator degrees of freedom."""
    if len(samples) < 2 or any(len(values) < 2 for values in samples):
        raise ValueError("Welch ANOVA requires at least two groups with two observations each.")
    if not all(isfinite(value) for values in samples for value in values):
        raise ValueError("Welch ANOVA requires finite observations.")
    variances = [variance(values) for values in samples]
    if any(value <= 0 for value in variances):
        raise ValueError("Welch ANOVA requires a positive sample variance in each group.")
    k = len(samples)
    weights = [len(values) / v for values, v in zip(samples, variances)]
    total_weight = fsum(weights)
    means = [fmean(values) for values in samples]
    weighted_mean = fsum(w * mean for w, mean in zip(weights, means)) / total_weight
    correction_sum = fsum((1 - w / total_weight) ** 2 / (len(values) - 1) for w, values in zip(weights, samples))
    numerator = fsum(w * (mean - weighted_mean) ** 2 for w, mean in zip(weights, means)) / (k - 1)
    correction = 1 + 2 * (k - 2) * correction_sum / (k * k - 1)
    statistic = numerator / correction
    df_numerator = k - 1
    df_denominator = (k * k - 1) / (3 * correction_sum)
    return {
        "f_statistic": statistic,
        "df_numerator": df_numerator,
        "df_denominator": df_denominator,
        "p_value_raw": float(f.sf(statistic, df_numerator, df_denominator)),
        "inverse_variance_weighted_mean": weighted_mean,
        "precision_weights": weights,
        "welch_correction": correction,
    }


def descriptive_eta_squared(samples: list[list[float]]) -> dict:
    """Use unweighted participant sums of squares for a descriptive sample effect."""
    all_values = [value for values in samples for value in values]
    grand_mean = fmean(all_values)
    means = [fmean(values) for values in samples]
    between = fsum(len(values) * (mean - grand_mean) ** 2 for values, mean in zip(samples, means))
    within = fsum(fsum((value - mean) ** 2 for value in values) for values, mean in zip(samples, means))
    total = fsum((value - grand_mean) ** 2 for value in all_values)
    if total <= 0:
        raise ValueError("Descriptive eta-squared requires nonzero total variation.")
    return {
        "grand_mean_mm": grand_mean,
        "sum_of_squares_unit": "mm^2",
        "between_group_sum_of_squares": between,
        "within_group_sum_of_squares": within,
        "total_sum_of_squares": total,
        "eta_squared": between / total,
        "definition": "Between-group sum of squares divided by total sum of squares in the observed sample",
        "interpretation": "A descriptive share of observed squared thickness variation associated with differences among group means. It depends on sample composition and is separate from Welch's precision-weighted F statistic.",
    }


def holm_adjust(raw_p_values: dict[str, float]) -> dict[str, float]:
    if not raw_p_values or any(not isfinite(p) or not 0 <= p <= 1 for p in raw_p_values.values()):
        raise ValueError("Holm adjustment requires finite p-values between zero and one.")
    adjusted = {}
    running_maximum = 0.0
    for rank, name in enumerate(sorted(raw_p_values, key=raw_p_values.get)):
        running_maximum = max(running_maximum, (len(raw_p_values) - rank) * raw_p_values[name])
        adjusted[name] = min(1.0, running_maximum)
    return {name: adjusted[name] for name in raw_p_values}


def build_report() -> dict:
    plan, groups, flow = load_status_groups()
    samples = [[row[LEFT_THICKNESS] for row in members] for members in groups.values()]
    result = welch_anova(samples)
    group_results = {}
    for (code, members), values in zip(groups.items(), samples):
        summary = numeric_summary(members, LEFT_THICKNESS)
        empirical_sd = sqrt(pvariance(values))
        skewness = fsum(((value - summary["mean"]) / empirical_sd) ** 3 for value in values) / len(values)
        group_results[str(code)] = {
            "label": CATEGORY_LABELS[STATUS][code],
            "summary": summary,
            "mean_interval": mean_interval(values, plan["confidence_level"]),
            "empirical_skewness": skewness,
        }
    variances = [variance(values) for values in samples]
    return {
        "stage": "Phase 3, Step 6",
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "analysis_plan_sha256": sha256(PLAN_FILE.read_bytes()).hexdigest(),
        "analysis_status": plan["analysis_status"],
        "sample": SAMPLE_NAME,
        "sample_size": sum(map(len, groups.values())),
        "outcome": {"column": LEFT_THICKNESS, "unit": "mm"},
        "grouping_column": STATUS,
        "group_order": list(groups),
        "groups": group_results,
        "hypothesis_test": plan["hypothesis_tests"][ESTIMAND],
        "method": {
            "name": "Welch one-way ANOVA",
            "formula": "w_i=n_i/s_i^2; C=sum((1-w_i/sum(w))^2/(n_i-1)); F=[sum(w_i*(mean_i-weighted_mean)^2)/(k-1)]/[1+2*(k-2)*C/(k^2-1)]; df2=(k^2-1)/(3*C)",
            "p_value": "Upper-tail F probability, scipy.stats.f.sf(F, k - 1, df2)",
            "scipy_version": version("scipy"),
            "mean_intervals": plan["estimands"][ESTIMAND]["standard_interval"],
            "interval_scope": plan["testing_policy"]["interval_scope"],
            "assumptions": "Independent participants and approximate normal-theory calibration within groups. Unequal population variances and unequal group sizes are allowed. The seven-person group has limited precision and weak protection from a large-sample approximation.",
        },
        "welch_anova": result,
        "descriptive_effect": descriptive_eta_squared(samples),
        "diagnostics": {
            "largest_to_smallest_group_size_ratio": max(map(len, samples)) / min(map(len, samples)),
            "largest_to_smallest_sample_variance_ratio": max(variances) / min(variances),
            "review": "Review group counts, ranges, quartiles, SDs, skewness, IQR flags and the group-specific normal Q-Q plots. Keep all measurements outside review fences.",
            "small_group": "The impairment-without-MCI group has seven participants. Its mean interval and Welch approximation are sensitive to departures from normality and individual extreme values.",
            "normality": "The empirical skewness and Q-Q plots describe tail asymmetry. Welch does not remove sensitivity to non-normal observations in small groups; the test is a model-based exploratory comparison.",
        },
        "interpretation": {
            "omnibus": "Rejection supports a difference in at least one population mean under the model assumptions. It does not identify a particular pair of groups.",
            "effect": "Report group means and intervals in millimeters alongside descriptive sample eta-squared.",
            "confidence_interval": plan["interpretation"]["confidence_interval"],
            "generalization": plan["sampling_model"]["scope"],
            "scope": "An unadjusted exploratory comparison; age, sex and education are addressed in Phase 4.",
        },
        "completed_test_family": str(FAMILY_FILE),
        "figures": [str(FIGURE_DIR / name) for name in (MEANS_FIGURE, DIAGNOSTICS_FIGURE)],
    }


def build_family_report(comparison: dict) -> dict:
    plan = load_registered_plan()
    contrast = json.loads(CONTRAST_FILE.read_text(encoding="utf-8"))
    if contrast != build_contrast_report():
        raise ValueError("The saved score contrast differs from its current inputs. Review Step 4 before adjustment.")
    if (
        comparison["participants_csv_sha256"] != contrast["participants_csv_sha256"]
        or comparison["analysis_plan_sha256"] != contrast["analysis_plan_sha256"]
    ):
        raise ValueError("The planned tests have different input data or analysis plans.")
    names = plan["testing_policy"]["family"]
    raw = {"score_difference_by_mri": contrast["estimate"]["p_value_raw"], ESTIMAND: comparison["welch_anova"]["p_value_raw"]}
    if set(names) != set(raw):
        raise ValueError("The computed tests differ from the registered Holm family.")
    adjusted = holm_adjust(raw)
    alpha = plan["testing_policy"]["familywise_alpha"]
    c = contrast["estimate"]
    a = comparison["welch_anova"]
    statistics = {"score_difference_by_mri": (c["t_statistic"], [c["degrees_of_freedom"]]), ESTIMAND: (a["f_statistic"], [a["df_numerator"], a["df_denominator"]])}
    return {
        "stage": "Phase 3, Step 6, completed planned test family",
        "participants_csv_sha256": comparison["participants_csv_sha256"],
        "analysis_plan_sha256": comparison["analysis_plan_sha256"],
        "inputs": {
            "score_contrast": str(CONTRAST_FILE),
            "score_contrast_sha256": sha256(CONTRAST_FILE.read_bytes()).hexdigest(),
            "thickness_comparison": str(REPORT_FILE),
            "thickness_comparison_sha256": sha256(encode_report(comparison)).hexdigest(),
        },
        "policy": plan["testing_policy"],
        "status": "Complete. This is the final multiplicity table for both planned tests; earlier step outputs retain their stage-specific pending status.",
        "adjustment": "Sort raw p-values ascending; multiply by the number of remaining tests; take cumulative maxima and cap at one.",
        "tests": {
            name: {
                "method": plan["hypothesis_tests"][name]["method"],
                "sample": plan["hypothesis_tests"][name]["sample"],
                "participants": contrast["sample_size"] if name == "score_difference_by_mri" else comparison["sample_size"],
                "null": plan["hypothesis_tests"][name]["null"],
                "alternative": plan["hypothesis_tests"][name]["alternative"],
                "statistic": statistics[name][0],
                "degrees_of_freedom": statistics[name][1],
                "p_value_raw": raw[name],
                "p_value_holm": adjusted[name],
                "reject_holm": adjusted[name] <= alpha,
            }
            for name in names
        },
        "interpretation": {
            "familywise_error": "Holm controls the familywise error rate under valid raw p-values, including dependent tests.",
            "confidence_intervals": "Reported 95% intervals are individual intervals; the Holm adjustment applies to the two test p-values.",
            "p_value": plan["interpretation"]["p_value"],
            "omnibus": comparison["interpretation"]["omnibus"],
        },
    }


def means_figure(report: dict):
    figure, axis = make_figure("Mean left thickness by cognitive status", "All matched participants   |   N = 1,237   |   Individual 95% t intervals", grid_axis="x")
    labels = []
    for position, code in enumerate(report["group_order"]):
        group = report["groups"][str(code)]
        interval = group["mean_interval"]
        axis.hlines(position, interval["lower"], interval["upper"], color=GRAYSCALE.outline, linewidth=1.7)
        axis.vlines([interval["lower"], interval["upper"]], position - 0.08, position + 0.08, color=GRAYSCALE.outline, linewidth=1.3)
        axis.plot(interval["mean"], position, marker="o", color=GRAYSCALE.ink, markersize=5)
        labels.append(f"{GROUP_LABELS[code]}\n(n = {interval['participants']:,})")
    axis.set_yticks(range(len(labels)), labels)
    axis.set(ylim=(len(labels) - 0.5, -0.5), xlabel="Left inferior parietal thickness (mm)")
    axis.xaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    axis.margins(x=0.08)
    return figure


def qq_figure(report: dict, groups: dict[int, list[dict]]):
    figure, axes = plt.subplots(2, 2, figsize=(8.6, 7.1), layout="constrained")
    figure.suptitle("Thickness distributions within cognitive groups", fontsize=16, fontweight="bold", color=GRAYSCALE.ink)
    for axis, code in zip(axes.flat, report["group_order"]):
        group = report["groups"][str(code)]
        summary = group["summary"]
        values = sorted(row[LEFT_THICKNESS] for row in groups[code])
        normal = NormalDist(summary["mean"], summary["sample_sd"])
        theoretical = [normal.inv_cdf((i + 0.5) / len(values)) for i in range(len(values))]
        lower, upper = min(min(values), min(theoretical)), max(max(values), max(theoretical))
        padding = 0.04 * (upper - lower)
        limits = (lower - padding, upper + padding)
        axis.plot(limits, limits, color=GRAYSCALE.outline, linestyle="--", linewidth=1)
        axis.scatter(theoretical, values, s=13, color=GRAYSCALE.ink, linewidths=0, alpha=0.8)
        axis.set(title=f"{GROUP_LABELS[code]} (n = {len(values):,})", xlim=limits, ylim=limits, xlabel="Fitted normal quantiles (mm)", ylabel="Observed thickness (mm)")
        axis.set_aspect("equal", adjustable="box")
        axis.text(0.04, 0.96, f"Skewness = {group['empirical_skewness']:.2f}", transform=axis.transAxes, va="top", fontsize=9, color=GRAYSCALE.muted)
        axis.spines[["top", "right"]].set_visible(False)
        axis.spines[["left", "bottom"]].set_color(GRAYSCALE.outline)
        axis.tick_params(colors=GRAYSCALE.muted)
        axis.grid(color=GRAYSCALE.grid, linewidth=0.6)
        axis.set_axisbelow(True)
        for coordinate in (axis.xaxis, axis.yaxis):
            coordinate.set_major_locator(MaxNLocator(nbins=4))
            coordinate.set_major_formatter(StrMethodFormatter("{x:.2f}"))
    return figure


def main() -> None:
    comparison = build_report()
    family = build_family_report(comparison)
    _, groups, _ = load_status_groups()
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "axes.labelcolor": GRAYSCALE.ink, "axes.labelpad": 9, "pdf.fonttype": 42, "figure.facecolor": "white"}):
        for name, figure in ((MEANS_FIGURE, means_figure(comparison)), (DIAGNOSTICS_FIGURE, qq_figure(comparison, groups))):
            save_pdf(figure, name, FIGURE_DIR)
            plt.close(figure)
    write_output(REPORT_FILE, encode_report(comparison))
    write_output(FAMILY_FILE, encode_report(family))


if __name__ == "__main__":
    main()
