"""Resample complete participants to describe uncertainty in both sample means."""

import json
from hashlib import sha256
from math import floor, isfinite, sqrt
from pathlib import Path
from platform import python_version
from random import Random
from statistics import NormalDist, fmean, stdev

from .build_dataset import write_output
from .phase1_figures import GRAYSCALE, make_figure, save_pdf
from .phase2_samples import load_phase2_samples
from .phase3_mean_intervals import (
    ESTIMANDS,
    PLAN_FILE,
    REPORT_FILE as MEANS_FILE,
    SAMPLE_NAME,
    load_registered_plan,
    mean_interval,
)

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

REPORT_FILE = Path("reports/phase3/sampling_distributions.json")
FIGURE_DIR = Path("figures/phase3")
FIGURES = {
    "mean_score": "01_mean_score_sampling.pdf",
    "mean_thickness": "02_mean_thickness_sampling.pdf",
}


def resample_means(rows: list[dict], columns: tuple[str, ...], replicates: int, seed: int) -> dict:
    """Use the same sampled participant rows for every measure in a replicate."""
    if len(rows) < 2 or replicates < 2:
        raise ValueError("Sampling distributions require at least two participants and two replicates.")
    if not columns or any(not isfinite(row[column]) for row in rows for column in columns):
        raise ValueError("Every bootstrap measure must have finite values for every participant.")

    generator = Random(seed)
    means = {column: [] for column in columns}
    for _ in range(replicates):
        draw = generator.choices(rows, k=len(rows))
        for column in columns:
            means[column].append(fmean(row[column] for row in draw))
    return means


def percentile_interval(values: list[float], confidence_level: float) -> dict:
    """Interpolate at rank (B - 1) * p, including both observed endpoints."""
    if len(values) < 2 or not 0 < confidence_level < 1 or not all(isfinite(v) for v in values):
        raise ValueError("A percentile interval requires finite replicates and a confidence level between zero and one.")
    ordered = sorted(values)

    def interpolate(probability: float) -> float:
        rank = (len(ordered) - 1) * probability
        lower_index = floor(rank)
        fraction = rank - lower_index
        upper_index = min(lower_index + 1, len(ordered) - 1)
        return ordered[lower_index] + fraction * (ordered[upper_index] - ordered[lower_index])

    tail = (1 - confidence_level) / 2
    lower, upper = interpolate(tail), interpolate(1 - tail)
    return {
        "confidence_level": confidence_level,
        "lower": lower,
        "upper": upper,
        "interval_width": upper - lower,
    }


def build_report() -> dict:
    plan = load_registered_plan()
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    plan_hash = sha256(PLAN_FILE.read_bytes()).hexdigest()
    standard = json.loads(MEANS_FILE.read_text(encoding="utf-8"))
    if (
        len(rows) != plan["samples"][SAMPLE_NAME]["included"]
        or flow["participants_csv_sha256"] != plan["inputs"]["participants_csv_sha256"]
        or standard["participants_csv_sha256"] != flow["participants_csv_sha256"]
        or standard["analysis_plan_sha256"] != plan_hash
        or standard["sample"] != SAMPLE_NAME
        or standard["sample_size"] != len(rows)
    ):
        raise ValueError("The saved mean intervals or sampling inputs differ from the registered Phase 3 plan.")
    for name in ESTIMANDS:
        specification = plan["estimands"][name]
        expected = {
            "column": specification["column"],
            "parameter": specification["parameter"],
            "unit": specification["unit"],
            **mean_interval([row[specification["column"]] for row in rows], plan["confidence_level"]),
        }
        if standard["estimates"][name] != expected:
            raise ValueError("The saved mean estimates differ from the current data. Rebuild the mean intervals.")

    settings = plan["bootstrap"]["whole_sample"]
    if settings["sample"] != SAMPLE_NAME or settings["statistics"] != list(ESTIMANDS):
        raise ValueError("The registered whole-sample bootstrap does not match the two planned means.")
    columns = tuple(plan["estimands"][name]["column"] for name in ESTIMANDS)
    replicates = plan["bootstrap"]["replicates"]
    draws = resample_means(rows, columns, replicates, settings["seed"])
    results = {}
    for name in ESTIMANDS:
        estimate = standard["estimates"][name]
        values = draws[estimate["column"]]
        bootstrap_mean = fmean(values)
        bootstrap_se = stdev(values)
        analytic_se = estimate["standard_error"]
        interval = percentile_interval(values, plan["confidence_level"])
        results[name] = {
            "column": estimate["column"],
            "unit": estimate["unit"],
            "observed_mean": estimate["mean"],
            "sample_sd": estimate["sample_sd"],
            "analytic_standard_error": analytic_se,
            "exact_empirical_bootstrap_standard_error": analytic_se * sqrt((len(rows) - 1) / len(rows)),
            "bootstrap_mean": bootstrap_mean,
            "estimated_bootstrap_bias": bootstrap_mean - estimate["mean"],
            "bootstrap_standard_error": bootstrap_se,
            "bootstrap_to_analytic_se_ratio": bootstrap_se / analytic_se if analytic_se else None,
            "bootstrap_mean_monte_carlo_standard_error": bootstrap_se / sqrt(replicates),
            "percentile_interval": interval,
            "t_interval": {key: estimate[key] for key in ("confidence_level", "lower", "upper", "interval_width")},
            "interval_comparison": {
                "lower_percentile_minus_t": interval["lower"] - estimate["lower"],
                "upper_percentile_minus_t": interval["upper"] - estimate["upper"],
                "percentile_to_t_width_ratio": interval["interval_width"] / estimate["interval_width"] if estimate["interval_width"] else None,
            },
            "figure": str(FIGURE_DIR / FIGURES[name]),
        }
    return {
        "stage": "Phase 3, Step 3",
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "analysis_plan_sha256": plan_hash,
        "mean_intervals_sha256": sha256(MEANS_FILE.read_bytes()).hexdigest(),
        "sample": SAMPLE_NAME,
        "sample_size": len(rows),
        "bootstrap": {
            "replicates": replicates,
            "participants_per_replicate": len(rows),
            "seed": settings["seed"],
            "generator": "Python random.Random, Mersenne Twister, using choices with equal weights",
            "python_version": python_version(),
            "row_order": "Eligible participant rows in the processed CSV order, identified by the input hash",
            "unit": "One complete participant row; all measures use the same sampled rows in each replicate",
            "draw_rule": settings["draw_rule"],
            "standard_error": "Sample standard deviation of the bootstrap means, denominator B - 1",
            "exact_empirical_standard_error": "Analytic SE times sqrt((n - 1) / n), for the empirical distribution with replacement",
            "percentile_interval": plan["bootstrap"]["confidence_interval"],
            "finite_replicates": "SE and percentile estimates have Monte Carlo variation from the finite number of replicates.",
            "replicate_alignment": "Element i in both bootstrap_means arrays comes from the same participant draw.",
        },
        "interpretation": {
            "sample_sd": "Spread of individual participants' measurements in the original units",
            "mean_offset": "For an arithmetic mean, the empirical-bootstrap expectation equals the observed mean. The reported bootstrap-mean offset reflects finite-replicate Monte Carlo variation.",
            "standard_error": "Estimated spread of the mean over repeated samples of the same size",
            "sampling_approximation": "Bootstrap distributions approximate the sampling distributions using the observed empirical distribution.",
            "normal_reference": "Reference curves are centered at the observed mean with SD equal to the analytic SE; they do not assume normal individual measurements.",
            "confidence_intervals": plan["interpretation"]["confidence_interval"],
            "interval_scope": plan["testing_policy"]["interval_scope"],
            "generalization": plan["sampling_model"]["scope"],
            "limits": "Agreement between bootstrap and t intervals in this sample does not establish their population coverage or remove sampling-frame limitations.",
        },
        "estimates": results,
        "bootstrap_means": {name: draws[plan["estimands"][name]["column"]] for name in ESTIMANDS},
    }


def sampling_figure(name: str, estimate: dict, values: list[float], bootstrap: dict):
    score = name == "mean_score"
    title = "Bootstrap distribution of the mean score" if score else "Bootstrap distribution of mean left thickness"
    figure, axis = make_figure(
        title,
        f"Normal cognition   |   n = {bootstrap['participants_per_replicate']:,}   |   {bootstrap['replicates']:,} participant resamples",
    )
    axis.hist(values, bins=40, density=True, color=GRAYSCALE.accent_light, edgecolor=GRAYSCALE.outline, linewidth=0.55, label="Bootstrap means")
    mean, se = estimate["observed_mean"], estimate["analytic_standard_error"]
    if se:
        normal = NormalDist(mean, se)
        lower = min(min(values), mean - 4 * se)
        upper = max(max(values), mean + 4 * se)
        x = [lower + (upper - lower) * i / 400 for i in range(401)]
        axis.plot(x, [normal.pdf(v) for v in x], color=GRAYSCALE.ink, linewidth=1.4, linestyle="--", label="Normal reference with analytic SE")
    axis.axvline(mean, color=GRAYSCALE.ink, linewidth=1.25, label="Observed mean")
    interval = estimate["percentile_interval"]
    for key in ("lower", "upper"):
        axis.axvline(interval[key], ymax=0.74, color=GRAYSCALE.muted, linewidth=1.1, linestyle=":", label="95% percentile limits" if key == "lower" else None)
    axis.set(
        xlabel="Mean correct animal names in 60 seconds" if score else "Mean left inferior parietal thickness (mm)",
        ylabel="Density (per correct response)" if score else "Density (per mm)",
        ylim=(0, None),
    )
    axis.set_ylim(top=axis.get_ylim()[1] * 1.30)
    axis.xaxis.set_major_locator(MaxNLocator(nbins=7))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.2f}" if score else "{x:.3f}"))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    precision = 5 if score else 6
    annotation = (
        f"Sample SD = {estimate['sample_sd']:.{precision}f}\n"
        f"Analytic SE = {se:.{precision}f}\n"
        f"Bootstrap SE = {estimate['bootstrap_standard_error']:.{precision}f}"
    )
    axis.text(0.98, 0.97, annotation, transform=axis.transAxes, ha="right", va="top", fontsize=9, color=GRAYSCALE.ink)
    axis.legend(loc="upper left", frameon=False, fontsize=8)
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
        for name in ESTIMANDS:
            figure = sampling_figure(name, report["estimates"][name], report["bootstrap_means"][name], report["bootstrap"])
            save_pdf(figure, FIGURES[name], FIGURE_DIR)
            plt.close(figure)
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()
