"""Describe animal-naming scores and left parietal thickness in paired participants."""

import json
from math import fsum, sqrt
from pathlib import Path
from statistics import NormalDist, pvariance, quantiles, variance

from .build_dataset import LEFT_THICKNESS, PAIR_FIELDS, write_output
from .descriptive_statistics import empirical_frequencies, numeric_summary
from .phase1_figures import GRAYSCALE, left_thickness, make_figure, save_pdf
from .phase2_samples import load_phase2_samples

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, PercentFormatter

SCORE, _ = PAIR_FIELDS
SAMPLE_NAME = "paired_normal"
REPORT_FILE = Path("reports/phase2/distributions.json")
FIGURE_DIR = Path("figures/phase2")


def summarize_distributions(rows: list[dict]) -> dict:
    """Use the same paired sample for both variables' empirical distributions."""
    summaries = {}
    for column in PAIR_FIELDS:
        values = [row[column] for row in rows]
        summary = numeric_summary(rows, column)
        empirical_variance = pvariance(values)
        sd = sqrt(empirical_variance)
        skewness = None
        if sd > 0:
            skewness = fsum(((value - summary["mean"]) / sd) ** 3 for value in values) / len(values)
        summaries[column] = {
            **summary,
            "empirical_variance": empirical_variance,
            "sample_variance": variance(values) if len(values) >= 2 else None,
            "empirical_skewness": skewness,
        }

    score_frequency = empirical_frequencies(rows, SCORE)
    for level in score_frequency["levels"]:
        level["probability"] = level["count"] / score_frequency["usable"]

    thickness = [row[LEFT_THICKNESS] for row in rows]
    thickness_summary = summaries[LEFT_THICKNESS]
    reference = None
    if thickness_summary["empirical_variance"] > 0:
        normal = NormalDist(thickness_summary["mean"], sqrt(thickness_summary["empirical_variance"]))
        empirical_percentiles = quantiles(thickness, n=100, method="inclusive")
        reference = {
            "mean_mm": normal.mean,
            "standard_deviation_mm": normal.stdev,
            "quantile_checks": [
                {
                    "percentile": percentile,
                    "empirical_mm": empirical_percentiles[percentile - 1],
                    "fitted_normal_mm": normal.inv_cdf(percentile / 100),
                }
                for percentile in (1, 5, 50, 95, 99)
            ],
        }

    return {
        "sample_size": len(rows),
        "units": {SCORE: "correct responses in 60 seconds", LEFT_THICKNESS: "mm"},
        "methods": {
            "empirical_variance": "Second central moment with denominator n",
            "sample_variance": "Sum of squared deviations from the sample mean divided by n - 1",
            "empirical_skewness": "Third central moment divided by the empirical standard deviation cubed",
            "empirical_quantiles": "Linear interpolation at zero-based rank (n - 1) * p",
            "outliers": (
                "Review values strictly outside the 1.5 IQR fences when n >= 4 and IQR > 0. "
                "Keep these values in the analysis."
            ),
            "normal_reference": "Normal reference matched to the empirical mean and empirical standard deviation",
            "qq_plotting_positions": "(i - 0.5) / n for sorted observations ranked i = 1 through n",
        },
        "summaries": summaries,
        "animal_score_frequency": score_frequency,
        "thickness_normal_reference": reference,
    }


def animal_score_pmf(frequency: dict):
    """Draw the empirical probability of each integer score."""
    probabilities = {level["value"]: level["probability"] for level in frequency["levels"]}
    scores = range(max(probabilities) + 1)
    figure, axis = make_figure(
        "Empirical distribution of animal-naming scores",
        f"Normal cognition   |   {frequency['usable']:,} paired participants",
        palette=GRAYSCALE,
    )
    axis.bar(scores, [probabilities.get(score, 0) for score in scores], width=0.82, color=GRAYSCALE.accent)
    axis.set(
        xlabel="Correct animal names in 60 seconds",
        ylabel="Empirical probability",
        xlim=(-0.5, max(scores) + 1.5),
        ylim=(0, None),
    )
    axis.xaxis.set_major_locator(MultipleLocator(5))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    return figure


def thickness_normal_qq(thickness: list[float], reference: dict | None):
    """Compare ordered measurements with corresponding fitted normal quantiles."""
    figure, axis = make_figure(
        "Normal Q-Q plot of left parietal thickness",
        f"Normal cognition   |   {len(thickness):,} paired participants",
        grid_axis="both",
        palette=GRAYSCALE,
    )
    axis.set(xlabel="Fitted normal quantiles (mm)", ylabel="Observed thickness quantiles (mm)")
    if reference is None:
        axis.text(
            0.5, 0.5, "A normal reference requires variation in thickness.",
            transform=axis.transAxes, ha="center", va="center",
        )
        axis.set_xticks([])
        axis.set_yticks([])
        return figure

    normal = NormalDist(reference["mean_mm"], reference["standard_deviation_mm"])
    observed = sorted(thickness)
    theoretical = [normal.inv_cdf((rank - 0.5) / len(observed)) for rank in range(1, len(observed) + 1)]
    lower = min(observed[0], theoretical[0])
    upper = max(observed[-1], theoretical[-1])
    pad = (upper - lower) * 0.035
    limits = (lower - pad, upper + pad)
    axis.plot(limits, limits, color=GRAYSCALE.outline, linewidth=1.1, linestyle="--", label="Normal reference")
    axis.scatter(theoretical, observed, color=GRAYSCALE.ink, s=8, linewidths=0, label="Observed")
    axis.set(xlim=limits, ylim=limits, aspect="equal")
    axis.xaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.legend(loc="lower right", frameon=False)
    return figure


def main() -> None:
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    thickness = [row[LEFT_THICKNESS] for row in rows]
    results = summarize_distributions(rows)
    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "sample": SAMPLE_NAME,
        **results,
    }
    output = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"

    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": GRAYSCALE.ink,
        "axes.labelpad": 9,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
    }):
        figures = {
            "01_animal_score_pmf.pdf": animal_score_pmf(results["animal_score_frequency"]),
            "02_thickness_histogram.pdf": left_thickness(thickness, palette=GRAYSCALE),
            "03_thickness_normal_qq.pdf": thickness_normal_qq(thickness, results["thickness_normal_reference"]),
        }
        for filename, figure in figures.items():
            save_pdf(figure, filename, FIGURE_DIR)
            plt.close(figure)
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()
