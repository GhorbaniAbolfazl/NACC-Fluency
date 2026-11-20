"""Describe the unadjusted association between paired thickness and fluency scores."""

import json
from pathlib import Path
from statistics import covariance, fmean

from .build_dataset import LEFT_THICKNESS, write_output
from .phase1_figures import GRAYSCALE, make_figure, save_pdf
from .phase2_distributions import FIGURE_DIR, SAMPLE_NAME, SCORE
from .phase2_samples import load_phase2_samples

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, StrMethodFormatter

REPORT_FILE = Path("reports/phase2/paired_covariance.json")
FIGURE_NAME = "04_thickness_and_animal_scores.pdf"


def summarize_covariance(thickness: list[float], scores: list[int]) -> dict:
    """Summarize centered cross-products using both n - 1 and n denominators."""
    n = len(thickness)
    if not n or n != len(scores):
        raise ValueError("Covariance requires paired thickness and score measurements.")

    sample_covariance = covariance(thickness, scores) if n > 1 else None
    empirical_covariance = sample_covariance * (n - 1) / n if sample_covariance is not None else 0.0
    direction = "zero"
    if empirical_covariance > 0:
        direction = "positive"
    elif empirical_covariance < 0:
        direction = "negative"

    return {
        "sample_size": n,
        "scope": "Unadjusted, descriptive association in participants with normal cognition",
        "variables": {
            "x": {"column": LEFT_THICKNESS, "unit": "mm", "mean": fmean(thickness)},
            "y": {"column": SCORE, "unit": "correct responses in 60 seconds", "mean": fmean(scores)},
        },
        "covariance": {
            "unit": "mm times correct responses",
            "sample": {"value": sample_covariance, "denominator": n - 1},
            "empirical": {"value": empirical_covariance, "denominator": n},
            "direction": direction,
        },
        "scatter": {
            "file": str(FIGURE_DIR / FIGURE_NAME),
            "observations": n,
            "distinct_coordinates": len(set(zip(thickness, scores))),
        },
    }


def paired_scatter(thickness: list[float], scores: list[int]):
    """Plot the unmodified measurements, with transparency to show overlap."""
    figure, axis = make_figure(
        "Left parietal thickness and animal-naming score",
        f"Normal cognition   |   {len(thickness):,} paired participants",
        grid_axis="both",
        palette=GRAYSCALE,
    )
    axis.scatter(thickness, scores, s=16, color=GRAYSCALE.ink, alpha=0.45, linewidths=0)
    axis.set(
        xlabel="Left inferior parietal cortical thickness (mm)",
        ylabel="Correct animal names in 60 seconds",
    )
    axis.xaxis.set_major_locator(MultipleLocator(0.2))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    axis.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=8))
    return figure


def main() -> None:
    samples, flow = load_phase2_samples()
    rows = samples[SAMPLE_NAME]
    thickness = [row[LEFT_THICKNESS] for row in rows]
    scores = [row[SCORE] for row in rows]
    report = {
        "participants_csv_sha256": flow["participants_csv_sha256"],
        "sample": SAMPLE_NAME,
        **summarize_covariance(thickness, scores),
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
        figure = paired_scatter(thickness, scores)
        save_pdf(figure, FIGURE_NAME, FIGURE_DIR)
        plt.close(figure)
    write_output(REPORT_FILE, output.encode("utf-8"))


if __name__ == "__main__":
    main()
