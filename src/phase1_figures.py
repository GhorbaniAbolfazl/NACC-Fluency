"""Draw the four phase-one figures as grayscale vector PDFs."""

import json
from io import BytesIO
from math import ceil, floor
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, StrMethodFormatter

from .build_dataset import LEFT_THICKNESS, NORMAL_COGNITION, STATUS, write_output
from .data_audit import CATEGORY_LABELS
from .descriptive_statistics import REPORT_FILE, load_participants

FIGURE_DIR = Path("figures/phase1")


class FigureColors(NamedTuple):
    ink: str
    muted: str
    accent: str
    accent_light: str
    neutral: str
    outline: str
    grid: str


GRAYSCALE = FigureColors(
    ink="#262626",
    muted="#525252",
    accent="#575757",
    accent_light="#c5c5c5",
    neutral="#cecece",
    outline="#5f5f5f",
    grid="#e5e5e5",
)


GROUP_LABELS = {
    1: "Normal cognition",
    2: "Impaired without MCI",
    3: "MCI",
    4: "Dementia",
    None: "Unavailable",
}
BOX_LABELS = {**GROUP_LABELS, 1: "Normal\ncognition", 2: "Impaired\nwithout MCI"}


def make_figure(title: str, subtitle: str, grid_axis: str = "y", *, palette: FigureColors = GRAYSCALE):
    figure, axis = plt.subplots(figsize=(7.8, 4.55), layout="constrained")
    figure.suptitle(title, color=palette.ink, fontsize=16, fontweight="bold")
    axis.set_title(subtitle, loc="left", pad=10, color=palette.muted, fontsize=10)
    axis.spines[["top", "right"]].set_visible(False)
    axis.spines[["left", "bottom"]].set_color(palette.outline)
    axis.tick_params(length=3, colors=palette.muted)
    axis.set_axisbelow(True)
    axis.grid(axis=grid_axis, color=palette.grid, linewidth=0.8)
    return figure, axis


def animal_scores(frequency: dict, *, palette: FigureColors = GRAYSCALE):
    counts = {level["value"]: level["count"] for level in frequency["levels"]}
    scores = range(0, max(counts) + 1)
    figure, axis = make_figure(
        "Animal-naming scores",
        f"Normal cognition   |   {frequency['usable']:,} scores available",
        palette=palette,
    )
    axis.bar(scores, [counts.get(score, 0) for score in scores], color=palette.accent, width=0.82)
    axis.set(xlabel="Correct animal names in 60 seconds", ylabel="Participants", xlim=(-0.5, max(scores) + 1.5))
    axis.set_ylim(bottom=0)
    axis.xaxis.set_major_locator(MultipleLocator(5))
    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    return figure


def left_thickness(values: list[float], *, palette: FigureColors = GRAYSCALE):
    # Anchor the bins at multiples of 0.05 mm and include every observed value.
    start = floor(min(values) * 20)
    end = max(start + 1, ceil(max(values) * 20))
    edges = [index / 20 for index in range(start, end + 1)]
    figure, axis = make_figure(
        "Left inferior parietal thickness",
        f"Normal cognition   |   {len(values):,} measurements   |   0.05 mm bins",
        palette=palette,
    )
    axis.hist(values, bins=edges, color=palette.accent_light, edgecolor=palette.accent, linewidth=0.85)
    axis.set(xlabel="Left inferior parietal thickness (mm)", ylabel="Participants", xlim=(edges[0], edges[-1]))
    axis.set_ylim(bottom=0)
    axis.xaxis.set_major_locator(MultipleLocator(0.2))
    axis.xaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    axis.yaxis.set_major_locator(MaxNLocator(integer=True))
    return figure


def cognitive_status(summary: dict, *, palette: FigureColors = GRAYSCALE):
    levels = [(level["code"], level["count"], level["percent_of_all"]) for level in summary["levels"]]
    if summary["unavailable"]:
        levels.append((None, summary["unavailable"], summary["unavailable_percent"]))
    figure, axis = make_figure(
        "Cognitive status in the study sample",
        f"All matched participants   |   N = {summary['total']:,}",
        grid_axis="x",
        palette=palette,
    )
    codes, counts, percentages = zip(*levels)
    bars = axis.barh(
        [GROUP_LABELS[code] for code in codes],
        counts,
        height=0.66,
        color=[palette.accent if code == NORMAL_COGNITION else palette.neutral for code in codes],
        edgecolor=palette.outline,
    )
    axis.bar_label(bars, labels=[f"{count:,} ({percent:.2f}%)" for count, percent in zip(counts, percentages)], padding=7)
    axis.invert_yaxis()
    axis.set(xlabel="Participants", xlim=(0, 1.32 * max(counts)))
    axis.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    return figure


def box_statistics(values: list[float], summary: dict) -> dict:
    review = summary["iqr_review"]
    lower = review["lower_fence"] if review else min(values)
    upper = review["upper_fence"] if review else max(values)
    inside = [value for value in values if lower <= value <= upper]
    return {
        "q1": summary["q1"],
        "med": summary["median"],
        "q3": summary["q3"],
        "whislo": min(inside, default=summary["q1"]),
        "whishi": max(inside, default=summary["q3"]),
        "fliers": [value for value in values if value < lower or value > upper],
    }


def thickness_by_status(rows: list[dict], summaries: dict, *, palette: FigureColors = GRAYSCALE):
    boxes = []
    colors = []
    labels = {**CATEGORY_LABELS[STATUS], None: "unavailable"}
    for code, label in labels.items():
        values = [row[LEFT_THICKNESS] for row in rows if row[STATUS] == code and row[LEFT_THICKNESS] is not None]
        if values:
            summary = summaries[label]["numeric"][LEFT_THICKNESS]
            box = box_statistics(values, summary)
            box["label"] = f"{BOX_LABELS[code]}\n(n = {len(values):,})"
            boxes.append(box)
            colors.append(palette.accent_light if code == NORMAL_COGNITION else palette.neutral)

    figure, axis = make_figure(
        "Left inferior parietal thickness by cognitive status",
        "All matched participants with usable thickness",
        palette=palette,
    )
    plots = axis.bxp(
        boxes,
        patch_artist=True,
        widths=0.55,
        boxprops={"edgecolor": palette.outline, "linewidth": 1.1},
        medianprops={"color": palette.ink, "linewidth": 1.7},
        whiskerprops={"color": palette.outline},
        capprops={"color": palette.outline},
        flierprops={"marker": "o", "markersize": 3.5, "markerfacecolor": "none", "markeredgecolor": palette.outline},
    )
    for patch, color in zip(plots["boxes"], colors):
        patch.set_facecolor(color)
    axis.set(xlabel="Cognitive status", ylabel="Left inferior parietal thickness (mm)")
    axis.yaxis.set_major_locator(MultipleLocator(0.2))
    axis.yaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    return figure


def save_pdf(figure, filename: str, directory: Path = FIGURE_DIR) -> None:
    with BytesIO() as buffer:
        figure.savefig(
            buffer,
            format="pdf",
            metadata={"Title": filename, "CreationDate": None, "ModDate": None},
        )
        write_output(directory / filename, buffer.getvalue())


def main() -> None:
    rows, flow = load_participants()
    descriptive = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    if descriptive["participants_csv_sha256"] != flow["participants_csv_sha256"]:
        raise ValueError("Descriptive statistics are out of date. Run python3 -m src.descriptive_statistics.")

    normal_thickness = [
        row[LEFT_THICKNESS] for row in rows if row[STATUS] == NORMAL_COGNITION and row[LEFT_THICKNESS] is not None
    ]
    animal_frequency = descriptive["discrete_frequencies"]["matched_normal"]["ANIMALS"]
    status_frequency = descriptive["samples"]["matched_all"]["categorical"][STATUS]
    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": GRAYSCALE.ink,
        "axes.labelpad": 9,
        "pdf.fonttype": 42,
        "figure.facecolor": "white",
    }):
        figures = {
            "01_animal_naming.pdf": animal_scores(animal_frequency, palette=GRAYSCALE),
            "02_left_parietal_thickness.pdf": left_thickness(normal_thickness, palette=GRAYSCALE),
            "03_cognitive_status.pdf": cognitive_status(status_frequency, palette=GRAYSCALE),
            "04_thickness_by_status.pdf": thickness_by_status(rows, descriptive["by_cognitive_status"], palette=GRAYSCALE),
        }
        for filename, figure in figures.items():
            save_pdf(figure, filename)
            plt.close(figure)


if __name__ == "__main__":
    main()
