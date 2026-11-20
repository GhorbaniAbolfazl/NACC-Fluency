"""Generate Phase 1 LaTeX table data and default English labels."""

import json
from pathlib import Path

from .build_dataset import LEFT_THICKNESS, write_output

PHASE_DIR = Path("reports/phase1")
LABELS = {
    "UDS": "UDS",
    "MRI": "MRI",
    "Participant": "Participant",
    "Matrix": "Adjacency matrix",
    "Region": "Brain region",
    "Animals": "Animal-naming score",
    "AnimalsShort": "Animal-naming score",
    "Vegetables": "Vegetable-naming score",
    "Age": "Age at visit",
    "Sex": "Recorded sex",
    "Education": "Years of education",
    "Thickness": "Left cortical thickness",
    "Matched": "Matched participant table",
    "StatusAndMRI": "Usable cognitive status and thickness",
    "PairedAll": "Paired score and thickness, all groups",
    "PairedNormal": "Paired score and thickness, normal cognition",
    "ModelAll": "Supplementary model, all groups",
    "ModelNormal": "Primary model, normal cognition",
    "ThicknessMM": "Left thickness (mm)",
    "AgeYears": "Age (years)",
    "EducationYears": "Education (years)",
    "Normal": "Normal cognition",
    "ImpairedWithoutMCI": "Impaired without MCI",
    "MCI": "MCI",
    "Dementia": "Dementia",
}


def label(key: str) -> str:
    return rf"\Table{key}{{}}"


def format_number(value: int | float | None, decimals: int = 0) -> str:
    if value is None:
        return r"\ReportMissing{}"
    return rf"\ReportNumber{{{value:,.{decimals}f}}}"


def command(name: str, body: str) -> str:
    return rf"\newcommand{{\{name}}}{{{body}}}"


def table_rows(name: str, rows: list[tuple[str, ...]]) -> str:
    body = "\n".join(" & ".join(row) + r" \\" for row in rows)
    return command(name, "%\n" + body + "\n")


def main() -> None:
    audit, flow, stats = (
        json.loads((PHASE_DIR / f"{filename}.json").read_text(encoding="utf-8"))
        for filename in ("data_audit", "sample_flow", "descriptive_statistics")
    )
    if stats["participants_csv_sha256"] != flow["participants_csv_sha256"]:
        raise ValueError("Descriptive report does not match the sample flow")

    linkage = flow["linkage"]
    samples = flow["samples"]
    matched_stats = stats["samples"]["matched_all"]
    normal_stats = stats["samples"]["matched_normal"]["numeric"]
    normal_exclusions = samples["model_normal"]["first_exclusion_reason"]
    score = normal_stats["ANIMALS"]
    thickness = normal_stats[LEFT_THICKNESS]
    score_review = score["iqr_review"]
    thickness_review = thickness["iqr_review"]
    scalars = {
        "NMatched": (linkage["matched_participants"], 0),
        "NNormal": (flow["matched_normal_cognition_before_eligibility"], 0),
        "NPairedNormal": (samples["paired_normal"]["included"], 0),
        "NModelNormal": (samples["model_normal"]["included"], 0),
        "NOnlyUDS": (linkage["uds_only"], 0),
        "NOnlyMRI": (linkage["mri_only"], 0),
        "NStatusEqual": (linkage["cognitive_status_agreement"]["equal"], 0),
        "NNoScoreNormal": (normal_exclusions.get("ANIMALS", 0), 0),
        "NNoEducationNormal": (normal_exclusions.get("EDUC", 0), 0),
        "NormalScoreMedian": (score["median"], 0),
        "NormalThicknessMedian": (thickness["median"], 3),
        "NormalScoreMin": (score["minimum"], 0),
        "NormalScoreMax": (score["maximum"], 0),
        "NormalThicknessMin": (thickness["minimum"], 3),
        "NormalThicknessMax": (thickness["maximum"], 3),
        "NormalScoreFlags": (score_review["below"] + score_review["above"] if score_review else None, 0),
        "NormalThicknessFlags": (
            thickness_review["below"] + thickness_review["above"] if thickness_review else None, 0
        ),
        "NormalAgeMean": (normal_stats["NACCAGE"]["mean"], 2),
        "MCIAgeMean": (stats["by_cognitive_status"]["MCI"]["numeric"]["NACCAGE"]["mean"], 2),
        "NWomen": (next(item["count"] for item in matched_stats["categorical"]["SEX"]["levels"] if item["code"] == 2), 0),
        "NSex": (matched_stats["categorical"]["SEX"]["usable"], 0),
    }

    sources = [
        (label("UDS"), format_number(audit["uds"]["rows"]), format_number(audit["uds"]["columns"]), label("Participant")),
        (label("MRI"), format_number(audit["mri"]["rows"]), format_number(audit["mri"]["columns"]), label("Participant")),
        (
            label("Matrix"),
            format_number(audit["adjacency_matrix"]["rows"]),
            format_number(audit["adjacency_matrix"]["data_columns"] + 1),
            label("Region"),
        ),
    ]

    quality_fields = {
        "Animals": audit["uds"]["selected_variables"]["ANIMALS"],
        "Vegetables": audit["uds"]["selected_variables"]["VEG"],
        "Age": audit["uds"]["selected_variables"]["NACCAGE"],
        "Sex": audit["uds"]["selected_variables"]["SEX"],
        "Education": audit["uds"]["selected_variables"]["EDUC"],
        "Thickness": audit["mri"]["selected_variables"][LEFT_THICKNESS],
    }
    quality = [
        (
            label(key),
            format_number(values["usable"]),
            format_number(values["unusable"]),
            format_number(values["unusable_percent"], 2) + r"\%",
        )
        for key, values in quality_fields.items()
    ]

    counts = [
        ("Matched", linkage["matched_participants"]),
        ("StatusAndMRI", samples["status_and_mri"]["included"]),
        ("PairedAll", samples["paired_all"]["included"]),
        ("PairedNormal", samples["paired_normal"]["included"]),
        ("ModelAll", samples["model_all"]["included"]),
        ("ModelNormal", samples["model_normal"]["included"]),
    ]

    descriptive = []
    fields = [
        ("AnimalsShort", "ANIMALS", 0),
        ("ThicknessMM", LEFT_THICKNESS, 3),
        ("AgeYears", "NACCAGE", 0),
        ("EducationYears", "EDUC", 0),
    ]
    for key, field, quartile_decimals in fields:
        values = normal_stats[field]
        mean_decimals = 3 if field == LEFT_THICKNESS else 2
        q1 = format_number(values["q1"], quartile_decimals)
        q3 = format_number(values["q3"], quartile_decimals)
        descriptive.append((
            label(key),
            format_number(values["usable"]),
            format_number(values["mean"], mean_decimals),
            format_number(values["sample_sd"], mean_decimals),
            format_number(values["median"], quartile_decimals),
            rf"\ReportRange{{{q1}}}{{{q3}}}" if values["usable"] else r"\ReportMissing{}",
            format_number(values["iqr"], quartile_decimals),
        ))

    groups = []
    for group_name, key in (
        ("normal cognition", "Normal"),
        ("impaired without MCI", "ImpairedWithoutMCI"),
        ("MCI", "MCI"),
        ("dementia", "Dementia"),
    ):
        group = stats["by_cognitive_status"][group_name]
        animal, thickness = (group["numeric"][field] for field in ("ANIMALS", LEFT_THICKNESS))
        groups.append((
            label(key),
            format_number(group["participants"]),
            format_number(animal["usable"]),
            format_number(animal["mean"], 2),
            format_number(thickness["usable"]),
            format_number(thickness["mean"], 3),
        ))

    lines = [
        "% Generated by python3 -m src.phase1_tables",
        "% Localized reports can override these labels and number formats before loading this file.",
        r"\providecommand{\ReportNumber}[1]{#1}",
        r"\providecommand{\ReportRange}[2]{#1--#2}",
        r"\providecommand{\ReportMissing}{Not available}",
    ]
    lines.extend(rf"\providecommand{{\Table{key}}}{{{text}}}" for key, text in LABELS.items())
    lines.extend(command(name, format_number(value, digits)) for name, (value, digits) in scalars.items())
    lines.extend([
        table_rows("SourceRows", sources),
        table_rows("QualityRows", quality),
        table_rows("FlowRows", [(label(key), format_number(count)) for key, count in counts]),
        table_rows("DescriptiveRows", descriptive),
        table_rows("GroupRows", groups),
    ])
    write_output(PHASE_DIR / "tables.tex", ("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
