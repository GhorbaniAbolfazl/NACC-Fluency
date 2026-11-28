"""Compare left and mean bilateral thickness on the same participants."""

import json
from collections import Counter
from dataclasses import replace
from decimal import Decimal
from hashlib import sha256
from pathlib import Path
from statistics import fmean, stdev

from .build_dataset import LEFT_THICKNESS, index_participants, write_output
from .data_audit import MRI_FILE, MRI_RULES, assess_value, read_participant_file
from .phase4_design import canonical_hash, prepare_designs
from .phase4_plan import fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary
from .phase5_plan import REPORT_FILE as PHASE5_PLAN, RIGHT_THICKNESS, load_registered_plan

REPORT_FILE = Path("reports/final/bilateral_sensitivity.json")
THICKNESS_TERM = "thickness_0_1mm"


def build_report() -> dict:
    """Keep outcome, adjustment and participants fixed while replacing thickness."""
    mri_path = Path("data") / MRI_FILE
    sources = (PHASE5_PLAN, mri_path, Path("src/phase5_bilateral.py"))
    hashes = {str(path): fingerprint(path) for path in sources}
    plan = load_registered_plan()
    phase4, designs = prepare_designs()
    reference = designs["integrated"]
    specification = plan["comparisons"]["bilateral_thickness"]
    require(list(reference.columns) == plan["reference"]["column_order"]
            and specification["alternative_components"] == [LEFT_THICKNESS, RIGHT_THICKNESS],
            "The common-design columns or the bilateral components differ from the registered plan.")
    thickness_index = reference.columns.index(THICKNESS_TERM)
    center = Decimal(str(plan["reference"]["coding"][THICKNESS_TERM]["reference_mm"]))
    unit = Decimal("0.1")
    require(specification["alternative_coding"] == f"(bilateral_mean_mm - {center}) / {unit}",
            "The bilateral coding differs from the registered center or unit.")
    _, raw_rows = read_participant_file(mri_path, MRI_RULES)
    mri = index_participants(raw_rows, MRI_FILE)

    included = []
    left_mm, right_mm, bilateral_mm, bilateral_column = [], [], [], []
    exclusion_reasons = Counter()
    for i, participant_id in enumerate(reference.participant_ids):
        require(participant_id in mri, "A reference participant is missing from MRI.")
        row = mri[participant_id]
        right, reason = assess_value(row[RIGHT_THICKNESS], MRI_RULES[RIGHT_THICKNESS])
        if right is None:
            exclusion_reasons[reason] += 1
            continue
        left, _ = assess_value(row[LEFT_THICKNESS], MRI_RULES[LEFT_THICKNESS])
        require(left is not None, "A reference row has unusable left thickness.")
        left_value, right_value = Decimal(str(left)), Decimal(str(right))
        require(float((left_value - center) / unit) == reference.matrix[i][thickness_index],
                "The raw left value and registered design disagree.")
        mean = (left_value + right_value) / 2
        included.append(i)
        left_mm.append(left)
        right_mm.append(right)
        bilateral_mm.append(float(mean))
        bilateral_column.append(float((mean - center) / unit))
    n_reference = len(reference.participant_ids)
    n_matched = len(included)
    require(n_reference == plan["reference"]["n"]
            and n_matched == plan["input_support"]["bilateral_support_in_reference_n"],
            "The common sample differs from the registered bilateral support.")
    left = replace(
        reference, name="left_on_bilateral_support",
        participant_ids=tuple(reference.participant_ids[i] for i in included),
        matrix=tuple(reference.matrix[i] for i in included),
        response=tuple(reference.response[i] for i in included),
    )
    bilateral = replace(
        left, name="bilateral_mean_on_same_support",
        matrix=tuple(row[:thickness_index] + (value,) + row[thickness_index + 1:]
                     for row, value in zip(left.matrix, bilateral_column, strict=True)),
    )
    confidence = phase4["estimation"]["confidence_level"]
    models = {}
    for design, measure in ((left, LEFT_THICKNESS), (bilateral, "bilateral_mean_mm")):
        fit = fit_ols(design)
        models[design.name] = {
            "thickness_measure": measure,
            "coefficients": coefficient_table(fit, confidence),
            "fit": fit_summary(fit),
            "matrix_sha256": canonical_hash(design.matrix),
        }
    left_results = models[left.name]
    bilateral_results = models[bilateral.name]
    left_slope = left_results["coefficients"][THICKNESS_TERM]["estimate"]
    bilateral_slope = bilateral_results["coefficients"][THICKNESS_TERM]["estimate"]
    report = {
        "stage": "Phase 5, Step 4",
        "source_files_sha256": hashes,
        "course_requirement": "general_3",
        "sample": {
            "reference_n": n_reference, "matched_n": n_matched,
            "excluded_for_unusable_right_thickness_n": n_reference - n_matched,
            "right_thickness_exclusion_reasons": dict(exclusion_reasons),
            "membership_sha256": sha256(("\n".join(left.participant_ids) + "\n").encode("utf-8")).hexdigest(),
            "response_sha256": canonical_hash(left.response),
        },
        "specification": {
            "formula": plan["reference"]["model"],
            "column_order": list(reference.columns),
            "thickness_column": THICKNESS_TERM,
            "left_coding": plan["reference"]["coding"][THICKNESS_TERM]["expression"],
            "bilateral_coding": specification["alternative_coding"],
            "bilateral_definition": specification["alternative_definition"],
            "only_design_change": specification["single_change"],
        },
        "thickness_summary_mm": {
            name: {"mean": fmean(values), "sample_sd": stdev(values)}
            for name, values in (("left", left_mm), ("right", right_mm), ("bilateral_mean", bilateral_mm))
        },
        "models": models,
        "bilateral_minus_left": {
            "thickness_estimate": bilateral_slope - left_slope,
            "R_squared": bilateral_results["fit"]["R_squared"] - left_results["fit"]["R_squared"],
            "in_sample_RMSE": bilateral_results["fit"]["in_sample_RMSE"] - left_results["fit"]["in_sample_RMSE"],
        },
        "interpretation": {
            "slope_unit": "Difference in correct animal responses per 0.1 mm of the specified thickness measure, holding age, recorded sex and education fixed.",
            "anatomical_scope": "The bilateral measure averages both hemispheres. This comparison does not estimate a right-only slope or a difference between hemisphere-specific slopes.",
            "interval_scope": "The individual nominal 95% intervals describe each slope. The observed shift between the slopes describes the two fits on this sample and has no interval here. The intervals cannot establish equivalence.",
            "fit_scope": "Fit metrics summarize the observations used for estimation. The intercepts hold different anatomical measures at the same numeric center.",
        },
        "privacy": "Only aggregate results and hashes are persisted.",
    }
    require(hashes == {str(path): fingerprint(path) for path in sources},
            "An input changed during the bilateral comparison.")
    return report


def main() -> None:
    data = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != data:
        write_output(REPORT_FILE, data)


if __name__ == "__main__":
    main()
