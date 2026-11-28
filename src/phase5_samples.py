"""Compare the registered samples and audit MRI coverage and visit information."""

import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
from statistics import fmean, stdev

from .build_dataset import (
    LEFT_THICKNESS, NORMAL_COGNITION, STATUS, build_dataset,
    index_participants, summarize_sample, write_output,
)
from .data_audit import (
    CATEGORY_LABELS, MRI_FILE, MRI_RULES, PARTICIPANT_ID, UDS_FILE, UDS_RULES,
    assess_value, group_mri, group_uds, read_participant_file,
)
from .phase4_design import FEATURE_ORDER, ModelDesign, encode_predictors
from .phase4_plan import REPORT_FILE as PHASE4_PLAN, fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary
from .phase5_plan import (
    PRIMARY_FILE, REPORT_FILE as PHASE5_PLAN, load_registered_plan, read_json,
)

REPORT_FILE = Path("reports/final/sample_sensitivity.json")
UDS_FIELDS = ("ANIMALS", "NACCAGE", "SEX", "EDUC")
THICKNESS_TERM = "thickness_0_1mm"


def composition(rows: list[dict]) -> dict:
    """Describe participants by cognitive status, recorded sex and UDS values."""
    uds_summary = {}
    for column in ("ANIMALS", "NACCAGE", "EDUC"):
        values = [row[column] for row in rows]
        uds_summary[column] = {"mean": fmean(values), "sample_sd": stdev(values)}
    return {
        "n": len(rows),
        "by_cognitive_status": {
            label: sum(row[STATUS] == code for row in rows)
            for code, label in CATEGORY_LABELS[STATUS].items()
        },
        "by_recorded_sex": {
            label: sum(row["SEX"] == code for row in rows)
            for code, label in CATEGORY_LABELS["SEX"].items()
        },
        "uds_summary": uds_summary,
    }


def audit_mri_and_visits(reference_ids: tuple[str, ...], joined: list[dict]) -> dict:
    """Define the UDS denominator before using any MRI measurements."""
    uds_header, uds_rows = read_participant_file(Path("data") / UDS_FILE, UDS_RULES)
    mri_header, mri_rows = read_participant_file(Path("data") / MRI_FILE, MRI_RULES)
    uds = index_participants(uds_rows, UDS_FILE)
    mri = index_participants(mri_rows, MRI_FILE)
    uds_groups = group_uds(uds_header)
    mri_groups = group_mri(mri_header)
    require(not uds_groups.get("uncategorized"), "Review new UDS columns for timing or visit metadata.")

    standardized = {
        participant_id: {field: assess_value(raw[field], UDS_RULES[field])[0]
                         for field in (STATUS, *UDS_FIELDS)}
        for participant_id, raw in uds.items()
    }
    uds_eligibility = summarize_sample(list(standardized.values()), "uds_only", UDS_FIELDS, True)
    before = {key for key, row in standardized.items() if row["eligible_uds_only"]}
    matched = before & mri.keys()
    after = {key for key in matched
             if assess_value(mri[key][LEFT_THICKNESS], MRI_RULES[LEFT_THICKNESS])[0] is not None}
    require(tuple(sorted(after)) == reference_ids,
            "The MRI-restricted UDS sample differs from the frozen reference.")
    require(before == after,
            "The frozen extract no longer supports the registered identical-sample UDS comparison.")
    joined_by_id = {row[PARTICIPANT_ID]: row for row in joined}
    require(all(standardized[key][field] == joined_by_id[key][field]
                for key in reference_ids for field in (STATUS, *UDS_FIELDS)),
            "The UDS-only fields differ between original UDS and joined records.")

    # The verified Phase 4 fit applies to both definitions because the design and response are identical.
    uds_only = read_json(PRIMARY_FILE)["models"]["uds_only"]
    return {
        "mri_completeness": {
            "within_extract": True,
            "normal_cognition_uds_eligible": uds_eligibility,
            "normal_cognition_before_uds_filters_n": sum(row[STATUS] == NORMAL_COGNITION
                                                        for row in standardized.values()),
            "with_linked_mri_n": len(matched),
            "with_linked_usable_left_thickness_n": len(after),
            "without_linked_mri_n": len(before - matched),
            "linked_but_without_usable_left_thickness_n": len(matched - after),
            "membership_identical": True,
            "uds_fields_identical": True,
            "coverage_fraction_among_uds_eligible": len(after) / len(before),
            "composition_before_and_after": composition([standardized[key] for key in reference_ids]),
            "uds_only_fit_for_both": {
                "formula": uds_only["formula"], "coefficients": uds_only["coefficients"],
                "fit": uds_only["fit"], "source": str(PRIMARY_FILE),
                "basis": "The revalidated reference fit applies to both definitions. Identical participants, outcomes and covariates yield zero changes in the estimates and fit metrics.",
            },
            "limit": "Complete coverage of this thickness variable within the teaching extract cannot measure the effect of MRI availability in the upstream NACC population or establish clinical scan quality.",
        },
        "visit_rule": {
            "uds": {"rows": len(uds_rows), "distinct_participants": len(uds), "columns_reviewed": len(uds_header),
                    "maximum_supplied_rows_per_participant": 1},
            "mri": {"rows": len(mri_rows), "distinct_participants": len(mri), "columns_reviewed": len(mri_header),
                    "maximum_supplied_rows_per_participant": 1,
                    "participant_fields": mri_groups["participant_fields"]},
            "matched_participants": len(uds.keys() & mri.keys()),
            "clinical_visit_date_available": False, "clinical_visit_number_available": False,
            "mri_scan_date_available": False, "source_visit_count_available": False,
            "mri_scan_identifier_available": False,
            "alternative_visit_selection_evaluable": False,
            "temporal_proximity_evaluable": False,
            "basis": "All UDS fields map to the existing documented content groups. MRI fields consist of the participant identifier, cognitive status and regional measurements. The supplied pairing is by participant identifier.",
            "limit": "One supplied pair per person establishes neither the number of original encounters nor first-visit, closest-scan or same-day status.",
        },
    }


def build_report() -> dict:
    """Change only the cognitive-status inclusion rule in the broader model."""
    sources = (PHASE5_PLAN, PHASE4_PLAN, PRIMARY_FILE, Path("data") / UDS_FILE,
               Path("data") / MRI_FILE, Path("src/phase5_samples.py"))
    hashes = {str(path): fingerprint(path) for path in sources}
    plan = load_registered_plan()
    phase4 = read_json(PHASE4_PLAN)
    rows, flow = build_dataset(Path("data"))
    normal = [row for row in rows if row["eligible_model_normal"]]
    broad = [row for row in rows if row["eligible_model_all"]]
    reference_ids = tuple(row[PARTICIPANT_ID] for row in normal)
    require(sha256(("\n".join(reference_ids) + "\n").encode("utf-8")).hexdigest()
            == plan["reference"]["membership_sha256"], "The reference membership changed.")
    require(len(normal) == plan["reference"]["n"] and len(broad) == plan["input_support"]["broader_n"],
            "A sample count differs from the registered Phase 5 plan.")
    require([row for row in broad if row[STATUS] == NORMAL_COGNITION] == normal,
            "The broad fit must retain every reference participant and the same reference values.")

    broad_composition = composition(broad)
    require(broad_composition["by_cognitive_status"] == {
        item["label"]: item["n"] for item in plan["input_support"]["broader_cognitive_status_counts"].values()
    }, "The broader cognitive-status composition changed.")
    audit = audit_mri_and_visits(reference_ids, rows)
    reference_model = read_json(PRIMARY_FILE)["models"]["integrated"]
    columns = tuple(plan["reference"]["column_order"])
    require(columns == FEATURE_ORDER, "The integrated predictor order changed.")
    design = ModelDesign(
        "broader_integrated", tuple(row[PARTICIPANT_ID] for row in broad), columns,
        tuple(encode_predictors(row, phase4) for row in broad),
        tuple(float(row["ANIMALS"]) for row in broad), tuple(reference_model["reference_row"]),
    )
    fit = fit_ols(design)
    coefficients = coefficient_table(fit, phase4["estimation"]["confidence_level"])
    normal_coefficients = reference_model["coefficients"]
    shift = coefficients[THICKNESS_TERM]["estimate"] - normal_coefficients[THICKNESS_TERM]["estimate"]
    status_totals = Counter(row[STATUS] for row in rows)
    report = {
        "stage": "Phase 5, Step 3",
        "source_files_sha256": hashes,
        "sample_definition": {
            "formula": plan["reference"]["model"], "column_order": list(columns),
            "predictor_coding": plan["reference"]["coding"],
            "linked_participants": flow["linkage"]["matched_participants"],
            "linked_by_cognitive_status": {label: status_totals[code]
                                           for code, label in CATEGORY_LABELS[STATUS].items()},
            "reference": {"eligibility": flow["samples"]["model_normal"],
                          "composition": composition(normal),
                          "coefficients": normal_coefficients, "fit": reference_model["fit"]},
            "broader": {"eligibility": flow["samples"]["model_all"], "composition": broad_composition,
                        "coefficients": coefficients, "fit": fit_summary(fit)},
            "additional_participants_n": len(broad) - len(normal),
            "thickness_estimate_change_broader_minus_reference": shift,
            "interpretation": "The broader common slope pools the recorded cognitive categories without adjusting for status. The difference is descriptive: the fits share the reference participants and target different population mixes.",
            "interval_scope": "Individual nominal 95% intervals within each fit. They are not intervals for the change between the overlapping fits.",
        },
        **audit,
        "course_requirements": {"general_2": "broader integrated fit", "track_c2": "within-extract MRI coverage audit and UDS-only comparison",
                                "track_c3": "visit and timing feasibility audit"},
        "privacy": "Only aggregate results and source hashes are persisted.",
    }
    require(hashes == {str(path): fingerprint(path) for path in sources},
            "An input changed during the sample-sensitivity analysis.")
    return report


def main() -> None:
    data = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != data:
        write_output(REPORT_FILE, data)


if __name__ == "__main__":
    main()
