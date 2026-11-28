"""Freeze the Phase 5 sensitivity comparisons before fitting alternative models."""

import json
from collections import Counter
from math import isclose
from pathlib import Path

from .build_dataset import LEFT_THICKNESS, STATUS, build_dataset, write_output
from .data_audit import CATEGORY_LABELS
from .phase4_design import canonical_hash, prepare_designs
from .phase4_plan import fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary

REPORT_FILE = Path("reports/final/analysis_plan.json")
PHASE4_REVIEW = Path("reports/phase4/regression_review.json")
PRIMARY_FILE = Path("reports/phase4/multiple_regression.json")
SIMPLE_FILE = Path("reports/phase4/simple_regression.json")
RIGHT_THICKNESS = "RH_INFERIORPARIETAL_AVGTH"
BOOTSTRAP_DRAWS = 10_000
BOOTSTRAP_SEED = 731_905
REFERENCE_REL_TOL = 2e-10
REFERENCE_ABS_TOL = 2e-11
INFLUENCE_MAX = 5
MCI_CODE = 3

REQUIREMENTS = {
    "general_1": (2, "Review at least three relevant studies in a structured comparison table."),
    "general_2": (3, "Repeat the reference specification on the broader eligible sample."),
    "general_3": (4, "Compare left-only and bilateral thickness on common support."),
    "general_4": (5, "Compare the complete reference fit with a temporary omission of the largest flagged Cook distances."),
    "general_5": (5, "Compare the standard and participant-bootstrap intervals for the primary coefficient."),
    "general_6": (6, "Compare normal cognition and MCI in one planned exploratory subgroup analysis."),
    "general_7": (7, "Integrate evidence from descriptive, probability, inferential and regression phases."),
    "general_8": (7, "Explain limitations and identify the claims supported by the data."),
    "track_c1": (7, "Present the final three-model comparison on the frozen common sample."),
    "track_c2": (3, "Audit normal-cognition UDS eligibility with and without the usable-MRI requirement."),
    "track_c3": (3, "Check availability of repeated visits and timing metadata for a visit-choice sensitivity analysis."),
    "track_c4": (6, "Compare normal cognition and MCI using one record per participant."),
    "track_c5": (7, "Summarize the incremental contribution and limits of integrating UDS and MRI."),
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def eligible_support() -> dict:
    """Record aggregate feasibility counts without modeling alternative samples."""
    rows, flow = build_dataset(Path("data"))
    normal = [row for row in rows if row["eligible_model_normal"]]
    broad = [row for row in rows if row["eligible_model_all"]]
    categories = Counter(row[STATUS] for row in broad)
    status_counts = {str(code): {"label": label, "n": categories[code]}
                     for code, label in CATEGORY_LABELS[STATUS].items()}
    require(sum(item["n"] for item in status_counts.values()) == len(broad),
            "Every participant in the planned broad sample must have a recorded cognitive category.")
    linkage = flow["linkage"]
    require(linkage["uds_only"] == linkage["mri_only"] == 0, "The expected fully matched teaching extract changed.")
    require(len(normal) == 842 and len(broad) == 1223 and categories[MCI_CODE] == 248,
            "The reference, broad or MCI sample count differs from the planned comparison.")
    return {
        "linked_rows": linkage["matched_participants"],
        "uds_only_rows": linkage["uds_only"],
        "mri_only_rows": linkage["mri_only"],
        "linked_rows_without_usable_left_thickness": linkage["matched_with_unusable_left_thickness"],
        "reference_n": len(normal),
        "broader_n": len(broad),
        "broader_cognitive_status_counts": status_counts,
        "bilateral_support_in_reference_n": sum(row[RIGHT_THICKNESS] is not None for row in normal),
        "normal_cognition_vs_mci_n": {"normal_cognition": len(normal), "mci": categories[MCI_CODE]},
        "basis": "Use the Phase 1 value rules on the source teaching extract. Usable supplied values do not establish clinical scan quality or whether a UDS value was imputed upstream.",
    }


def verify_reference() -> tuple[dict, dict]:
    """Compare only the saved Phase 4 fits using their frozen designs."""
    phase4, designs = prepare_designs()
    multiple = read_json(PRIMARY_FILE)
    simple = read_json(SIMPLE_FILE)

    def compare_numeric(actual: dict, expected: dict, context: str) -> None:
        require(actual.keys() == expected.keys(), f"Reference keys changed: {context}.")
        for key, value in actual.items():
            recorded = expected[key]
            if isinstance(value, dict):
                compare_numeric(value, recorded, f"{context}.{key}")
            elif isinstance(value, float):
                require(isclose(value, recorded, rel_tol=REFERENCE_REL_TOL, abs_tol=REFERENCE_ABS_TOL),
                        f"Reference value changed: {context}.{key}.")
            else:
                require(value == recorded, f"Reference value changed: {context}.{key}.")

    for name, design in designs.items():
        saved = simple if name == "mri_only" else multiple["models"][name]
        design_hashes = saved["design"] if name == "mri_only" else saved
        require(design_hashes["matrix_sha256"] == canonical_hash(design.matrix)
                and design_hashes["response_sha256"] == canonical_hash(design.response)
                and design_hashes["membership_sha256"] == phase4["sample"]["membership_sha256"],
                f"The reference design or sample changed: {name}.")
        fit = fit_ols(design)
        compare_numeric(fit_summary(fit), saved["fit"], f"{name}.fit")
        compare_numeric(coefficient_table(fit, phase4["estimation"]["confidence_level"]),
                        saved["coefficients"], f"{name}.coefficients")
    return phase4, multiple["primary_test"]


def build_plan() -> dict:
    """Verify the frozen Phase 4 evidence and register Phase 5 choices."""
    review = read_json(PHASE4_REVIEW)
    fingerprints = {**review["source_files_sha256"], str(PHASE4_REVIEW): fingerprint(PHASE4_REVIEW),
                    "src/phase5_plan.py": fingerprint(Path("src/phase5_plan.py"))}
    require(review["stage"] == "Phase 4, Step 8" and len(review["coverage"]) == 17,
            "A completed Phase 4 review is required.")
    require(fingerprints == {name: fingerprint(Path(name)) for name in fingerprints},
            "Phase 4 artifacts or inputs differ from the frozen review.")
    phase4, primary = verify_reference()
    support = eligible_support()
    require(primary["model"] == "integrated" and primary["term"] == "thickness_0_1mm"
            and review["participants"] == support["reference_n"],
            "The recorded reference differs from the integrated Phase 4 thickness estimand.")
    plan = {
        "plan_version": 1,
        "stage": "Phase 5, Step 1",
        "status": "Reference reproduced and comparisons specified before Phase 5 sensitivity fits.",
        "analysis_character": "Protocol for robustness work on the already examined teaching extract. The Phase 4 reference and diagnostics were known when the Phase 5 comparisons were specified.",
        "source_files_sha256": fingerprints,
        "reference": {
            "question": phase4["question"],
            "sample": phase4["sample"]["name"],
            "n": support["reference_n"],
            "membership_sha256": phase4["sample"]["membership_sha256"],
            "model": phase4["models"]["integrated"]["formula"],
            "column_order": phase4["models"]["integrated"]["columns"],
            "response": phase4["variables"]["response"],
            "coding": phase4["variables"]["explanatory"],
            "estimand": primary["parameter"],
            "estimate": primary["estimate"],
            "model_based": primary["model_based"],
            "hc3_companion": primary["hc3_companion"],
            "source_review": str(PHASE4_REVIEW),
            "verification": {
                "file_checks": "Match the recorded Phase 4 artifact and input hashes and require identical sample and design hashes.",
                "numerical_checks": "Refit only the three registered Phase 4 designs. Compare all coefficient estimates, both uncertainty summaries and fit metrics.",
                "relative_tolerance": REFERENCE_REL_TOL,
                "absolute_tolerance": REFERENCE_ABS_TOL,
                "runtime_policy": "The Phase 4 review preserves its original environment. The Phase 5 reference check compares numerical results at the declared tolerance rather than requiring identical environment metadata.",
            },
            "retention": "Keep the Phase 4 integrated fit on all 842 eligible participants as the reference for the final report.",
        },
        "input_support": support,
        "comparisons": {
            "broader_sample": {
                "baseline": "eligible_model_normal",
                "alternative": "eligible_model_all",
                "single_change": "Allow the other recorded cognitive-status categories under the same score, thickness and covariate eligibility rules.",
                "fit": "Use the reference outcome, units, covariates, centering and integrated additive OLS formula in both samples.",
                "report": "Counts by cognitive status, coefficients, individual 95% intervals and descriptive fit. For the broader fit, the common slope pools the recorded statuses without adjusting for status itself.",
                "interpretation": "The broader fit describes a different population mix. Differences may reflect between-status structure as well as the altered inclusion rule.",
            },
            "mri_completeness": {
                "same_target": "Normal cognition with usable ANIMALS, NACCAGE, SEX and EDUC.",
                "baseline_membership": "Derive UDS eligibility from the original UDS file before requiring a matching usable MRI value, using the Phase 1 rules.",
                "alternative_membership": "Require a matched row with usable left thickness in addition to the same UDS eligibility.",
                "comparable_model": "Compare UDS-only on both sets. Evaluate membership by ID internally and summarize n, membership differences and composition in aggregate. Use the integrated fit to describe the MRI-eligible set.",
                "current_support": "The known supplied extract has a usable left value for every matched participant and no UDS-only participants. The audit should verify whether the two eligible sets are identical.",
                "interpretation": "An identical eligible set makes the within-extract UDS-only comparison degenerate. Selection into the teaching extract and effects among people with missing MRI outside it remain unidentifiable.",
            },
            "visit_rule": {
                "availability": "Alternative visit selection is not evaluable on the supplied files.",
                "evidence": phase4["sampling_model"]["visit_and_scan_timing"],
                "action": "Document the one supplied pair per participant and the missing visit/scan metadata. Apply a versioned plan if authoritative repeated records and timing become available.",
            },
            "bilateral_thickness": {
                "sample": "Intersect the reference sample with usable right thickness. For both fits use the exact same ordered participants.",
                "baseline_column": LEFT_THICKNESS,
                "alternative_components": [LEFT_THICKNESS, RIGHT_THICKNESS],
                "alternative_definition": "Arithmetic mean of usable left and right inferior parietal thickness in mm, with each hemisphere weighted equally.",
                "alternative_coding": "(bilateral_mean_mm - 2.344) / 0.1",
                "single_change": "Replace left thickness with bilateral mean thickness. Hold ANIMALS, all UDS covariates and the OLS estimator fixed.",
                "interpretation": "Compare slopes per 0.1 mm in their respective anatomical measures. Use 2.344 mm as a centering constant in both fits and avoid comparing their intercepts as if they represented the same anatomical quantity.",
            },
            "influence": {
                "baseline": "Full integrated Phase 4 fit on the reference sample.",
                "screen": "Integrated model Cook distance strictly greater than 4/n, calculated once on the full design using p=5 including the intercept.",
                "order": "Descending full-fit Cook distance. Use ascending participant ID only to break exact ties and keep IDs in memory.",
                "maximum_joint_omissions": INFLUENCE_MAX,
                "comparison": "Temporarily omit the first min(5, number flagged) observations together and refit the same formula once. If no case passes the screen, report zero eligible omissions.",
                "historical_context": "The Phase 4 diagnostic already includes all flagged individual omissions and a larger 70-observation union omission.",
                "report": "Counts, anonymous diagnostic summaries, coefficient shifts and both usual-formula intervals. Discuss the limits of post-selection interval coverage. Keep these fits outside the primary inference.",
            },
            "bootstrap": {
                "target": "The same adjusted left-thickness slope on the same 842-person reference sample.",
                "draws": BOOTSTRAP_DRAWS,
                "seed": BOOTSTRAP_SEED,
                "generator": "numpy.random.Generator with numpy.random.PCG64",
                "resampling_unit": "Whole participant row, including ANIMALS, thickness and all UDS covariates.",
                "draw_size": support["reference_n"],
                "replacement": True,
                "row_order": "Use the registered ascending participant-ID order before drawing indices.",
                "fit": "Equal-weight OLS using the same fixed coding and a full-rank least-squares solve. Repeated indices represent bootstrap multiplicity.",
                "interval": {"method": "percentile", "confidence_level": 0.95, "quantiles": [0.025, 0.975],
                             "quantile_method": "linear"},
                "invalid_fit_policy": "Count non-finite or rank-deficient fits explicitly. If any occur, stop before publishing the interval and document the issue. Keep the planned draw sequence fixed for review.",
                "comparison": "Compare the percentile interval and empirical bootstrap SE with the Phase 4 model-based interval and SE. Include HC3 as existing reference context. Report shifts in endpoints and width.",
                "limits": "Row resampling assumes independent participants and represents variation within the modeled extract. It cannot add source-imputation, selection or visit-alignment uncertainty that is absent from the inputs.",
            },
            "cognitive_subgroup": {
                "groups": {"normal_cognition": 1, "mci": MCI_CODE},
                "sample": "Participants in eligible_model_all with the two listed cognitive codes. Use one supplied row per participant.",
                "rationale": "Evaluate whether the thickness slope differs across two clinically distinct cognitive groups with enough observations for this limited comparison.",
                "within_group_fit": "Use the reference integrated formula separately in each group, with the same original units and coding.",
                "direct_comparison": "Fit the pooled fully interacted specification: ANIMALS ~ G * (thickness_0_1mm + age_c + female + education_c), with normal cognition as the reference and G=1 for MCI. The G-by-thickness term equals the MCI slope minus the normal-cognition slope.",
                "uncertainty": "Report model-based and HC3 individual 95% intervals for the two slopes. Report the interaction estimate and an HC3 individual 95% interval using the pooled residual degrees of freedom as an approximate t reference.",
                "verification": "Require both subgroup designs and the pooled design to have full column rank. Confirm that group slopes from the fully interacted fit agree with their separate-model estimates.",
                "scope": "Exploratory comparison of associations across the supplied groups. The difference and interval are the evidence for the comparison. Avoid conclusions from separate significance patterns.",
            },
        },
        "interpretation_policy": {
            "comparison_scale": "Report effect estimates, absolute change when comparable, individual 95% intervals, sample sizes and model fit in original units.",
            "selection": "Use the fixed set of comparisons regardless of sign, interval inclusion of zero or nominal p-value. Comparisons sharing records are dependent.",
            "robustness": "Evaluate magnitude, direction and uncertainty together. A sign change near zero or non-rejection is insufficient for an equivalence claim. No practical-equivalence margin has been supplied.",
            "exploratory_outcome_candidate": "VEG remains an optional candidate outside this locked comparison set. It requires an expressly versioned decision about score eligibility and sample comparability before fitting.",
        },
        "literature": {"minimum_studies": 3, "step": 2,
                       "table_fields": ["data and population", "question", "statistical method", "key result", "comparability with the project"]},
        "requirement_map": {key: {"step": step, "planned_evidence": description}
                            for key, (step, description) in REQUIREMENTS.items()},
        "work_sequence": [
            {"step": 1, "task": "Validate the reference and freeze the Phase 5 specification"},
            {"step": 2, "task": "Review at least three relevant studies"},
            {"step": 3, "task": "Assess sample definitions, MRI completeness and visit metadata"},
            {"step": 4, "task": "Assess the bilateral alternative on matched support"},
            {"step": 5, "task": "Assess Cook-distance sensitivity and participant bootstrap intervals"},
            {"step": 6, "task": "Compare the two specified cognitive groups"},
            {"step": 7, "task": "Synthesize evidence and verify the final reports, notebook and figures"},
        ],
        "delivery": {
            "report_directory": "reports/final",
            "sensitivity_and_limitations": "Separate LaTeX source and compiled document inside the final report directory.",
            "code_directory": "src",
            "notebook": "notebooks/phase5.ipynb",
            "figure_directory": "figures/phase5",
            "tag": "phase-5-final",
            "repository_and_course_submission": "Git publication, tagging and CW packaging are manual.",
        },
        "privacy": "Save aggregate counts, statistics and checksums. Keep participant IDs, model rows and individual bootstrap indices local in memory.",
        "course_guide": {"path": "docs/amarpr.pdf", "sections": [9, 10, 11], "pages": [15, 16, 17, 18, 19]},
    }
    require(fingerprints == {name: fingerprint(Path(name)) for name in fingerprints},
            "A reference file or the planning code changed during validation.")
    return plan


def load_registered_plan() -> dict:
    """Reject a stale or modified plan before downstream sensitivity fits."""
    saved = read_json(REPORT_FILE)
    require(saved == build_plan(), "The Phase 5 plan differs from its sources. Review and version the plan before fitting.")
    return saved


def main() -> None:
    plan = build_plan()
    if REPORT_FILE.exists():
        require(read_json(REPORT_FILE) == plan,
                "An existing Phase 5 plan differs. Review and version it explicitly before replacing it.")
        return
    write_output(REPORT_FILE, (json.dumps(plan, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
