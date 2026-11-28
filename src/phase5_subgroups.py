"""Compare the adjusted thickness slopes in normal cognition and MCI."""

import json
from hashlib import sha256
from math import isclose
from pathlib import Path

from .build_dataset import STATUS, build_dataset, write_output
from .data_audit import PARTICIPANT_ID
from .phase4_design import FEATURE_ORDER, ModelDesign, canonical_hash, encode_predictors
from .phase4_plan import REPORT_FILE as PHASE4_PLAN, fingerprint, require
from .phase4_regression import coefficient_table, fit_ols, fit_summary
from .phase5_plan import (
    PRIMARY_FILE, REFERENCE_ABS_TOL, REFERENCE_REL_TOL,
    REPORT_FILE as PHASE5_PLAN, load_registered_plan, read_json,
)

REPORT_FILE = Path("reports/final/cognitive_subgroups.json")
THICKNESS_TERM = "thickness_0_1mm"
GROUP_TERM = "mci"
INTERACTION = f"{GROUP_TERM}_x_{THICKNESS_TERM}"


def close(actual: float, expected: float) -> bool:
    return isclose(actual, expected, rel_tol=REFERENCE_REL_TOL, abs_tol=REFERENCE_ABS_TOL)


def build_report() -> dict:
    sources = (PHASE5_PLAN, PHASE4_PLAN, PRIMARY_FILE, Path("src/phase5_subgroups.py"))
    hashes = {str(path): fingerprint(path) for path in sources}
    plan = load_registered_plan()
    phase4 = read_json(PHASE4_PLAN)
    specification = plan["comparisons"]["cognitive_subgroup"]
    codes = specification["groups"]
    expected_sizes = plan["input_support"]["normal_cognition_vs_mci_n"]
    reference = read_json(PRIMARY_FILE)["models"]["integrated"]
    columns = tuple(plan["reference"]["column_order"])
    require(columns == FEATURE_ORDER and codes["normal_cognition"] != codes["mci"],
            "The registered integrated columns or group specification changed.")

    rows, _ = build_dataset(Path("data"))
    eligible = [row for row in rows if row["eligible_model_all"] and row[STATUS] in codes.values()]
    ids = tuple(row[PARTICIPANT_ID] for row in eligible)
    require(ids == tuple(sorted(set(ids))), "Each participant must enter once in ascending ID order.")
    matrix = tuple(encode_predictors(row, phase4) for row in eligible)
    response = tuple(float(row["ANIMALS"]) for row in eligible)
    profile = tuple(reference["reference_row"])
    confidence = phase4["estimation"]["confidence_level"]

    fits, groups = {}, {}
    for name, code in codes.items():
        indices = tuple(i for i, row in enumerate(eligible) if row[STATUS] == code)
        require(len(indices) == expected_sizes[name], "A subgroup size differs from the registered plan.")
        design = ModelDesign(name, tuple(ids[i] for i in indices), columns,
                             tuple(matrix[i] for i in indices), tuple(response[i] for i in indices), profile)
        membership = sha256(("\n".join(design.participant_ids) + "\n").encode("utf-8")).hexdigest()
        if name == "normal_cognition":
            require(membership == plan["reference"]["membership_sha256"]
                    and canonical_hash(design.matrix) == reference["matrix_sha256"]
                    and canonical_hash(design.response) == reference["response_sha256"],
                    "The normal-cognition subgroup differs from the frozen reference.")
        fit = fit_ols(design)
        fits[name] = fit
        groups[name] = {"cognitive_status_code": code, "participants": fit.participants,
                        "membership_sha256": membership, "coefficients": coefficient_table(fit, confidence),
                        "fit": fit_summary(fit)}

    # Multiplying the fixed predictors by G allows every slope and the intercept to vary by group.
    pooled_columns = columns + (GROUP_TERM,) + tuple(f"{GROUP_TERM}_x_{term}" for term in columns[1:])
    interaction_rows = []
    for row, predictors in zip(eligible, matrix, strict=True):
        g = float(row[STATUS] == codes["mci"])
        interaction_rows.append(predictors + (g,) + tuple(g * x for x in predictors[1:]))
    pooled = fit_ols(ModelDesign("fully_interacted", ids, pooled_columns, tuple(interaction_rows),
                                 response, profile + (0.0,) * len(columns)))
    p = len(columns)
    normal, mci = fits["normal_cognition"], fits["mci"]
    require(all(close(value, expected) for value, expected in zip(pooled.coefficients[:p],
                                                                   normal.coefficients, strict=True)),
            "The pooled model does not reproduce the normal-cognition fit.")
    require(all(close(main + interaction, expected)
                for main, interaction, expected in zip(pooled.coefficients[:p], pooled.coefficients[p:],
                                                       mci.coefficients, strict=True)),
            "The pooled model does not reproduce the MCI fit.")
    require(close(pooled.sse, normal.sse + mci.sse), "The pooled residual sum of squares does not partition by group.")
    j = columns.index(THICKNESS_TERM)
    contrast = coefficient_table(pooled, confidence)[INTERACTION]
    require(close(contrast["estimate"], mci.coefficients[j] - normal.coefficients[j])
            and close(pooled.hc3_covariance[p + j][p + j],
                      normal.hc3_covariance[j][j] + mci.hc3_covariance[j][j]),
            "The interaction disagrees with the separate slopes or their HC3 variance sum.")

    report = {
        "stage": "Phase 5, Step 6", "course_requirements": ["general_6", "track_c4"],
        "source_files_sha256": hashes,
        "specification": {
            "within_group_formula": plan["reference"]["model"],
            "pooled_formula": "ANIMALS ~ G * (thickness_0_1mm + age_c + female + education_c)",
            "group_indicator": {"0": "normal cognition", "1": "MCI"},
            "coding": plan["reference"]["coding"],
            "sample_rule": specification["sample"],
        },
        "groups": groups,
        "direct_comparison": {
            "estimate_meaning": "MCI minus normal cognition in correct animal responses per 0.1 mm of left inferior parietal thickness, with the same UDS predictors held fixed within each group.",
            "term": INTERACTION, "estimate": contrast["estimate"], "hc3": contrast["hc3"],
            "pooled_column_order": list(pooled.columns), "pooled_fit": fit_summary(pooled),
            "reference_distribution": "Approximate t reference with the fully interacted model's residual degrees of freedom.",
        },
        "validation": {
            "both_subgroup_designs_and_pooled_design_full_rank": True,
            "pooled_coefficients_reproduce_separate_fits": True,
            "pooled_SSE_equals_sum_of_subgroup_SSE": True,
            "interaction_HC3_variance_equals_sum_of_subgroup_HC3_variances": True,
        },
        "interpretation": {
            "scope": specification["scope"],
            "limits": "Use the direct interaction interval to describe the uncertainty of the difference. Inference is nominal and conditional on the supplied extract and fitted linear models. No equivalence margin is specified.",
        },
        "privacy": "Only aggregate summaries and hashes are persisted.",
    }
    require(hashes == {str(path): fingerprint(path) for path in sources},
            "An input changed during the subgroup comparison.")
    return report


def main() -> None:
    data = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != data:
        write_output(REPORT_FILE, data)


if __name__ == "__main__":
    main()
