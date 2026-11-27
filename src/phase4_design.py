"""Prepare and validate the three common-sample Phase 4 designs without fitting."""

import json
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from hashlib import sha256
from importlib.metadata import version
from math import isfinite
from pathlib import Path
from platform import python_version
from statistics import fmean, stdev
from sys import float_info

from scipy.linalg import qr, svd

from .build_dataset import LEFT_THICKNESS, write_output
from .data_audit import PARTICIPANT_ID
from .phase4_plan import (
    REPORT_FILE as PLAN_FILE, fingerprint, load_registered_plan, require, validate_inputs,
)

REPORT_FILE = Path("reports/phase4/design_matrices.json")
FEATURE_ORDER = ("intercept", "thickness_0_1mm", "age_c", "female", "education_c")


@dataclass(frozen=True)
class ModelDesign:
    """Immutable participant-aligned values, retained in memory rather than exported."""

    name: str
    participant_ids: tuple[str, ...]
    columns: tuple[str, ...]
    matrix: tuple[tuple[float, ...], ...]
    response: tuple[float, ...]
    reference_row: tuple[float, ...]


def canonical_hash(values: object) -> str:
    """Hash compact UTF-8 JSON with explicit float values and no nonfinite numbers."""
    text = json.dumps(values, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    return sha256(text.encode("utf-8")).hexdigest()


def encode_predictors(row: dict, plan: dict) -> tuple[float, ...]:
    """Apply the registered coding directly; never evaluate an expression string."""
    spec = plan["variables"]["explanatory"]
    require(row["SEX"] in (1, 2), "Recorded sex must be 1 or 2 before indicator coding.")
    references = {
        "age_c": ("NACCAGE", spec["age_c"]["reference_years"], Decimal(1)),
        "education_c": ("EDUC", spec["education_c"]["reference_years"], Decimal(1)),
        "thickness_0_1mm": (LEFT_THICKNESS, spec["thickness_0_1mm"]["reference_mm"], Decimal("0.1")),
    }
    values = {"intercept": 1.0, "female": float(spec["female"]["code_map"][str(row["SEX"])])}
    for term, (column, reference, unit) in references.items():
        raw = Decimal(str(row[column]))
        require(raw.is_finite(), f"A nonfinite value cannot enter {term}.")
        values[term] = float((raw - Decimal(str(reference))) / unit)
    require(all(isfinite(value) for value in values.values()), "Every encoded predictor must be finite.")
    return tuple(values[name] for name in FEATURE_ORDER)


def prepare_designs() -> tuple[dict, dict[str, ModelDesign]]:
    """Use the frozen plan and identical participant order for every model."""
    plan = load_registered_plan()
    rows, flow, inputs = validate_inputs()
    require(inputs == plan["inputs"]["source_files_sha256"], "Inputs changed after the registered plan was loaded.")
    require(flow["samples"]["model_normal"]["included"] == len(rows) == 842,
            "The common modeling sample must contain the registered 842 participants.")
    ids = tuple(row[PARTICIPANT_ID] for row in rows)
    require(ids == tuple(sorted(ids)) and len(set(ids)) == len(ids),
            "The common participant order must be unique and ascending by ID.")
    membership = sha256(("\n".join(ids) + "\n").encode("utf-8")).hexdigest()
    require(membership == plan["sample"]["membership_sha256"], "Model membership differs from the registered plan.")
    response = tuple(float(row["ANIMALS"]) for row in rows)
    require(all(isfinite(value) and value.is_integer() and 0 <= value <= 77 for value in response),
            "Only registered usable integer animal-naming scores may enter the response vector.")
    require(len(set(response)) > 1, "The response must vary before fitting or calculating R-squared.")
    encoded = tuple(encode_predictors(row, plan) for row in rows)
    reference = encode_predictors(plan["prediction"]["profile"]["values"], plan)
    require(reference == tuple(float(plan["prediction"]["design_row"][name]) for name in FEATURE_ORDER),
            "The recorded reference profile does not match its registered design row.")
    designs = {}
    for name, specification in plan["models"].items():
        columns = tuple(specification["columns"])
        require(len(columns) == len(set(columns)) and columns[0] == "intercept"
                and set(columns) <= set(FEATURE_ORDER), "Unexpected or repeated design columns.")
        indices = tuple(FEATURE_ORDER.index(column) for column in columns)
        matrix = tuple(tuple(row[index] for index in indices) for row in encoded)
        designs[name] = ModelDesign(name, ids, columns, matrix, response,
                                    tuple(reference[index] for index in indices))
    # Projection back to the supplied units checks the sign, centering and 0.1-mm scale.
    by_term = {name: index for index, name in enumerate(FEATURE_ORDER)}
    spec = plan["variables"]["explanatory"]
    for raw, values in zip(rows, encoded):
        for term, column, scale, ref in (
            ("age_c", "NACCAGE", 1.0, spec["age_c"]["reference_years"]),
            ("education_c", "EDUC", 1.0, spec["education_c"]["reference_years"]),
            ("thickness_0_1mm", LEFT_THICKNESS, 0.1, spec["thickness_0_1mm"]["reference_mm"]),
        ):
            recovered = values[by_term[term]] * scale + ref
            tolerance = 8 * float_info.epsilon * max(1.0, abs(raw[column]))
            require(abs(recovered - raw[column]) <= tolerance, f"Coding does not recover the supplied {column} value.")
        require(int(values[by_term["female"]]) == raw["SEX"] - 1, "Recorded-sex indicator direction is reversed.")
    require(inputs == {path: fingerprint(Path(path)) for path in inputs}, "Inputs changed during design preparation.")
    return plan, designs


def matrix_diagnostics(design: ModelDesign) -> dict:
    """Check dimensions, finite values, SVD/QR rank and predictor-only conditioning."""
    n, p = len(design.participant_ids), len(design.columns)
    require(n > p >= 2 and len(design.matrix) == len(design.response) == n,
            "Design dimensions and participant/response lengths must agree, with positive residual df.")
    require(len(set(design.participant_ids)) == n and all(design.participant_ids), "Each design requires unique nonblank IDs.")
    require(len(set(design.columns)) == p and design.columns[0] == "intercept", "Design columns must be unique with an explicit intercept first.")
    require(all(len(row) == p and all(isfinite(value) for value in row) for row in design.matrix),
            "Every design row must have the expected finite columns.")
    require(len(design.reference_row) == p and all(isfinite(value) for value in design.reference_row),
            "The profile row must match the model columns.")
    require(all(isfinite(value) and float(value).is_integer() and 0 <= value <= 77 for value in design.response),
            "Response values must remain finite usable integer counts.")
    require(all(row[0] == 1.0 for row in design.matrix), "The intercept column must contain ones.")
    columns = list(zip(*design.matrix))
    require(all(len(set(values)) > 1 for values in columns[1:]), "A non-intercept predictor is constant.")

    def factorize(matrix: tuple[tuple[float, ...], ...]) -> dict:
        singular = [float(value) for value in svd(matrix, full_matrices=False, compute_uv=False,
                                                 check_finite=True, lapack_driver="gesvd")]
        tolerance = max(n, p) * float_info.epsilon * singular[0]
        svd_rank = sum(value > tolerance for value in singular)
        _, r, pivots = qr(matrix, mode="economic", pivoting=True, check_finite=True)
        diagonal = [abs(float(r[index, index])) for index in range(p)]
        qr_rank = sum(value > tolerance for value in diagonal)
        require(svd_rank == qr_rank == p, f"{design.name} is rank deficient or its numerical rank is ambiguous; no column is silently removed.")
        return {"svd_rank": svd_rank, "pivoted_qr_rank": qr_rank, "singular_values": singular,
                "rank_tolerance": tolerance, "qr_absolute_diagonal": diagonal,
                "qr_pivot_order": [design.columns[int(index)] for index in pivots],
                "condition_number_2norm": singular[0] / singular[-1]}

    raw = factorize(design.matrix)
    centers = [0.0, *(fmean(values) for values in columns[1:])]
    scales = [1.0, *(stdev(values) for values in columns[1:])]
    require(all(scale > 0 and isfinite(scale) for scale in scales), "Predictor standardization requires positive sample SDs.")
    standardized = tuple(tuple((value - centers[j]) / scales[j] for j, value in enumerate(row))
                         for row in design.matrix)
    standardized_result = factorize(standardized)
    return {
        "rows": n, "columns": p, "full_column_rank": True, "residual_df": n - p,
        "raw_design": raw, "standardized_design": standardized_result,
        "standardization": "For conditioning diagnostics only: retain intercept=1; center each other column at its mean and divide by its sample SD (n-1). The model matrix and its registered coefficient units are unchanged.",
        "standardization_centers": dict(zip(design.columns, centers)),
        "standardization_divisors": dict(zip(design.columns, scales)),
        "column_summaries": {name: {"minimum": min(values), "maximum": max(values),
                                     "mean": fmean(values), "sample_sd": stdev(values)}
                             for name, values in zip(design.columns, columns)},
    }


def _build_report(plan: dict, designs: dict[str, ModelDesign]) -> dict:
    """Publish aggregate checks and hashes, never participant-level matrices or scores."""
    tracked = {**plan["inputs"]["source_files_sha256"], str(PLAN_FILE): fingerprint(PLAN_FILE),
               "src/phase4_design.py": fingerprint(Path("src/phase4_design.py")),
               "requirements.txt": fingerprint(Path("requirements.txt"))}
    require(set(designs) == set(plan["models"]), "The three prepared models must match the registered model set.")
    first = next(iter(designs.values()))
    expected_n = plan["sample"]["counts"]["model_normal"]
    require(len(first.participant_ids) == expected_n, "Design sample size differs from the plan.")
    reports = {}
    for name, design in designs.items():
        spec = plan["models"][name]
        require(design.participant_ids == first.participant_ids and design.response == first.response,
                "All models must use identical IDs, response values and row order.")
        require(design.name == name and design.columns == tuple(spec["columns"]), "Prepared model columns differ from the registered order.")
        diagnostics = matrix_diagnostics(design)
        require(diagnostics["residual_df"] == spec["residual_df_if_full_rank"], "Residual degrees of freedom differ from the registered model.")
        require(design.reference_row == tuple(float(plan["prediction"]["design_row"][column]) for column in design.columns),
                "A model's reference row differs from the registered coding.")
        reports[name] = {
            "role": spec["role"], "column_order": list(design.columns),
            "matrix_sha256": canonical_hash(design.matrix), "response_sha256": canonical_hash(design.response),
            "ordered_ids_sha256": sha256(("\n".join(design.participant_ids) + "\n").encode("utf-8")).hexdigest(),
            "reference_row": list(design.reference_row),
            "rows_equal_to_reference": sum(row == design.reference_row for row in design.matrix),
            **diagnostics,
        }
    integrated = designs["integrated"]
    for design in designs.values():
        indices = [integrated.columns.index(column) for column in design.columns]
        require(design.matrix == tuple(tuple(row[j] for j in indices) for row in integrated.matrix),
                "A model column differs from the corresponding integrated-model column.")
    require(all(report["ordered_ids_sha256"] == plan["sample"]["membership_sha256"] for report in reports.values()),
            "Prepared participant membership differs from the registered fingerprint.")
    require(reports["integrated"]["rows_equal_to_reference"] >= 1,
            "The reference profile must still be jointly observed in the common sample.")
    female = [row[integrated.columns.index("female")] for row in integrated.matrix]
    counts = Counter(female)
    require(set(counts) == {0.0, 1.0} and counts[0.0] == plan["sample"]["sex_counts"]["male_code_1"]
            and counts[1.0] == plan["sample"]["sex_counts"]["female_code_2"], "The sex indicator counts differ from the registered source codes.")
    require(tracked == {path: fingerprint(Path(path)) for path in tracked}, "Inputs changed while the design report was being built.")
    return {
        "stage": "Phase 4, Step 2", "status": "Common-sample model matrices prepared and validated; no regression fitted.",
        "sample": "model_normal", "participants": expected_n, "source_files_sha256": tracked,
        "environment": {"python": python_version(), "scipy": version("scipy"), "numpy_runtime": version("numpy"),
                        "numeric_type": "IEEE-754 binary64 Python floats passed to scipy.linalg", "machine_epsilon": float_info.epsilon},
        "coding": {"method": "Registered reference values with Decimal subtraction and 0.1-mm scaling before float conversion; no evaluation of expression strings.",
                   "reference_profile": plan["prediction"]["profile"]["values"],
                   "sex_indicator_counts": {"male_0": counts[0.0], "female_1": counts[1.0]},
                   "round_trip": "All age, education and thickness values recover within 8*epsilon*max(1,abs(source value)); sex recovers exactly.",
                   "intercepts": "Every model has a leading column of ones.",
                   "thickness": "Continuous centered thickness per 0.1 mm; the Phase 3 cutoff is not used."},
        "response": {"column": "ANIMALS", "unit": "Correct responses in 60 seconds", "participants": len(first.response),
                     "all_finite_usable_integers": True, "minimum": min(first.response), "maximum": max(first.response),
                     "same_values_and_order_for_every_model": True, "response_sha256": canonical_hash(first.response)},
        "models": reports,
        "validation": {"passed": True, "checks": [
            "Frozen plan and completed Phase 3 input fingerprints are current.",
            "Participant table and eligibility reproduce from the raw CSVs in memory.",
            "The same 842 unique sorted participant IDs and response order enter all three models.",
            "Registered centering, sex coding and 0.1-mm thickness units are preserved and invertible within floating-point tolerance.",
            "Each model has exactly the registered columns, order and intercept.",
            "All design values are finite, each explanatory column varies, and usable response values are integer counts.",
            "Raw and standardized matrices have full column rank by SVD and pivoted QR with a recorded tolerance.",
            "Residual degrees of freedom agree with n minus the registered parameter count.",
            "Common columns are exactly identical across the three model matrices.",
            "The frozen prediction profile maps to the registered row and is jointly observed in the integrated sample.",
            "Inputs and plan remained unchanged during preparation."],
            "rank_method": "SVD with scipy.linalg.svd(compute_uv=False, lapack_driver='gesvd') cross-checked by economic pivoted scipy.linalg.qr; tolerance=max(n,p)*epsilon*largest singular value.",
            "conditioning_limit": "Condition numbers describe the predictor matrix, not validity of regression assumptions or predictive performance. VIF, residual checks and influence analyses belong to later steps."},
        "privacy": {"persisted": "Aggregate summaries, model column names, reference rows, singular values and fingerprints only.",
                    "not_exported": ["Participant IDs", "Participant-level design rows", "Participant-level response values"],
                    "in_memory_api": "prepare_designs() returns immutable ModelDesign objects; later steps load them with load_prepared_designs()."},
        "hash_encoding": "Matrix and response hashes use compact UTF-8 JSON, ensure_ascii=True, allow_nan=False, separators=(',',':'), Python float values and retained order. IDs use the registered sorted newline encoding, not compact JSON.",
        "scope": {"not_performed": ["Regression fitting", "Score-thickness correlation", "Coefficient estimates or tests", "Residual or VIF diagnostics", "Predictions or intervals"],
                  "limitations": "Numerical readiness does not verify participant independence, sampling design, source imputation, scan quality or clinical/MRI temporal proximity. The Step 1 limitations remain in force.",
                  "next_step": "Fit the prespecified MRI-only simple regression and calculate Pearson correlation on these same 842 participants in Step 3."},
        "sources": {"svd": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.svd.html",
                    "pivoted_qr": "https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html"},
    }


def build_report() -> dict:
    plan, designs = prepare_designs()
    return _build_report(plan, designs)


def load_prepared_designs() -> tuple[dict, dict[str, ModelDesign]]:
    """Require a reproducible saved preparation report before any model fitting."""
    saved = json.loads(REPORT_FILE.read_text(encoding="utf-8"))
    plan, designs = prepare_designs()
    require(saved == _build_report(plan, designs),
            "The saved Phase 4 design report is stale. Review changes and rebuild it before fitting.")
    return plan, designs


def main() -> None:
    output = (json.dumps(build_report(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if REPORT_FILE.exists() and REPORT_FILE.read_bytes() == output:
        return
    write_output(REPORT_FILE, output)


if __name__ == "__main__":
    main()
