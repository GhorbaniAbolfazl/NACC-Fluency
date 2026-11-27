"""Audit the Phase 4 calculations and the public course deliverables."""

import json
from io import BytesIO
from pathlib import Path

import nbformat
import numpy as np
from scipy.stats import t

from . import phase4_comparison as comparison
from . import phase4_diagnostics as diagnostics
from . import phase4_simple as simple
from . import phase4_tables as tables
from .build_dataset import write_output
from .phase1_figures import GRAYSCALE, plt
from .phase4_diagnostic_figures import build_figures
from .phase4_plan import fingerprint, require
from .phase4_prediction import REPORT_FILE as PREDICTION_FILE

REPORT_FILE = Path("reports/phase4/regression_review.json")
NOTEBOOK = Path("notebooks/phase4.ipynb")
REPORT = Path("reports/phase4.pdf")
REPORT_SOURCE = Path("reports/phase4/phase4.tex")


def check_close(actual, expected, message):
    require(np.allclose(actual, expected, rtol=2e-10, atol=2e-11), message)


def check_regression(results: dict, context: dict) -> None:
    """Use a separate SVD solution to check the QR fits and both covariance methods."""
    for name, design in context["designs"].items():
        matrix = np.asarray(design.matrix)
        y = np.asarray(design.response)
        n, p = matrix.shape
        u, singular, vh = np.linalg.svd(matrix, full_matrices=False)
        inverse_map = (vh.T / singular) @ u.T
        beta = inverse_map @ y
        residual = y - matrix @ beta
        leverage = np.sum(u**2, axis=1)
        sse = float(residual @ residual)
        mse = sse / (n - p)
        covariance = mse * (inverse_map @ inverse_map.T)
        hc3 = (inverse_map * (residual / (1 - leverage))**2) @ inverse_map.T
        model = results["simple"] if name == "mri_only" else results["multiple"]["models"][name]
        fit = context["fits"][name]
        check_close(beta, fit.coefficients, f"Independent SVD coefficient mismatch: {name}.")
        check_close([sse, mse], [model["fit"]["SSE"], fit.residual_sd**2], f"Independent residual variance mismatch: {name}.")
        check_close(covariance, model["covariance"]["model_based"], f"Classical covariance mismatch: {name}.")
        check_close(hc3, model["covariance"]["hc3"], f"HC3 covariance mismatch: {name}.")
        check_close(leverage, fit.leverages, f"Leverage mismatch: {name}.")
        studentized = residual / np.sqrt(mse * (1 - leverage))
        cook = studentized**2 * leverage / (p * (1 - leverage))
        check_close(studentized, context["arrays"][name]["studentized"], f"Studentized residual mismatch: {name}.")
        check_close(cook, context["arrays"][name]["cooks_distance"], f"Cook distance mismatch: {name}.")
        if name == "integrated":
            profile = results["prediction"]["intervals"]
            x0 = np.asarray(design.reference_row)
            center = float(x0 @ beta)
            mean_variance = float(x0 @ covariance @ x0)
            critical = float(t.ppf(.975, n - p))
            check_close(center, profile["point_estimate"], "Profile center mismatch.")
            for label, variance in (("mean_confidence_interval", mean_variance),
                                    ("individual_prediction_interval", mean_variance + mse)):
                half_width = critical * np.sqrt(variance)
                bounds = profile[label]
                check_close([center - half_width, center + half_width], [bounds["lower"], bounds["upper"]],
                            f"Profile interval mismatch: {label}.")


def figure_groups(results: dict, context: dict):
    """Draw figures with their original settings and keep row-level arrays in memory."""
    style = {"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "figure.facecolor": "white"}
    with plt.rc_context({**style, "axes.labelcolor": GRAYSCALE.ink, "axes.labelpad": 9}):
        yield [(Path(results["simple"]["figure"]["path"]),
                simple.regression_figure(results["simple"], context["designs"]["mri_only"]))]
    with plt.rc_context(style):
        yield [(Path(results["comparison"]["figure"]["path"]), comparison.comparison_figure(results["comparison"]))]
    with plt.rc_context({**style, "font.size": 9}):
        paths = [Path(item["path"]) for item in results["diagnostic"]["figures"]]
        yield list(zip(paths, build_figures(results["diagnostic"], context)))


def check_figures(results: dict, context: dict) -> list[str]:
    checked = []
    for group in figure_groups(results, context):
        for path, figure in group:
            try:
                with BytesIO() as buffer:
                    figure.savefig(buffer, format="pdf",
                                   metadata={"Title": path.name, "CreationDate": None, "ModDate": None})
                    require(buffer.getvalue() == path.read_bytes(), f"Saved figure differs from the regenerated PDF: {path}")
            finally:
                plt.close(figure)
            checked.append(str(path))
    return checked


def audit_analysis() -> tuple[dict, dict, list[str]]:
    """Verify the analytical chain without depending on assembled deliverables."""
    results = tables.load_results()
    diagnostic, context = diagnostics.analyze()
    require(diagnostic == results["diagnostic"], "Step 6 changed during review.")
    check_regression(results, context)
    require(results["prediction"]["primary_test"] == results["comparison"]["primary_test"]
            == results["multiple"]["primary_test"], "Primary test changed across steps.")
    require(set(results["comparison"]["collinearity"]) == {"uds_only", "integrated"},
            "Both multiple-model collinearity summaries are required.")
    paths = check_figures(results, context)
    require(len(paths) == 5 and len(set(paths)) == 5, "Expected five distinct figures.")
    return results, context, paths


def requirement_coverage(plan: dict) -> dict:
    sections = {
        "general_1": ("design",), "general_2": ("simple",), "general_3": ("simple",),
        "general_4": ("multiple",), "general_5": ("compare",), "general_6": ("design", "multiple"),
        "general_7": ("diagnostics",), "general_8": ("profile",), "general_9": ("profile",),
        "general_10": ("simple",), "general_11": ("compare",), "track_c1": ("design",),
        "track_c2": ("design", "multiple"), "track_c3": ("simple", "diagnostics"),
        "track_c4": ("multiple", "compare"), "track_c5": ("design",), "track_c6": ("diagnostics",),
    }
    require(set(sections) == set(plan["requirement_map"]) - {"status"}, "Requirement coverage differs from the frozen plan.")
    coverage = {key: {"status": "addressed", "sections": list(labels), "scope": plan["requirement_map"][key]}
                for key, labels in sections.items()}
    coverage["track_c5"]["status"] = "partially_verifiable"
    coverage["track_c5"]["limitation"] = plan["sampling_model"]["visit_and_scan_timing"]["track_c5_status"]
    return coverage


def build_review() -> dict:
    results, _, figures = audit_analysis()
    require(tables.TABLE_FILE.read_bytes() == tables.render_tables(results), "The shared numeric tables are stale.")
    document = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(document)
    code = [cell for cell in document.cells if cell.cell_type == "code"]
    require(bool(code) and all(cell.execution_count is not None for cell in code), "Execute the public notebook.")
    outputs = [output for cell in code for output in cell.outputs]
    require(all(output.output_type != "error" for output in outputs), "The public notebook contains an error.")
    require(sum("image/png" in output.get("data", {}) for output in outputs) == 5, "The public notebook must display all five figures.")
    coverage = requirement_coverage(results["plan"])
    source = REPORT_SOURCE.read_text(encoding="utf-8")
    for sections in (item["sections"] for item in coverage.values()):
        require(all(r"\label{sec:" + label + "}" in source for label in sections), "A report evidence section is missing.")
    require(all(Path(name).name in source for name in figures), "The report source must include all five figures.")
    pdf = REPORT.read_bytes()
    require(pdf.startswith(b"%PDF-") and pdf.rstrip().endswith(b"%%EOF"), "The compiled public report is missing or incomplete.")
    files = {**results["prediction"]["source_files_sha256"]}
    for path in (PREDICTION_FILE, tables.TABLE_FILE, NOTEBOOK, REPORT, REPORT_SOURCE, Path("reports/template.tex"),
                 Path("src/phase4_tables.py"), Path("src/phase4_review.py")):
        files[str(path)] = fingerprint(path)
    require(files == {name: fingerprint(Path(name)) for name in files}, "Files changed during final review.")
    return {
        "stage": "Phase 4, Step 8", "sample": results["prediction"]["sample"],
        "participants": results["prediction"]["participants"],
        "analysis_validation": [
            "Steps 1-7 regenerate and match saved aggregate artifacts.",
            "An independent SVD checks all three QR fits, model-based and HC3 covariances, and influence coordinates.",
            "The same SVD independently reproduces both fixed-profile intervals.",
            "All five vector PDFs reproduce byte for byte.",
            "The primary test agrees across regression, comparison and prediction.",
        ],
        "delivery_validation": [
            "Shared LaTeX numeric tables reproduce from the validated analysis.",
            "The public notebook is valid, executed without recorded errors and contains five figure outputs.",
            "The public LaTeX source includes all five figures and the named evidence sections.",
            "The compiled public PDF has the expected PDF header and end marker.",
        ],
        "coverage": coverage,
        "source_files_sha256": files,
        "delivery": {"report": str(REPORT), "notebook": str(NOTEBOOK),
                     "manual": "Repository publication, phase-4 tagging and course packaging are manual."},
    }


def main() -> None:
    content = (json.dumps(build_review(), indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != content:
        write_output(REPORT_FILE, content)


if __name__ == "__main__":
    main()
