"""Verify the final submission files against the reproduced Phase 5 synthesis."""

import json
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import nbformat

from .build_dataset import write_output
from .phase1_figures import plt, save_pdf
from .phase4_plan import fingerprint, require
from .phase5_plan import read_json
from .phase5_summary import (
    DIRECTORY, FIGURE_DIR, REPORT_FILE as SYNTHESIS, TABLE_FILE, build_figures,
    build_synthesis, compare, render_tables,
)

REPORT_FILE = DIRECTORY / "final_review.json"
NOTEBOOK = Path("notebooks/phase5.ipynb")
DOCUMENTS = tuple(DIRECTORY / name for name in ("final_report", "sensitivity_and_limitations"))
INCLUDE = re.compile(r"\\input\{([^}]+)\}")


def latex_sources(path: Path, sources: set[Path]) -> str:
    """Expand local includes for coverage and dependency checks."""
    require(path.is_file(), f"Missing LaTeX source {path}.")
    require(path not in sources, f"Repeated or cyclic LaTeX input {path}.")
    sources.add(path)
    return INCLUDE.sub(lambda match: latex_sources(Path(match[1]), sources), path.read_text(encoding="utf-8"))


def build_review() -> dict:
    result = build_synthesis()
    compare(read_json(SYNTHESIS), result, "final synthesis")
    require(TABLE_FILE.read_bytes() == render_tables(result), "The final LaTeX tables are out of date.")
    tracked = {SYNTHESIS, TABLE_FILE, NOTEBOOK, Path("src/phase5_review.py"), Path("src/phase5_summary.py")}
    figure_paths = []
    with TemporaryDirectory() as directory:
        for filename, figure in build_figures(result):
            try:
                save_pdf(figure, filename, Path(directory))
            finally:
                plt.close(figure)
            target = FIGURE_DIR / filename
            require(target.read_bytes() == (Path(directory) / filename).read_bytes(),
                    f"The final vector figure does not reproduce: {target}.")
            figure_paths.append(target)
    tracked.update(figure_paths)

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(notebook)
    code = [cell for cell in notebook.cells if cell.cell_type == "code"]
    require(code and all(cell.execution_count is not None for cell in code), "Execute the final notebook.")
    outputs = [output for cell in code for output in cell.outputs]
    require(all(output.output_type != "error" for output in outputs), "The notebook contains a cell error.")
    require(sum("image/png" in output.get("data", {}) for output in outputs) == len(figure_paths),
            "The notebook does not display both final figures.")
    require(any("text/markdown" in output.get("data", {}) for output in outputs),
            "The executed notebook is missing its result tables.")

    required_sections = (
        "Abstract", "Question and study design", "Data and descriptive statistics",
        "Probability and distributions", "Statistical inference", "Regression and correlation",
        "Sensitivity and robustness", "Discussion and limitations", "Conclusion", "References and appendix",
    )
    for stem in DOCUMENTS:
        tex, pdf = stem.with_suffix(".tex"), stem.with_suffix(".pdf")
        dependencies = set()
        expanded = latex_sources(tex, dependencies)
        needed = required_sections if stem.name == "final_report" else required_sections[6:8]
        require(all(r"\section{" + section + "}" in expanded for section in needed),
                f"A required section is missing from {tex}.")
        if stem.name == "final_report":
            require(all(str(path) in expanded for path in figure_paths), "The report does not include both final figures.")
        graphics = set(Path(match) for match in re.findall(
            r"\\finalfigure\{[^}]+\}\{([^}]+)\}", expanded))
        require(graphics and all(path.is_file() for path in graphics), f"A report figure is missing in {tex}.")
        content = pdf.read_bytes()
        require(content.startswith(b"%PDF-") and content.rstrip().endswith(b"%%EOF"),
                f"The compiled report is missing or incomplete: {pdf}.")
        tracked.update(dependencies | graphics | {pdf})
    return {
        "stage": "Phase 5, Step 7",
        "participants_in_primary_analysis": result["primary"]["n"],
        "coverage": result["coverage"],
        "visit_selection": "Applicability audited. Repeated records and timing fields are unavailable.",
        "verification": result["checks"] + [
            "Shared LaTeX tables regenerate from the validated synthesis.",
            "Both new grayscale vector figures reproduce within the project environment.",
            "The final notebook is valid, with every code cell executed, aggregate tables and two figure outputs.",
            "The two compiled documents are present, with all required sections and resolved source/figure dependencies.",
        ],
        "files_sha256": {str(path): fingerprint(path) for path in sorted(tracked)},
        "course_submission": {
            "tag": "phase-5-final",
            "manual_steps": "Publish the repository, apply the tag and package the course submission manually.",
            "local_data": "Keep supplied CSVs and the participant-level derived table local.",
        },
    }


def main() -> None:
    report = build_review()
    content = (json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()
    if not REPORT_FILE.exists() or REPORT_FILE.read_bytes() != content:
        write_output(REPORT_FILE, content)


if __name__ == "__main__":
    main()
