#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 02 — Implementation trace.

Purpose
-------
Trace how the frozen v2 implementation propagates through the scientific
pipeline without modifying any v2 files.

This audit:

1. Locates definitions of PsiCDM background functions.
2. Searches the frozen source tree for references to:
       H_psicdm
       E_psicdm
       dL_psicdm
       DM_psicdm
       DH_psicdm
3. Searches for the published v2 functional form:
       eps0 / (1 + z)
4. Searches for alternative PsiCDM implementations.
5. Inspects which analysis scripts import/call the frozen implementation.
6. Produces a machine-readable textual trace report.

IMPORTANT:
This script performs source-code tracing only.
It does NOT modify the frozen repository and does NOT decide whether
individual likelihood implementations are scientifically correct.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

AUDIT_ROOT = Path(__file__).resolve().parent
V2_ROOT = AUDIT_ROOT.parent / "psi-continuum-v2"
PACKAGE_ROOT = V2_ROOT / "psi_continuum_v2"
RESULTS_DIR = AUDIT_ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if not V2_ROOT.exists():
    raise RuntimeError(
        f"Frozen v2 repository not found:\n{V2_ROOT}"
    )

if not PACKAGE_ROOT.exists():
    raise RuntimeError(
        f"Frozen v2 package not found:\n{PACKAGE_ROOT}"
    )


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

TARGET_FUNCTIONS = [
    "E_psicdm",
    "H_psicdm",
    "dL_psicdm",
    "DM_psicdm",
    "DH_psicdm",
]

MODEL_TERMS = [
    "PsiCDMParams",
    "eps0",
]

# Direct textual patterns that may indicate the published v2 equation.
PUBLISHED_PATTERNS = [
    r"eps0\s*/\s*\(\s*1(?:\.0)?\s*\+\s*z\s*\)",
    r"eps0\s*/\s*\(\s*1\s*\+\s*z\s*\)",
    r"1(?:\.0)?\s*\+\s*eps0\s*/\s*\(\s*1(?:\.0)?\s*\+\s*z\s*\)",
]

# Files/directories that are not scientific source.
IGNORE_PARTS = {
    "__pycache__",
    ".git",
    ".pytest_cache",
}


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def python_files():
    """Return all relevant Python source files in the frozen package."""
    files = []

    for path in PACKAGE_ROOT.rglob("*.py"):
        if any(part in IGNORE_PARTS for part in path.parts):
            continue
        files.append(path)

    return sorted(files)


def rel(path: Path) -> str:
    """Path relative to frozen repository."""
    try:
        return str(path.relative_to(V2_ROOT))
    except ValueError:
        return str(path)


def read_source(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def line_number(source: str, index: int) -> int:
    return source.count("\n", 0, index) + 1


def find_text_occurrences(path: Path, term: str):
    source = read_source(path)
    results = []

    for match in re.finditer(re.escape(term), source):
        results.append(
            (
                line_number(source, match.start()),
                source.splitlines()[
                    line_number(source, match.start()) - 1
                ].strip(),
            )
        )

    return results


def get_function_definitions(path: Path):
    """Parse top-level and nested function definitions via AST."""
    source = read_source(path)

    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    definitions = []

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            definitions.append(
                (
                    node.name,
                    node.lineno,
                )
            )

    return definitions


def get_imports_and_calls(path: Path):
    """
    Return imported names and function call names using AST.

    This is used only for tracing. It is not a full static call graph.
    """
    source = read_source(path)

    try:
        tree = ast.parse(source)
    except SyntaxError:
        return [], []

    imports = []
    calls = []

    for node in ast.walk(tree):

        if isinstance(node, ast.ImportFrom):
            module = node.module or ""

            for alias in node.names:
                imports.append(
                    (
                        node.lineno,
                        module,
                        alias.name,
                        alias.asname,
                    )
                )

        elif isinstance(node, ast.Import):
            for alias in node.names:
                imports.append(
                    (
                        node.lineno,
                        alias.name,
                        None,
                        alias.asname,
                    )
                )

        elif isinstance(node, ast.Call):

            if isinstance(node.func, ast.Name):
                calls.append(
                    (
                        node.lineno,
                        node.func.id,
                    )
                )

            elif isinstance(node.func, ast.Attribute):
                calls.append(
                    (
                        node.lineno,
                        node.func.attr,
                    )
                )

    return imports, calls


# ----------------------------------------------------------------------
# Main audit
# ----------------------------------------------------------------------

def run_audit():

    lines = []

    def out(text=""):
        print(text)
        lines.append(str(text))

    py_files = python_files()

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 02: IMPLEMENTATION TRACE")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"Package root      : {PACKAGE_ROOT}")
    out(f"Python files      : {len(py_files)}")
    out()

    # ------------------------------------------------------------------
    # TEST 1 — Locate definitions
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 1 — PsiCDM function definitions")
    out("-" * 78)

    definitions_found = {}

    for path in py_files:

        definitions = get_function_definitions(path)

        for name, lineno in definitions:

            if name in TARGET_FUNCTIONS:

                definitions_found.setdefault(name, []).append(
                    (path, lineno)
                )

                out(
                    f"{name:<15} "
                    f"{rel(path)}:{lineno}"
                )

    out()

    for name in TARGET_FUNCTIONS:

        count = len(definitions_found.get(name, []))

        out(
            f"{name:<15} definitions = {count}"
        )

    out()

    # ------------------------------------------------------------------
    # TEST 2 — Search all references
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 2 — References to PsiCDM background functions")
    out("-" * 78)
    out()

    reference_map = {}

    for target in TARGET_FUNCTIONS:

        out(f"[{target}]")

        matches = []

        for path in py_files:

            occurrences = find_text_occurrences(path, target)

            for lineno, text in occurrences:
                matches.append(
                    (path, lineno, text)
                )

        reference_map[target] = matches

        if not matches:
            out("  NO REFERENCES FOUND")
        else:
            for path, lineno, text in matches:
                out(
                    f"  {rel(path)}:{lineno}"
                )
                out(
                    f"      {text}"
                )

        out()

    # ------------------------------------------------------------------
    # TEST 3 — Search published equation
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 3 — Search for published eps0/(1+z) implementation")
    out("-" * 78)
    out()

    published_matches = []

    for path in py_files:

        source = read_source(path)

        for pattern in PUBLISHED_PATTERNS:

            for match in re.finditer(pattern, source):

                lineno = line_number(
                    source,
                    match.start(),
                )

                text = source.splitlines()[
                    lineno - 1
                ].strip()

                item = (
                    path,
                    lineno,
                    text,
                )

                if item not in published_matches:
                    published_matches.append(item)

    if published_matches:

        out("Potential published-form matches found:")

        for path, lineno, text in published_matches:
            out(
                f"  {rel(path)}:{lineno}"
            )
            out(
                f"      {text}"
            )

    else:
        out(
            "NO direct implementation of "
            "eps0/(1+z) was found in package Python source."
        )

    out()

    # ------------------------------------------------------------------
    # TEST 4 — Alternative PsiCDM model code
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 4 — Files containing PsiCDMParams / eps0")
    out("-" * 78)
    out()

    model_files = []

    for path in py_files:

        source = read_source(path)

        if any(term in source for term in MODEL_TERMS):
            model_files.append(path)
            out(rel(path))

    out()
    out(
        f"Files containing PsiCDM model terms: "
        f"{len(model_files)}"
    )
    out()

    # ------------------------------------------------------------------
    # TEST 5 — Analysis pipeline imports/calls
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 5 — Analysis scripts using frozen PsiCDM background")
    out("-" * 78)
    out()

    analysis_root = PACKAGE_ROOT / "analysis"

    analysis_users = []

    if analysis_root.exists():

        for path in sorted(analysis_root.rglob("*.py")):

            if any(part in IGNORE_PARTS for part in path.parts):
                continue

            imports, calls = get_imports_and_calls(path)

            relevant_imports = []

            for item in imports:

                lineno, module, name, alias = item

                text = " ".join(
                    str(x)
                    for x in (module, name, alias)
                    if x is not None
                )

                if (
                    "psicdm" in text.lower()
                    or "PsiCDM" in text
                ):
                    relevant_imports.append(item)

            relevant_calls = [
                (lineno, name)
                for lineno, name in calls
                if name in TARGET_FUNCTIONS
            ]

            if relevant_imports or relevant_calls:

                analysis_users.append(path)

                out(rel(path))

                for lineno, module, name, alias in relevant_imports:

                    if name is None:
                        import_text = module
                    else:
                        import_text = f"{module}.{name}"

                    if alias:
                        import_text += f" as {alias}"

                    out(
                        f"    IMPORT line {lineno}: "
                        f"{import_text}"
                    )

                for lineno, name in relevant_calls:

                    out(
                        f"    CALL   line {lineno}: "
                        f"{name}(...)"
                    )

                out()

    out(
        f"Analysis scripts connected to PsiCDM: "
        f"{len(analysis_users)}"
    )
    out()

    # ------------------------------------------------------------------
    # TEST 6 — Direct source trace of main implementation
    # ------------------------------------------------------------------

    out("-" * 78)
    out("TEST 6 — Main psicdm.py source")
    out("-" * 78)
    out()

    main_model_file = (
        PACKAGE_ROOT
        / "cosmology"
        / "background"
        / "psicdm.py"
    )

    if main_model_file.exists():

        source_lines = read_source(
            main_model_file
        ).splitlines()

        for number, text in enumerate(
            source_lines,
            start=1,
        ):
            out(
                f"{number:4d}: {text}"
            )

    else:
        out(
            "ERROR: expected psicdm.py not found."
        )

    out()

    # ------------------------------------------------------------------
    # Verdict
    # ------------------------------------------------------------------

    out("=" * 78)
    out("AUDIT 02 TRACE SUMMARY")
    out("=" * 78)

    unique_main_definitions = all(
        len(definitions_found.get(name, [])) == 1
        for name in TARGET_FUNCTIONS
    )

    out(
        "Unique definitions for all target background functions: "
        f"{'YES' if unique_main_definitions else 'NO'}"
    )

    out(
        "Direct published eps0/(1+z) implementation found: "
        f"{'YES' if published_matches else 'NO'}"
    )

    out(
        "Analysis scripts connected to PsiCDM implementation: "
        f"{len(analysis_users)}"
    )

    out()

    if (
        unique_main_definitions
        and not published_matches
        and len(analysis_users) > 0
    ):
        out("AUDIT 02 VERDICT: TRACE CONFIRMED")
        out(
            "The frozen package exposes a single traced PsiCDM "
            "background implementation to the analysis layer, "
            "and no direct alternative implementation of the "
            "published eps0/(1+z) form was detected."
        )
    else:
        out("AUDIT 02 VERDICT: MANUAL REVIEW REQUIRED")
        out(
            "The source trace is not sufficiently unique for an "
            "automatic propagation conclusion."
        )

    out()
    out(
        "NOTE: This verdict concerns implementation tracing only. "
        "Scientific correctness of individual likelihoods is NOT "
        "assessed in Audit 02."
    )

    out("=" * 78)

    output_file = RESULTS_DIR / "02_trace.txt"

    output_file.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print()
    print(f"Report written to: {output_file}")


if __name__ == "__main__":
    run_audit()
