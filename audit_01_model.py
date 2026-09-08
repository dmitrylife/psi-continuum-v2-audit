#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 01 — Published model equation vs frozen implementation.

This script DOES NOT modify psi-continuum-v2.

Published v2 specification:

    H_Psi(z) = H_Lambda(z) * (1 + eps0 / (1 + z))

The script imports the frozen v2 implementation and compares its numerical
output against the published equation.

Status:
    PASS  — implementation agrees with the published equation
    FAIL  — implementation disagrees with the published equation
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------

AUDIT_ROOT = Path(__file__).resolve().parent
V2_ROOT = AUDIT_ROOT.parent / "psi-continuum-v2"
RESULTS_DIR = AUDIT_ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if not V2_ROOT.exists():
    raise RuntimeError(
        f"Frozen v2 repository not found:\n{V2_ROOT}"
    )

# Import frozen package without installing/changing it.
sys.path.insert(0, str(V2_ROOT))


# ----------------------------------------------------------------------
# Frozen v2 imports
# ----------------------------------------------------------------------

from psi_continuum_v2.cosmology.background.lcdm import H_lcdm
from psi_continuum_v2.cosmology.background.psicdm import H_psicdm
from psi_continuum_v2.cosmology.models.lcdm_params import LCDMParams
from psi_continuum_v2.cosmology.models.psicdm_params import PsiCDMParams


# ----------------------------------------------------------------------
# Audit configuration
# ----------------------------------------------------------------------

H0 = 70.0
OM0 = 0.3
N = 1.0

Z_VALUES = np.array([
    0.0,
    0.1,
    0.5,
    1.0,
    2.0,
    5.0,
    10.0,
])

EPS_VALUES = [
    -0.10,
    -0.05,
    0.00,
    0.03,
    0.05,
    0.075,
    0.10,
]

# Numerical tolerance for equality with the published specification.
RTOL = 1e-10
ATOL = 1e-12


# ----------------------------------------------------------------------
# Published v2 specification
# ----------------------------------------------------------------------

def H_psicdm_spec(z, lcdm_params, eps0):
    """
    Published Psi-Continuum v2 equation:

        H_Psi(z) = H_Lambda(z) * (1 + eps0 / (1 + z))
    """
    z = np.asarray(z, dtype=float)

    h_lambda = np.asarray(
        H_lcdm(z, lcdm_params),
        dtype=float,
    )

    return h_lambda * (1.0 + eps0 / (1.0 + z))


# ----------------------------------------------------------------------
# Main audit
# ----------------------------------------------------------------------

def run_audit():
    lines = []

    def out(text=""):
        print(text)
        lines.append(str(text))

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 01: PUBLISHED MODEL EQUATION")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"H0                : {H0}")
    out(f"Omega_m           : {OM0}")
    out(f"n                 : {N}")
    out()

    out("Published specification:")
    out("H_Psi(z) = H_Lambda(z) * (1 + eps0 / (1 + z))")
    out()

    lcdm_params = LCDMParams(
        H0=H0,
        Om0=OM0,
    )

    overall_pass = True

    # --------------------------------------------------------------
    # Test 1 — exact LCDM recovery at eps0 = 0
    # --------------------------------------------------------------

    out("-" * 78)
    out("TEST 1 — LCDM recovery at eps0 = 0")
    out("-" * 78)

    psi_zero = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=0.0,
        n=N,
    )

    h_lcdm = np.asarray(
        H_lcdm(Z_VALUES, lcdm_params),
        dtype=float,
    )

    h_psi_zero = np.asarray(
        H_psicdm(Z_VALUES, psi_zero),
        dtype=float,
    )

    lcdm_recovery = np.allclose(
        h_psi_zero,
        h_lcdm,
        rtol=RTOL,
        atol=ATOL,
    )

    out(f"Result: {'PASS' if lcdm_recovery else 'FAIL'}")
    out(f"max |H_Psi/H_LCDM - 1| = "
        f"{np.max(np.abs(h_psi_zero / h_lcdm - 1.0)):.12e}")
    out()

    if not lcdm_recovery:
        overall_pass = False

    # --------------------------------------------------------------
    # Test 2 — published equation
    # --------------------------------------------------------------

    out("-" * 78)
    out("TEST 2 — Frozen implementation vs published equation")
    out("-" * 78)
    out()

    header = (
        f"{'eps0':>8} "
        f"{'z':>8} "
        f"{'H_code':>15} "
        f"{'H_spec':>15} "
        f"{'rel.error':>15}"
    )

    out(header)
    out("-" * len(header))

    max_rel_error = 0.0
    worst_case = None

    for eps0 in EPS_VALUES:

        psi_params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N,
        )

        h_code = np.asarray(
            H_psicdm(Z_VALUES, psi_params),
            dtype=float,
        )

        h_spec = np.asarray(
            H_psicdm_spec(
                Z_VALUES,
                lcdm_params,
                eps0,
            ),
            dtype=float,
        )

        for z, hc, hs in zip(Z_VALUES, h_code, h_spec):

            rel_error = abs(hc - hs) / abs(hs)

            if rel_error > max_rel_error:
                max_rel_error = rel_error
                worst_case = (eps0, z, hc, hs)

            out(
                f"{eps0:8.3f} "
                f"{z:8.3f} "
                f"{hc:15.8f} "
                f"{hs:15.8f} "
                f"{rel_error:15.8e}"
            )

    equation_pass = max_rel_error <= (RTOL + ATOL)

    out()
    out(f"Result: {'PASS' if equation_pass else 'FAIL'}")
    out(f"Maximum relative error = {max_rel_error:.12e}")

    if worst_case is not None:
        eps0, z, hc, hs = worst_case

        out(
            "Worst case: "
            f"eps0={eps0}, z={z}, "
            f"H_code={hc:.12f}, "
            f"H_spec={hs:.12f}"
        )

    if not equation_pass:
        overall_pass = False

    out()

    # --------------------------------------------------------------
    # Test 3 — present-day normalization
    # --------------------------------------------------------------

    out("-" * 78)
    out("TEST 3 — z = 0 behaviour")
    out("-" * 78)
    out()

    out(
        f"{'eps0':>8} "
        f"{'H_code(0)':>15} "
        f"{'H_spec(0)':>15} "
        f"{'spec/H0':>12}"
    )

    out("-" * 56)

    for eps0 in EPS_VALUES:

        psi_params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N,
        )

        h_code_0 = float(
            np.asarray(H_psicdm(0.0, psi_params))
        )

        h_spec_0 = H0 * (1.0 + eps0)

        out(
            f"{eps0:8.3f} "
            f"{h_code_0:15.8f} "
            f"{h_spec_0:15.8f} "
            f"{h_spec_0 / H0:12.8f}"
        )

    out()
    out("Published v2 requires:")
    out("H_Psi(0) = H0 * (1 + eps0)")
    out()

    # --------------------------------------------------------------
    # Final verdict
    # --------------------------------------------------------------

    out("=" * 78)

    if overall_pass:
        out("AUDIT 01 VERDICT: PASS")
        out(
            "Frozen implementation is numerically consistent "
            "with the published v2 model equation."
        )
    else:
        out("AUDIT 01 VERDICT: FAIL")
        out(
            "Frozen implementation is NOT numerically consistent "
            "with the published v2 model equation."
        )

    out("=" * 78)

    # --------------------------------------------------------------
    # Save report
    # --------------------------------------------------------------

    output_file = RESULTS_DIR / "01_model.txt"

    output_file.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print()
    print(f"Report written to: {output_file}")

    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(run_audit())
