#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 03 — Distances and background observables.

Purpose
-------
Test the frozen v2 distance pipeline at two distinct levels:

A. INTERNAL CONSISTENCY
   Verify that the frozen implementation satisfies its own identities:
       D_H = c / H
       D_M = d_L / (1+z)
       d_L = (1+z) * integral[c/H(z) dz]

B. PUBLISHED-SPEC CONSISTENCY
   Independently construct distances from the published v2 equation:

       H_Psi(z) = H_Lambda(z) * (1 + eps0/(1+z))

   and compare them with the frozen implementation.

The frozen repository is never modified.
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

sys.path.insert(0, str(V2_ROOT))


# ----------------------------------------------------------------------
# Frozen imports
# ----------------------------------------------------------------------

from psi_continuum_v2.cosmology.constants import C_LIGHT

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
)

from psi_continuum_v2.cosmology.background.psicdm import (
    H_psicdm,
    dL_psicdm,
    DM_psicdm,
    DH_psicdm,
)

from psi_continuum_v2.cosmology.models.lcdm_params import (
    LCDMParams,
)

from psi_continuum_v2.cosmology.models.psicdm_params import (
    PsiCDMParams,
)


# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

H0 = 70.0
OM0 = 0.3
N = 1.0

Z_VALUES = np.array([
    0.01,
    0.05,
    0.10,
    0.30,
    0.50,
    1.00,
    1.50,
    2.00,
    3.00,
])

EPS_VALUES = [
    0.00,
    0.03,
    0.05,
    0.075,
    -0.05,
]

# Frozen dL uses nz=4000 by default.
# We use a substantially finer grid for independent integration.
SPEC_NZ = 100_000

# Tight tolerance for algebraic identities.
IDENTITY_RTOL = 1e-10

# Numerical integration comparison tolerance.
INTEGRATION_RTOL = 5e-5


# ----------------------------------------------------------------------
# Published specification
# ----------------------------------------------------------------------

def H_spec(z, lcdm_params, eps0):
    """
    Published v2 background equation:

        H_Psi(z) = H_Lambda(z) * (1 + eps0/(1+z))
    """
    z = np.asarray(z, dtype=float)

    return (
        np.asarray(H_lcdm(z, lcdm_params), dtype=float)
        * (1.0 + eps0 / (1.0 + z))
    )


def cumulative_trapezoid_manual(y, x):
    """
    Independent cumulative trapezoidal integration.

    No scipy dependency is used.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    result = np.zeros_like(y)

    dx = np.diff(x)

    result[1:] = np.cumsum(
        0.5 * (y[:-1] + y[1:]) * dx
    )

    return result


def distances_from_H(z_values, H_function, nz=SPEC_NZ):
    """
    Flat-background distances derived independently from H(z).

    Returns:
        dL, DM, DH
    """
    z_values = np.asarray(z_values, dtype=float)

    z_max = float(np.max(z_values))

    z_grid = np.linspace(
        0.0,
        z_max,
        nz,
    )

    H_grid = np.asarray(
        H_function(z_grid),
        dtype=float,
    )

    inv_H = 1.0 / H_grid

    integral = cumulative_trapezoid_manual(
        inv_H,
        z_grid,
    )

    DM_grid = C_LIGHT * integral

    DM = np.interp(
        z_values,
        z_grid,
        DM_grid,
    )

    dL = (1.0 + z_values) * DM

    DH = C_LIGHT / np.asarray(
        H_function(z_values),
        dtype=float,
    )

    return dL, DM, DH


# ----------------------------------------------------------------------
# Utility
# ----------------------------------------------------------------------

def relative_error(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    scale = np.maximum(
        np.abs(b),
        np.finfo(float).tiny,
    )

    return np.abs(a - b) / scale


def max_rel(a, b):
    return float(
        np.max(relative_error(a, b))
    )


# ----------------------------------------------------------------------
# Audit
# ----------------------------------------------------------------------

def run_audit():

    lines = []

    def out(text=""):
        print(text)
        lines.append(str(text))

    lcdm = LCDMParams(
        H0=H0,
        Om0=OM0,
    )

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 03: DISTANCES AND BACKGROUND OBSERVABLES")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"H0                : {H0}")
    out(f"Omega_m           : {OM0}")
    out(f"n                 : {N}")
    out(f"Independent grid  : {SPEC_NZ}")
    out()

    # ==================================================================
    # TEST 1
    # ==================================================================

    out("-" * 78)
    out("TEST 1 — Frozen algebraic distance identities")
    out("-" * 78)
    out()

    identity_pass = True

    for eps0 in EPS_VALUES:

        psi = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N,
        )

        H = np.asarray(
            H_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        dL = np.asarray(
            dL_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        DM = np.asarray(
            DM_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        DH = np.asarray(
            DH_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        dm_identity = max_rel(
            DM,
            dL / (1.0 + Z_VALUES),
        )

        dh_identity = max_rel(
            DH,
            C_LIGHT / H,
        )

        passed = (
            dm_identity < IDENTITY_RTOL
            and dh_identity < IDENTITY_RTOL
        )

        if not passed:
            identity_pass = False

        out(
            f"eps0={eps0:+.3f}  "
            f"DM identity={dm_identity:.12e}  "
            f"DH identity={dh_identity:.12e}  "
            f"{'PASS' if passed else 'FAIL'}"
        )

    out()
    out(
        f"TEST 1 RESULT: "
        f"{'PASS' if identity_pass else 'FAIL'}"
    )
    out()

    # ==================================================================
    # TEST 2
    # ==================================================================

    out("-" * 78)
    out("TEST 2 — Frozen dL vs independent integration of frozen H(z)")
    out("-" * 78)
    out()

    integration_pass = True

    for eps0 in EPS_VALUES:

        psi = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N,
        )

        frozen_dL = np.asarray(
            dL_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        independent_dL, _, _ = distances_from_H(
            Z_VALUES,
            lambda z: H_psicdm(z, psi),
        )

        error = max_rel(
            frozen_dL,
            independent_dL,
        )

        passed = error < INTEGRATION_RTOL

        if not passed:
            integration_pass = False

        out(
            f"eps0={eps0:+.3f}  "
            f"max relative error={error:.12e}  "
            f"{'PASS' if passed else 'FAIL'}"
        )

    out()
    out(
        f"TEST 2 RESULT: "
        f"{'PASS' if integration_pass else 'FAIL'}"
    )
    out()

    # ==================================================================
    # TEST 3
    # ==================================================================

    out("-" * 78)
    out("TEST 3 — LCDM null limit for distances")
    out("-" * 78)
    out()

    psi0 = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=0.0,
        n=N,
    )

    frozen_dL_0 = np.asarray(
        dL_psicdm(Z_VALUES, psi0),
        dtype=float,
    )

    frozen_DM_0 = np.asarray(
        DM_psicdm(Z_VALUES, psi0),
        dtype=float,
    )

    frozen_DH_0 = np.asarray(
        DH_psicdm(Z_VALUES, psi0),
        dtype=float,
    )

    lcdm_dL, lcdm_DM, lcdm_DH = distances_from_H(
        Z_VALUES,
        lambda z: H_lcdm(z, lcdm),
    )

    null_dL = max_rel(
        frozen_dL_0,
        lcdm_dL,
    )

    null_DM = max_rel(
        frozen_DM_0,
        lcdm_DM,
    )

    null_DH = max_rel(
        frozen_DH_0,
        lcdm_DH,
    )

    null_pass = (
        null_dL < INTEGRATION_RTOL
        and null_DM < INTEGRATION_RTOL
        and null_DH < INTEGRATION_RTOL
    )

    out(
        f"dL max relative error = {null_dL:.12e}"
    )
    out(
        f"DM max relative error = {null_DM:.12e}"
    )
    out(
        f"DH max relative error = {null_DH:.12e}"
    )
    out()

    out(
        f"TEST 3 RESULT: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )
    out()

    # ==================================================================
    # TEST 4
    # ==================================================================

    out("-" * 78)
    out("TEST 4 — Frozen distances vs published v2 specification")
    out("-" * 78)
    out()

    spec_match_pass = True

    for eps0 in EPS_VALUES:

        psi = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N,
        )

        frozen_dL = np.asarray(
            dL_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        frozen_DM = np.asarray(
            DM_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        frozen_DH = np.asarray(
            DH_psicdm(Z_VALUES, psi),
            dtype=float,
        )

        spec_dL, spec_DM, spec_DH = distances_from_H(
            Z_VALUES,
            lambda z: H_spec(
                z,
                lcdm,
                eps0,
            ),
        )

        err_dL = max_rel(
            frozen_dL,
            spec_dL,
        )

        err_DM = max_rel(
            frozen_DM,
            spec_DM,
        )

        err_DH = max_rel(
            frozen_DH,
            spec_DH,
        )

        passed = (
            err_dL < INTEGRATION_RTOL
            and err_DM < INTEGRATION_RTOL
            and err_DH < INTEGRATION_RTOL
        )

        if not passed:
            spec_match_pass = False

        out(
            f"eps0={eps0:+.3f}  "
            f"dL={err_dL:.8e}  "
            f"DM={err_DM:.8e}  "
            f"DH={err_DH:.8e}  "
            f"{'PASS' if passed else 'FAIL'}"
        )

    out()
    out(
        f"TEST 4 RESULT: "
        f"{'PASS' if spec_match_pass else 'FAIL'}"
    )
    out()

    # ==================================================================
    # TEST 5
    # ==================================================================

    out("-" * 78)
    out("TEST 5 — Point-by-point comparison for eps0 = +0.03")
    out("-" * 78)
    out()

    eps0 = 0.03

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=eps0,
        n=N,
    )

    frozen_dL = np.asarray(
        dL_psicdm(Z_VALUES, psi),
        dtype=float,
    )

    frozen_DM = np.asarray(
        DM_psicdm(Z_VALUES, psi),
        dtype=float,
    )

    frozen_DH = np.asarray(
        DH_psicdm(Z_VALUES, psi),
        dtype=float,
    )

    spec_dL, spec_DM, spec_DH = distances_from_H(
        Z_VALUES,
        lambda z: H_spec(
            z,
            lcdm,
            eps0,
        ),
    )

    out(
        f"{'z':>6} "
        f"{'dL_code':>13} "
        f"{'dL_spec':>13} "
        f"{'dL_rel':>11} "
        f"{'DH_code':>13} "
        f"{'DH_spec':>13} "
        f"{'DH_rel':>11}"
    )

    out("-" * 88)

    for i, z in enumerate(Z_VALUES):

        dL_err = relative_error(
            frozen_dL[i],
            spec_dL[i],
        )

        DH_err = relative_error(
            frozen_DH[i],
            spec_DH[i],
        )

        out(
            f"{z:6.2f} "
            f"{frozen_dL[i]:13.5f} "
            f"{spec_dL[i]:13.5f} "
            f"{float(dL_err):11.5e} "
            f"{frozen_DH[i]:13.5f} "
            f"{spec_DH[i]:13.5f} "
            f"{float(DH_err):11.5e}"
        )

    out()

    # ==================================================================
    # Summary
    # ==================================================================

    out("=" * 78)
    out("AUDIT 03 SUMMARY")
    out("=" * 78)

    out(
        "Frozen algebraic identities: "
        f"{'PASS' if identity_pass else 'FAIL'}"
    )

    out(
        "Frozen numerical integration: "
        f"{'PASS' if integration_pass else 'FAIL'}"
    )

    out(
        "LCDM null limit: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )

    out(
        "Published-v2 distance consistency: "
        f"{'PASS' if spec_match_pass else 'FAIL'}"
    )

    out()

    internal_pass = (
        identity_pass
        and integration_pass
        and null_pass
    )

    if internal_pass:
        out(
            "INTERNAL FROZEN PIPELINE: PASS"
        )
        out(
            "The frozen distance machinery is internally "
            "consistent with the H(z) implementation it receives."
        )
    else:
        out(
            "INTERNAL FROZEN PIPELINE: FAIL"
        )

    out()

    if spec_match_pass:
        out(
            "PUBLISHED SPECIFICATION: PASS"
        )
    else:
        out(
            "PUBLISHED SPECIFICATION: FAIL"
        )
        out(
            "Distances generated by the frozen pipeline differ "
            "from distances implied by the published v2 H(z) equation."
        )

    out()
    out(
        "Interpretation: internal numerical consistency and "
        "consistency with the published scientific model are "
        "separate audit questions."
    )

    out("=" * 78)

    output_file = RESULTS_DIR / "03_distances.txt"

    output_file.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print()
    print(
        f"Report written to: {output_file}"
    )


if __name__ == "__main__":
    run_audit()
