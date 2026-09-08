#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 04 — Pantheon+SH0ES HF supernova likelihood.

This audit performs two logically separate tests:

A. FROZEN REPRODUCTION
   Reproduce the original frozen v2 Pantheon+ HF calculation using:
       - the original loader,
       - the original covariance,
       - the original mu_from_dL,
       - the original chi2_sn_full_cov,
       - the original 81-point eps0 grid.

B. PUBLISHED-SPEC COUNTERFACTUAL
   Keep the SAME data and SAME likelihood, changing only the background
   distance relation to that implied by the published v2 equation:

       H_Psi(z) = H_Lambda(z) * (1 + eps0/(1+z))

The frozen repository is never modified.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


# ======================================================================
# Paths
# ======================================================================

AUDIT_ROOT = Path(__file__).resolve().parent
V2_ROOT = AUDIT_ROOT.parent / "psi-continuum-v2"
RESULTS_DIR = AUDIT_ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if not V2_ROOT.exists():
    raise RuntimeError(
        f"Frozen v2 repository not found:\n{V2_ROOT}"
    )

sys.path.insert(0, str(V2_ROOT))


# ======================================================================
# Frozen v2 imports
# ======================================================================

from psi_continuum_v2.utils import get_data_path

from psi_continuum_v2.cosmology.data_loaders.pantheonplus_loader import (
    load_pantheonplus_hf,
)

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
    dL_lcdm,
    mu_from_dL,
)

from psi_continuum_v2.cosmology.background.psicdm import (
    dL_psicdm,
)

from psi_continuum_v2.cosmology.models.lcdm_params import (
    LCDMParams,
)

from psi_continuum_v2.cosmology.models.psicdm_params import (
    PsiCDMParams,
)

from psi_continuum_v2.cosmology.likelihoods.sn_likelihood import (
    chi2_sn_full_cov,
)

from psi_continuum_v2.cosmology.constants import (
    C_LIGHT,
)


# ======================================================================
# Configuration
# ======================================================================

H0 = 70.0
OM0 = 0.3
N = 1.0

# EXACT grid used by frozen sn_test_psicdm_pplus.py
EPS_VALUES = np.linspace(-0.1, 0.1, 81)

# Independent integration grid for published specification
SPEC_NZ = 100_000

# Known frozen outputs already produced by v2
REFERENCE_LCDM = 2609.19713024
REFERENCE_BEST_EPS = -0.10
REFERENCE_BEST_CHI2 = 2593.27825030

CHI2_TOL = 1e-3


# ======================================================================
# Published v2 equation
# ======================================================================

def H_published(z, lcdm_params, eps0):
    """
    Published v2 defining equation:

        H_Psi(z) = H_Lambda(z) * (1 + eps0/(1+z))
    """
    z = np.asarray(z, dtype=float)

    H_lambda = np.asarray(
        H_lcdm(z, lcdm_params),
        dtype=float,
    )

    return H_lambda * (
        1.0 + eps0 / (1.0 + z)
    )


def cumulative_trapezoid_manual(y, x):
    """
    Independent cumulative trapezoidal integration.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)

    result = np.zeros_like(y)

    result[1:] = np.cumsum(
        0.5
        * (y[:-1] + y[1:])
        * np.diff(x)
    )

    return result


def dL_published(z_values, lcdm_params, eps0):
    """
    Luminosity distance implied by the published v2 equation.

        D_M(z) = c integral_0^z dz'/H_Psi(z')
        d_L(z) = (1+z) D_M(z)
    """
    z_values = np.asarray(z_values, dtype=float)

    z_max = float(np.max(z_values))

    z_grid = np.linspace(
        0.0,
        z_max,
        SPEC_NZ,
    )

    H_grid = H_published(
        z_grid,
        lcdm_params,
        eps0,
    )

    integral = cumulative_trapezoid_manual(
        1.0 / H_grid,
        z_grid,
    )

    DM_grid = C_LIGHT * integral

    DM = np.interp(
        z_values,
        z_grid,
        DM_grid,
    )

    return (
        (1.0 + z_values)
        * DM
    )


# ======================================================================
# Helpers
# ======================================================================

def nearest_index(values, target):
    return int(
        np.argmin(
            np.abs(
                np.asarray(values)
                - target
            )
        )
    )


# ======================================================================
# Main audit
# ======================================================================

def run_audit():

    lines = []

    def out(text=""):
        print(text)
        lines.append(str(text))

    # ------------------------------------------------------------------
    # Load data EXACTLY as frozen pipeline
    # ------------------------------------------------------------------

    sn = load_pantheonplus_hf(
        get_data_path("pantheon_plus")
    )

    z = np.asarray(
        sn["z"],
        dtype=float,
    )

    mu_obs = np.asarray(
        sn["mu"],
        dtype=float,
    )

    cov = np.asarray(
        sn["cov"],
        dtype=float,
    )

    lcdm_params = LCDMParams(
        H0=H0,
        Om0=OM0,
    )

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 04: PANTHEON+SH0ES HF SUPERNOVA LIKELIHOOD")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"Data directory    : {get_data_path('pantheon_plus')}")
    out(f"N_SN              : {sn['N']}")
    out(f"covariance shape  : {cov.shape}")
    out(f"H0                : {H0}")
    out(f"Omega_m           : {OM0}")
    out(f"n                 : {N}")
    out(f"eps grid points   : {len(EPS_VALUES)}")
    out(
        f"eps range         : "
        f"{EPS_VALUES[0]:+.4f} .. "
        f"{EPS_VALUES[-1]:+.4f}"
    )
    out()

    # ==================================================================
    # TEST 1 — Exact frozen LCDM calculation
    # ==================================================================

    out("-" * 78)
    out("TEST 1 — Reproduce frozen LCDM Pantheon+ result")
    out("-" * 78)
    out()

    dL_l = dL_lcdm(
        z,
        lcdm_params,
    )

    mu_l = mu_from_dL(
        dL_l
    )

    chi2_lcdm = float(
        chi2_sn_full_cov(
            mu_obs,
            mu_l,
            cov,
        )
    )

    lcdm_difference = (
        chi2_lcdm
        - REFERENCE_LCDM
    )

    lcdm_pass = (
        abs(lcdm_difference)
        < CHI2_TOL
    )

    out(
        f"computed chi2_LCDM  = "
        f"{chi2_lcdm:.12f}"
    )

    out(
        f"reference chi2_LCDM = "
        f"{REFERENCE_LCDM:.12f}"
    )

    out(
        f"difference          = "
        f"{lcdm_difference:+.12e}"
    )

    out()

    out(
        f"TEST 1 RESULT: "
        f"{'PASS' if lcdm_pass else 'FAIL'}"
    )

    out()

    # ==================================================================
    # TEST 2 — Frozen eps0=0 null limit
    # ==================================================================

    out("-" * 78)
    out("TEST 2 — Frozen PsiCDM eps0=0 null limit")
    out("-" * 78)
    out()

    psi0 = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=0.0,
        n=N,
    )

    dL_psi0 = dL_psicdm(
        z,
        psi0,
    )

    mu_psi0 = mu_from_dL(
        dL_psi0
    )

    chi2_psi0 = float(
        chi2_sn_full_cov(
            mu_obs,
            mu_psi0,
            cov,
        )
    )

    null_difference = (
        chi2_psi0
        - chi2_lcdm
    )

    null_pass = (
        abs(null_difference)
        < CHI2_TOL
    )

    out(
        f"chi2_LCDM        = "
        f"{chi2_lcdm:.12f}"
    )

    out(
        f"chi2_Psi(eps=0) = "
        f"{chi2_psi0:.12f}"
    )

    out(
        f"difference       = "
        f"{null_difference:+.12e}"
    )

    out()

    out(
        f"TEST 2 RESULT: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )

    out()

    # ==================================================================
    # TEST 3 — Exact frozen PsiCDM scan
    # ==================================================================

    out("-" * 78)
    out("TEST 3 — Reproduce frozen PsiCDM SN scan")
    out("-" * 78)
    out()

    chi2_frozen = np.empty(
        len(EPS_VALUES),
        dtype=float,
    )

    for i, eps0 in enumerate(EPS_VALUES):

        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=float(eps0),
            n=N,
        )

        dL = dL_psicdm(
            z,
            params,
        )

        mu = mu_from_dL(
            dL
        )

        chi2_frozen[i] = (
            chi2_sn_full_cov(
                mu_obs,
                mu,
                cov,
            )
        )

    frozen_best_idx = int(
        np.argmin(
            chi2_frozen
        )
    )

    frozen_best_eps = float(
        EPS_VALUES[
            frozen_best_idx
        ]
    )

    frozen_best_chi2 = float(
        chi2_frozen[
            frozen_best_idx
        ]
    )

    frozen_delta = (
        frozen_best_chi2
        - chi2_lcdm
    )

    frozen_boundary = (
        frozen_best_idx == 0
        or
        frozen_best_idx
        == len(EPS_VALUES) - 1
    )

    frozen_scan_pass = (
        abs(
            frozen_best_eps
            - REFERENCE_BEST_EPS
        ) < 1e-12
        and
        abs(
            frozen_best_chi2
            - REFERENCE_BEST_CHI2
        ) < CHI2_TOL
    )

    out(
        f"computed best eps0 = "
        f"{frozen_best_eps:+.6f}"
    )

    out(
        f"computed best chi2 = "
        f"{frozen_best_chi2:.12f}"
    )

    out(
        f"computed Delta chi2 = "
        f"{frozen_delta:+.12f}"
    )

    out(
        f"minimum boundary   = "
        f"{'YES' if frozen_boundary else 'NO'}"
    )

    out()

    out(
        f"reference best eps0 = "
        f"{REFERENCE_BEST_EPS:+.6f}"
    )

    out(
        f"reference best chi2 = "
        f"{REFERENCE_BEST_CHI2:.12f}"
    )

    out()

    out(
        f"TEST 3 RESULT: "
        f"{'PASS' if frozen_scan_pass else 'FAIL'}"
    )

    out()

    # ==================================================================
    # Guard
    # ==================================================================

    reproduction_pass = (
        lcdm_pass
        and null_pass
        and frozen_scan_pass
    )

    if not reproduction_pass:

        out("=" * 78)
        out("AUDIT STOPPED")
        out("=" * 78)
        out()
        out(
            "Frozen SN pipeline was not reproduced."
        )
        out(
            "Published-spec comparison is therefore not interpreted."
        )

        report_file = (
            RESULTS_DIR
            / "04_sn.txt"
        )

        report_file.write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8",
        )

        print()
        print(
            f"Report written to: {report_file}"
        )

        return

    # ==================================================================
    # TEST 4 — Published v2 equation
    # ==================================================================

    out("-" * 78)
    out("TEST 4 — Published-v2 equation SN scan")
    out("-" * 78)
    out()

    chi2_published = np.empty(
        len(EPS_VALUES),
        dtype=float,
    )

    for i, eps0 in enumerate(EPS_VALUES):

        dL = dL_published(
            z,
            lcdm_params,
            float(eps0),
        )

        mu = mu_from_dL(
            dL
        )

        chi2_published[i] = (
            chi2_sn_full_cov(
                mu_obs,
                mu,
                cov,
            )
        )

    published_best_idx = int(
        np.argmin(
            chi2_published
        )
    )

    published_best_eps = float(
        EPS_VALUES[
            published_best_idx
        ]
    )

    published_best_chi2 = float(
        chi2_published[
            published_best_idx
        ]
    )

    published_delta = (
        published_best_chi2
        - chi2_lcdm
    )

    published_boundary = (
        published_best_idx == 0
        or
        published_best_idx
        == len(EPS_VALUES) - 1
    )

    out(
        f"best eps0        = "
        f"{published_best_eps:+.6f}"
    )

    out(
        f"best chi2        = "
        f"{published_best_chi2:.12f}"
    )

    out(
        f"Delta chi2       = "
        f"{published_delta:+.12f}"
    )

    out(
        f"minimum boundary = "
        f"{'YES' if published_boundary else 'NO'}"
    )

    out()

    # ==================================================================
    # TEST 5 — Selected points
    # ==================================================================

    out("-" * 78)
    out("TEST 5 — Frozen vs published-spec chi2 at selected eps0")
    out("-" * 78)
    out()

    selected_eps = [
        -0.100,
        -0.050,
        0.000,
        0.030,
        0.050,
        0.075,
        0.100,
    ]

    out(
        f"{'eps0':>8} "
        f"{'chi2_frozen':>16} "
        f"{'chi2_pub':>16} "
        f"{'pub-frozen':>16} "
        f"{'Dchi2_fr':>14} "
        f"{'Dchi2_pub':>14}"
    )

    out("-" * 92)

    for target in selected_eps:

        i = nearest_index(
            EPS_VALUES,
            target,
        )

        eps0 = float(
            EPS_VALUES[i]
        )

        cf = float(
            chi2_frozen[i]
        )

        cp = float(
            chi2_published[i]
        )

        out(
            f"{eps0:+8.3f} "
            f"{cf:16.6f} "
            f"{cp:16.6f} "
            f"{cp-cf:+16.6f} "
            f"{cf-chi2_lcdm:+14.6f} "
            f"{cp-chi2_lcdm:+14.6f}"
        )

    out()

    # ==================================================================
    # TEST 6 — Curve comparison
    # ==================================================================

    out("-" * 78)
    out("TEST 6 — Statistical impact of model mismatch")
    out("-" * 78)
    out()

    max_curve_difference = float(
        np.max(
            np.abs(
                chi2_published
                - chi2_frozen
            )
        )
    )

    max_curve_idx = int(
        np.argmax(
            np.abs(
                chi2_published
                - chi2_frozen
            )
        )
    )

    same_curve = bool(
        np.allclose(
            chi2_published,
            chi2_frozen,
            rtol=0.0,
            atol=CHI2_TOL,
        )
    )

    same_best = (
        abs(
            published_best_eps
            - frozen_best_eps
        ) < 1e-12
    )

    out(
        f"Frozen best eps0     = "
        f"{frozen_best_eps:+.6f}"
    )

    out(
        f"Published best eps0  = "
        f"{published_best_eps:+.6f}"
    )

    out(
        f"Best-eps shift       = "
        f"{published_best_eps-frozen_best_eps:+.6f}"
    )

    out()

    out(
        f"Frozen best Delta chi2    = "
        f"{frozen_delta:+.6f}"
    )

    out(
        f"Published best Delta chi2 = "
        f"{published_delta:+.6f}"
    )

    out()

    out(
        f"Maximum |chi2_pub-chi2_frozen| = "
        f"{max_curve_difference:.6f}"
    )

    out(
        f"at eps0 = "
        f"{EPS_VALUES[max_curve_idx]:+.6f}"
    )

    out()

    out(
        f"Same best eps0 = "
        f"{'YES' if same_best else 'NO'}"
    )

    out(
        f"Same chi2 curve = "
        f"{'YES' if same_curve else 'NO'}"
    )

    out()

    # ==================================================================
    # Save full scan
    # ==================================================================

    scan_file = (
        RESULTS_DIR
        / "04_sn_scan.csv"
    )

    scan_data = np.column_stack(
        (
            EPS_VALUES,
            chi2_frozen,
            chi2_published,
            chi2_frozen - chi2_lcdm,
            chi2_published - chi2_lcdm,
            chi2_published - chi2_frozen,
        )
    )

    np.savetxt(
        scan_file,
        scan_data,
        delimiter=",",
        header=(
            "eps0,"
            "chi2_frozen,"
            "chi2_published,"
            "delta_chi2_frozen,"
            "delta_chi2_published,"
            "published_minus_frozen"
        ),
        comments="",
        fmt="%.12f",
    )

    # ==================================================================
    # Final summary
    # ==================================================================

    out("=" * 78)
    out("AUDIT 04 SUMMARY")
    out("=" * 78)
    out()

    out(
        "Frozen LCDM reproduction: "
        f"{'PASS' if lcdm_pass else 'FAIL'}"
    )

    out(
        "Frozen eps0=0 null limit: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )

    out(
        "Frozen SN scan reproduction: "
        f"{'PASS' if frozen_scan_pass else 'FAIL'}"
    )

    out()

    out(
        "FROZEN SN PIPELINE: "
        f"{'PASS' if reproduction_pass else 'FAIL'}"
    )

    out()

    if same_curve:

        out(
            "PUBLISHED-SPEC IMPACT ON SN: "
            "NONE DETECTED"
        )

    else:

        out(
            "PUBLISHED-SPEC IMPACT ON SN: "
            "DETECTED"
        )

        out(
            "Changing only the background relation from the "
            "frozen implementation to the published v2 equation "
            "changes the Pantheon+ chi2(eps0) curve."
        )

    out()

    if frozen_boundary:

        out(
            "WARNING: frozen SN-only minimum is "
            "on the eps0 scan boundary."
        )

    if published_boundary:

        out(
            "WARNING: published-spec SN-only minimum is "
            "on the eps0 scan boundary."
        )

    out()

    out(
        "NOTE: This audit preserves the frozen v2 SN data, "
        "covariance matrix, distance-modulus conversion, and "
        "chi2 definition. Only the background distance relation "
        "is changed in the published-spec comparison."
    )

    out("=" * 78)

    report_file = (
        RESULTS_DIR
        / "04_sn.txt"
    )

    report_file.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )

    print()
    print(
        f"Report written to: {report_file}"
    )

    print(
        f"Scan written to:   {scan_file}"
    )


if __name__ == "__main__":
    run_audit()
