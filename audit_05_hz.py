#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 05 — H(z) / Cosmic Chronometers.

Purpose
-------
1. Reproduce the frozen v2 H(z) likelihood and eps0 scan exactly.
2. Verify the eps0=0 null limit.
3. Replace only the frozen PsiCDM H(z) relation by the published
   v2 defining equation:

       H_Psi(z) = H_Lambda(z) * (1 + eps0/(1+z))

4. Compare the resulting chi2(eps0) curves.

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

from psi_continuum_v2.cosmology.data_loaders.hz_loader import (
    load_hz_compilation,
)

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
)

from psi_continuum_v2.cosmology.background.psicdm import (
    H_psicdm,
)

from psi_continuum_v2.cosmology.models.lcdm_params import (
    LCDMParams,
)

from psi_continuum_v2.cosmology.models.psicdm_params import (
    PsiCDMParams,
)

from psi_continuum_v2.cosmology.likelihoods.hz_likelihood import (
    chi2_hz,
)


# ======================================================================
# Configuration
# ======================================================================

H0 = 70.0
OM0 = 0.3
N_POWER = 1.0

# EXACT grid used by frozen analysis/tests/hz_test_psicdm.py
EPS_VALUES = np.linspace(-0.2, 0.2, 81)

# Previously generated frozen reference values.
# These are used only as reproduction checks.
REFERENCE_N = 32
REFERENCE_LCDM = 11.78688559
REFERENCE_BEST_EPS = 0.055
REFERENCE_BEST_CHI2 = 11.74585998

CHI2_TOL = 1e-6


# ======================================================================
# Published v2 model
# ======================================================================

def H_published(z, params):
    """
    Published v2 defining equation:

        H_Psi(z) =
            H_Lambda(z) * (1 + eps0/(1+z))

    params is PsiCDMParams so that this function has the same call
    signature expected by frozen chi2_hz().
    """

    z = np.asarray(z, dtype=float)

    lcdm_params = LCDMParams(
        H0=params.H0,
        Om0=params.Om0,
    )

    return (
        H_lcdm(z, lcdm_params)
        * (1.0 + params.eps0 / (1.0 + z))
    )


# ======================================================================
# Independent diagonal H(z) chi2
# ======================================================================

def chi2_hz_independent(hzdata, H_th):
    """
    Independent implementation of the likelihood documented by frozen v2:

        chi2 = sum[(H_obs - H_th)^2 / sigma_H^2]
    """

    H_obs = np.asarray(
        hzdata["Hz"],
        dtype=float,
    )

    sigma = np.asarray(
        hzdata["sigma_Hz"],
        dtype=float,
    )

    H_th = np.asarray(
        H_th,
        dtype=float,
    )

    return float(
        np.sum(
            (H_obs - H_th) ** 2
            / sigma**2
        )
    )


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
    # Load frozen dataset
    # ------------------------------------------------------------------

    hzdata = load_hz_compilation(
        get_data_path("hz")
    )

    z = np.asarray(
        hzdata["z"],
        dtype=float,
    )

    H_obs = np.asarray(
        hzdata["Hz"],
        dtype=float,
    )

    sigma_H = np.asarray(
        hzdata["sigma_Hz"],
        dtype=float,
    )

    lcdm_params = LCDMParams(
        H0=H0,
        Om0=OM0,
    )

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 05: H(z) / COSMIC CHRONOMETERS")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"Data directory    : {get_data_path('hz')}")
    out(f"N_Hz              : {hzdata['N']}")
    out(f"H0                : {H0}")
    out(f"Omega_m           : {OM0}")
    out(f"n                 : {N_POWER}")
    out(f"eps grid points   : {len(EPS_VALUES)}")
    out(
        f"eps range         : "
        f"{EPS_VALUES[0]:+.4f} .. "
        f"{EPS_VALUES[-1]:+.4f}"
    )
    out()

    out(
        f"z range           : "
        f"{np.min(z):.6f} .. {np.max(z):.6f}"
    )

    out(
        f"H_obs range       : "
        f"{np.min(H_obs):.6f} .. {np.max(H_obs):.6f}"
    )

    out()

    # ==================================================================
    # TEST 1 — Dataset and LCDM reproduction
    # ==================================================================

    out("-" * 78)
    out("TEST 1 — Reproduce frozen LCDM H(z) result")
    out("-" * 78)
    out()

    chi2_lcdm = float(
        chi2_hz(
            hzdata,
            H_lcdm,
            lcdm_params,
        )
    )

    H_l = H_lcdm(
        z,
        lcdm_params,
    )

    chi2_lcdm_ind = (
        chi2_hz_independent(
            hzdata,
            H_l,
        )
    )

    n_pass = (
        hzdata["N"]
        == REFERENCE_N
    )

    likelihood_pass = (
        abs(
            chi2_lcdm
            - chi2_lcdm_ind
        )
        < 1e-12
    )

    reference_pass = (
        abs(
            chi2_lcdm
            - REFERENCE_LCDM
        )
        < CHI2_TOL
    )

    test1_pass = (
        n_pass
        and likelihood_pass
        and reference_pass
    )

    out(
        f"N_Hz                 = "
        f"{hzdata['N']}"
    )

    out(
        f"frozen chi2_LCDM     = "
        f"{chi2_lcdm:.12f}"
    )

    out(
        f"independent chi2     = "
        f"{chi2_lcdm_ind:.12f}"
    )

    out(
        f"reference chi2_LCDM  = "
        f"{REFERENCE_LCDM:.12f}"
    )

    out(
        f"frozen-independent   = "
        f"{chi2_lcdm-chi2_lcdm_ind:+.12e}"
    )

    out(
        f"frozen-reference     = "
        f"{chi2_lcdm-REFERENCE_LCDM:+.12e}"
    )

    out()

    out(
        f"TEST 1 RESULT: "
        f"{'PASS' if test1_pass else 'FAIL'}"
    )

    out()

    # ==================================================================
    # TEST 2 — eps0=0 null limit
    # ==================================================================

    out("-" * 78)
    out("TEST 2 — Frozen and published eps0=0 null limits")
    out("-" * 78)
    out()

    psi0 = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=0.0,
        n=N_POWER,
    )

    H_frozen_0 = H_psicdm(
        z,
        psi0,
    )

    H_published_0 = H_published(
        z,
        psi0,
    )

    max_frozen_null = float(
        np.max(
            np.abs(
                H_frozen_0
                - H_l
            )
        )
    )

    max_published_null = float(
        np.max(
            np.abs(
                H_published_0
                - H_l
            )
        )
    )

    chi2_frozen_0 = float(
        chi2_hz(
            hzdata,
            H_psicdm,
            psi0,
        )
    )

    chi2_published_0 = float(
        chi2_hz(
            hzdata,
            H_published,
            psi0,
        )
    )

    null_pass = (
        max_frozen_null < 1e-12
        and
        max_published_null < 1e-12
        and
        abs(
            chi2_frozen_0
            - chi2_lcdm
        ) < 1e-12
        and
        abs(
            chi2_published_0
            - chi2_lcdm
        ) < 1e-12
    )

    out(
        f"max |H_frozen(eps=0)-H_LCDM|    = "
        f"{max_frozen_null:.12e}"
    )

    out(
        f"max |H_published(eps=0)-H_LCDM| = "
        f"{max_published_null:.12e}"
    )

    out(
        f"chi2_LCDM            = "
        f"{chi2_lcdm:.12f}"
    )

    out(
        f"chi2_frozen(eps=0)   = "
        f"{chi2_frozen_0:.12f}"
    )

    out(
        f"chi2_published(eps=0)= "
        f"{chi2_published_0:.12f}"
    )

    out()

    out(
        f"TEST 2 RESULT: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )

    out()

    # ==================================================================
    # TEST 3 — Frozen scan reproduction
    # ==================================================================

    out("-" * 78)
    out("TEST 3 — Reproduce frozen H(z) eps0 scan")
    out("-" * 78)
    out()

    frozen_chi2 = np.empty(
        len(EPS_VALUES),
        dtype=float,
    )

    for i, eps0 in enumerate(EPS_VALUES):

        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=float(eps0),
            n=N_POWER,
        )

        frozen_chi2[i] = (
            chi2_hz(
                hzdata,
                H_psicdm,
                params,
            )
        )

    frozen_best_idx = int(
        np.argmin(
            frozen_chi2
        )
    )

    frozen_best_eps = float(
        EPS_VALUES[
            frozen_best_idx
        ]
    )

    frozen_best_chi2 = float(
        frozen_chi2[
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
    # Reproduction guard
    # ==================================================================

    reproduction_pass = (
        test1_pass
        and null_pass
        and frozen_scan_pass
    )

    if not reproduction_pass:

        out("=" * 78)
        out("AUDIT STOPPED")
        out("=" * 78)
        out()
        out(
            "Frozen H(z) pipeline was not reproduced."
        )
        out(
            "Published-spec comparison is therefore not interpreted."
        )

        report_file = (
            RESULTS_DIR
            / "05_hz.txt"
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
    # TEST 4 — Published equation scan
    # ==================================================================

    out("-" * 78)
    out("TEST 4 — Published-v2 equation H(z) scan")
    out("-" * 78)
    out()

    published_chi2 = np.empty(
        len(EPS_VALUES),
        dtype=float,
    )

    for i, eps0 in enumerate(EPS_VALUES):

        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=float(eps0),
            n=N_POWER,
        )

        published_chi2[i] = (
            chi2_hz(
                hzdata,
                H_published,
                params,
            )
        )

    published_best_idx = int(
        np.argmin(
            published_chi2
        )
    )

    published_best_eps = float(
        EPS_VALUES[
            published_best_idx
        ]
    )

    published_best_chi2 = float(
        published_chi2[
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
    # TEST 5 — Selected eps values
    # ==================================================================

    out("-" * 78)
    out("TEST 5 — Frozen vs published H(z) likelihood")
    out("-" * 78)
    out()

    out(
        f"{'eps0':>8} "
        f"{'chi2_frozen':>16} "
        f"{'chi2_pub':>16} "
        f"{'pub-frozen':>16} "
        f"{'Dchi2_fr':>14} "
        f"{'Dchi2_pub':>14}"
    )

    out("-" * 92)

    selected_eps = [
        -0.200,
        -0.100,
        -0.050,
        0.000,
        0.030,
        0.050,
        0.055,
        0.075,
        0.100,
        0.200,
    ]

    for target in selected_eps:

        i = nearest_index(
            EPS_VALUES,
            target,
        )

        eps0 = float(
            EPS_VALUES[i]
        )

        cf = float(
            frozen_chi2[i]
        )

        cp = float(
            published_chi2[i]
        )

        out(
            f"{eps0:+8.3f} "
            f"{cf:16.8f} "
            f"{cp:16.8f} "
            f"{cp-cf:+16.8f} "
            f"{cf-chi2_lcdm:+14.8f} "
            f"{cp-chi2_lcdm:+14.8f}"
        )

    out()

    # ==================================================================
    # TEST 6 — Pointwise H(z) differences
    # ==================================================================

    out("-" * 78)
    out("TEST 6 — Pointwise H(z) model differences")
    out("-" * 78)
    out()

    diagnostic_eps = [
        0.030,
        0.055,
        0.075,
        -0.050,
    ]

    for target in diagnostic_eps:

        i = nearest_index(
            EPS_VALUES,
            target,
        )

        eps0 = float(
            EPS_VALUES[i]
        )

        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N_POWER,
        )

        H_f = H_psicdm(
            z,
            params,
        )

        H_p = H_published(
            z,
            params,
        )

        rel = (
            (H_p - H_f)
            / H_f
        )

        out(
            f"eps0={eps0:+.3f}: "
            f"max |H_pub-H_frozen| = "
            f"{np.max(np.abs(H_p-H_f)):.8f}, "
            f"max |relative difference| = "
            f"{np.max(np.abs(rel)):.8e}"
        )

    out()

    # ==================================================================
    # TEST 7 — H(0) semantics
    # ==================================================================

    out("-" * 78)
    out("TEST 7 — H(0) normalization semantics")
    out("-" * 78)
    out()

    selected_h0_eps = [
        -0.10,
        0.00,
        0.055,
        0.10,
    ]

    out(
        f"{'eps0':>8} "
        f"{'H0_frozen':>14} "
        f"{'H0_published':>16}"
    )

    out("-" * 42)

    for eps0 in selected_h0_eps:

        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=eps0,
            n=N_POWER,
        )

        H0_f = float(
            np.asarray(
                H_psicdm(
                    np.array([0.0]),
                    params,
                )
            )[0]
        )

        H0_p = float(
            np.asarray(
                H_published(
                    np.array([0.0]),
                    params,
                )
            )[0]
        )

        out(
            f"{eps0:+8.3f} "
            f"{H0_f:14.6f} "
            f"{H0_p:16.6f}"
        )

    out()

    # ==================================================================
    # TEST 8 — Statistical impact
    # ==================================================================

    out("-" * 78)
    out("TEST 8 — Statistical impact of implementation mismatch")
    out("-" * 78)
    out()

    curve_difference = (
        published_chi2
        - frozen_chi2
    )

    max_curve_idx = int(
        np.argmax(
            np.abs(
                curve_difference
            )
        )
    )

    max_curve_difference = float(
        np.abs(
            curve_difference[
                max_curve_idx
            ]
        )
    )

    same_curve = bool(
        np.allclose(
            frozen_chi2,
            published_chi2,
            rtol=0.0,
            atol=CHI2_TOL,
        )
    )

    same_best = bool(
        abs(
            frozen_best_eps
            - published_best_eps
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
        f"{frozen_delta:+.12f}"
    )

    out(
        f"Published best Delta chi2 = "
        f"{published_delta:+.12f}"
    )

    out()

    out(
        f"Maximum |chi2_pub-chi2_frozen| = "
        f"{max_curve_difference:.12f}"
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
    # Save scan
    # ==================================================================

    scan_file = (
        RESULTS_DIR
        / "05_hz_scan.csv"
    )

    scan = np.column_stack(
        (
            EPS_VALUES,
            frozen_chi2,
            published_chi2,
            frozen_chi2 - chi2_lcdm,
            published_chi2 - chi2_lcdm,
            published_chi2 - frozen_chi2,
        )
    )

    np.savetxt(
        scan_file,
        scan,
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
    # Summary
    # ==================================================================

    out("=" * 78)
    out("AUDIT 05 SUMMARY")
    out("=" * 78)
    out()

    out(
        "Frozen H(z) dataset/likelihood reproduction: "
        f"{'PASS' if test1_pass else 'FAIL'}"
    )

    out(
        "eps0=0 null limit: "
        f"{'PASS' if null_pass else 'FAIL'}"
    )

    out(
        "Frozen H(z) scan reproduction: "
        f"{'PASS' if frozen_scan_pass else 'FAIL'}"
    )

    out()

    out(
        "FROZEN H(z) PIPELINE: "
        f"{'PASS' if reproduction_pass else 'FAIL'}"
    )

    out()

    if same_curve:

        out(
            "PUBLISHED-SPEC IMPACT ON H(z): "
            "NONE DETECTED"
        )

    else:

        out(
            "PUBLISHED-SPEC IMPACT ON H(z): "
            "DETECTED"
        )

        out(
            "Replacing the frozen PsiCDM H(z) relation by the "
            "published v2 defining equation changes the "
            "H(z) chi2(eps0) curve."
        )

    out()

    if frozen_boundary:

        out(
            "WARNING: frozen H(z) minimum lies on the scan boundary."
        )

    if published_boundary:

        out(
            "WARNING: published-spec H(z) minimum lies on the scan boundary."
        )

    out()

    out(
        "NOTE: The H(z) dataset, observational uncertainties, "
        "diagonal chi2 definition, H0, Omega_m, n, and eps0 scan "
        "are identical in the frozen and published-spec calculations. "
        "Only the PsiCDM H(z) relation is changed."
    )

    out("=" * 78)

    report_file = (
        RESULTS_DIR
        / "05_hz.txt"
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
