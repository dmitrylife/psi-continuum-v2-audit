#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 04b — SN calibration / H0 semantics.

Purpose
-------
Audit 04 showed that replacing the frozen PsiCDM background with the
literal published-v2 equation produces a very large change in the
Pantheon+SH0ES HF chi2.

This audit determines how much of that change is caused by:

    (1) an overall distance-modulus / absolute-calibration shift,

versus

    (2) genuine redshift-dependent shape changes.

We preserve:
    - frozen Pantheon+SH0ES HF loader,
    - frozen MU_SH0ES observations,
    - frozen STAT+SYS covariance,
    - frozen mu_from_dL conversion.

We compare two statistics:

RAW:
    chi2 = r^T C^-1 r

PROFILED:
    chi2_prof = min_Delta [(r - Delta*1)^T C^-1 (r - Delta*1)]

where Delta is one global additive distance-modulus offset.

The analytic best-fit offset is

    Delta_best = (1^T C^-1 r) / (1^T C^-1 1)

and

    chi2_prof =
        r^T C^-1 r
        - (1^T C^-1 r)^2 / (1^T C^-1 1)

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

from psi_continuum_v2.cosmology.constants import (
    C_LIGHT,
)


# ======================================================================
# Configuration
# ======================================================================

H0 = 70.0
OM0 = 0.3
N = 1.0

# Same grid as original frozen SN test / Audit 04
EPS_VALUES = np.linspace(-0.1, 0.1, 81)

SPEC_NZ = 100_000


# ======================================================================
# Published-v2 background
# ======================================================================

def H_published(z, lcdm_params, eps0):

    z = np.asarray(z, dtype=float)

    return (
        H_lcdm(z, lcdm_params)
        * (1.0 + eps0 / (1.0 + z))
    )


def cumulative_trapezoid_manual(y, x):

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

    z_values = np.asarray(
        z_values,
        dtype=float,
    )

    z_max = float(
        np.max(z_values)
    )

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

    DM_grid = (
        C_LIGHT * integral
    )

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
# SN statistics
# ======================================================================

def chi2_raw(mu_obs, mu_th, inv_cov):

    r = (
        np.asarray(mu_obs)
        - np.asarray(mu_th)
    )

    return float(
        r @ inv_cov @ r
    )


def chi2_profile_offset(
    mu_obs,
    mu_th,
    inv_cov,
):
    """
    Profile over one additive global magnitude offset.

    Residual:
        r = mu_obs - mu_th

    Model:
        r -> r - Delta

    Returns:
        chi2_profiled
        Delta_best
    """

    r = (
        np.asarray(mu_obs, dtype=float)
        - np.asarray(mu_th, dtype=float)
    )

    ones = np.ones_like(r)

    inv_r = inv_cov @ r
    inv_1 = inv_cov @ ones

    A = float(
        r @ inv_r
    )

    B = float(
        ones @ inv_r
    )

    C = float(
        ones @ inv_1
    )

    delta_best = (
        B / C
    )

    chi2_prof = (
        A - B * B / C
    )

    return (
        float(chi2_prof),
        float(delta_best),
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
# Main
# ======================================================================

def run_audit():

    lines = []

    def out(text=""):
        print(text)
        lines.append(str(text))

    # ------------------------------------------------------------------
    # Data
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

    # Invert once.
    inv_cov = np.linalg.inv(cov)

    lcdm_params = LCDMParams(
        H0=H0,
        Om0=OM0,
    )

    out("=" * 78)
    out("PSI-CONTINUUM v2 — INDEPENDENT AUDIT")
    out("AUDIT 04b: SN CALIBRATION / H0 SEMANTICS")
    out("=" * 78)
    out()

    out(f"Frozen repository : {V2_ROOT}")
    out(f"N_SN              : {sn['N']}")
    out(f"covariance shape  : {cov.shape}")
    out(f"H0 baseline       : {H0}")
    out(f"Omega_m           : {OM0}")
    out(f"eps grid          : {EPS_VALUES[0]:+.3f} .. {EPS_VALUES[-1]:+.3f}")
    out()

    # ==================================================================
    # TEST 1 — LCDM raw vs profiled
    # ==================================================================

    out("-" * 78)
    out("TEST 1 — LCDM raw and offset-profiled likelihood")
    out("-" * 78)
    out()

    dL_l = dL_lcdm(
        z,
        lcdm_params,
    )

    mu_l = mu_from_dL(
        dL_l
    )

    chi2_l_raw = chi2_raw(
        mu_obs,
        mu_l,
        inv_cov,
    )

    (
        chi2_l_prof,
        delta_l,
    ) = chi2_profile_offset(
        mu_obs,
        mu_l,
        inv_cov,
    )

    out(
        f"LCDM raw chi2       = "
        f"{chi2_l_raw:.12f}"
    )

    out(
        f"LCDM profiled chi2  = "
        f"{chi2_l_prof:.12f}"
    )

    out(
        f"LCDM best offset    = "
        f"{delta_l:+.12f} mag"
    )

    out(
        f"chi2 removed by offset = "
        f"{chi2_l_raw-chi2_l_prof:.12f}"
    )

    out()

    # ==================================================================
    # TEST 2 — Scan frozen and published models
    # ==================================================================

    out("-" * 78)
    out("TEST 2 — Raw and offset-profiled eps0 scans")
    out("-" * 78)
    out()

    n_eps = len(EPS_VALUES)

    frozen_raw = np.empty(n_eps)
    frozen_prof = np.empty(n_eps)
    frozen_offset = np.empty(n_eps)

    published_raw = np.empty(n_eps)
    published_prof = np.empty(n_eps)
    published_offset = np.empty(n_eps)

    for i, eps0 in enumerate(EPS_VALUES):

        # --------------------------------------------------------------
        # Frozen implementation
        # --------------------------------------------------------------

        psi_params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=float(eps0),
            n=N,
        )

        dL_f = dL_psicdm(
            z,
            psi_params,
        )

        mu_f = mu_from_dL(
            dL_f
        )

        frozen_raw[i] = chi2_raw(
            mu_obs,
            mu_f,
            inv_cov,
        )

        (
            frozen_prof[i],
            frozen_offset[i],
        ) = chi2_profile_offset(
            mu_obs,
            mu_f,
            inv_cov,
        )

        # --------------------------------------------------------------
        # Published v2 equation
        # --------------------------------------------------------------

        dL_p = dL_published(
            z,
            lcdm_params,
            float(eps0),
        )

        mu_p = mu_from_dL(
            dL_p
        )

        published_raw[i] = chi2_raw(
            mu_obs,
            mu_p,
            inv_cov,
        )

        (
            published_prof[i],
            published_offset[i],
        ) = chi2_profile_offset(
            mu_obs,
            mu_p,
            inv_cov,
        )

    # ==================================================================
    # TEST 3 — Minima
    # ==================================================================

    out("-" * 78)
    out("TEST 3 — Raw versus profiled minima")
    out("-" * 78)
    out()

    fr_raw_i = int(
        np.argmin(frozen_raw)
    )

    fr_prof_i = int(
        np.argmin(frozen_prof)
    )

    pub_raw_i = int(
        np.argmin(published_raw)
    )

    pub_prof_i = int(
        np.argmin(published_prof)
    )

    out("FROZEN IMPLEMENTATION")
    out(
        f"  raw best eps0       = "
        f"{EPS_VALUES[fr_raw_i]:+.6f}"
    )
    out(
        f"  raw best chi2       = "
        f"{frozen_raw[fr_raw_i]:.12f}"
    )
    out(
        f"  profiled best eps0  = "
        f"{EPS_VALUES[fr_prof_i]:+.6f}"
    )
    out(
        f"  profiled best chi2  = "
        f"{frozen_prof[fr_prof_i]:.12f}"
    )
    out(
        f"  profiled offset     = "
        f"{frozen_offset[fr_prof_i]:+.12f} mag"
    )
    out()

    out("PUBLISHED V2 EQUATION")
    out(
        f"  raw best eps0       = "
        f"{EPS_VALUES[pub_raw_i]:+.6f}"
    )
    out(
        f"  raw best chi2       = "
        f"{published_raw[pub_raw_i]:.12f}"
    )
    out(
        f"  profiled best eps0  = "
        f"{EPS_VALUES[pub_prof_i]:+.6f}"
    )
    out(
        f"  profiled best chi2  = "
        f"{published_prof[pub_prof_i]:.12f}"
    )
    out(
        f"  profiled offset     = "
        f"{published_offset[pub_prof_i]:+.12f} mag"
    )
    out()

    # ==================================================================
    # TEST 4 — Relative to profiled LCDM
    # ==================================================================

    out("-" * 78)
    out("TEST 4 — Profiled Delta chi2 relative to profiled LCDM")
    out("-" * 78)
    out()

    delta_fr_prof = (
        frozen_prof
        - chi2_l_prof
    )

    delta_pub_prof = (
        published_prof
        - chi2_l_prof
    )

    out(
        f"Frozen best profiled Delta chi2    = "
        f"{delta_fr_prof[fr_prof_i]:+.12f}"
    )

    out(
        f"Published best profiled Delta chi2 = "
        f"{delta_pub_prof[pub_prof_i]:+.12f}"
    )

    out()

    # ==================================================================
    # TEST 5 — Selected epsilon values
    # ==================================================================

    out("-" * 78)
    out("TEST 5 — Selected eps0 values")
    out("-" * 78)
    out()

    out(
        f"{'eps':>7} "
        f"{'fr_raw':>12} "
        f"{'fr_prof':>12} "
        f"{'pub_raw':>12} "
        f"{'pub_prof':>12} "
        f"{'pub_off':>11}"
    )

    out("-" * 73)

    selected = [
        -0.100,
        -0.050,
        0.000,
        0.030,
        0.050,
        0.055,
        0.075,
        0.100,
    ]

    for target in selected:

        i = nearest_index(
            EPS_VALUES,
            target,
        )

        out(
            f"{EPS_VALUES[i]:+7.3f} "
            f"{frozen_raw[i]:12.3f} "
            f"{frozen_prof[i]:12.3f} "
            f"{published_raw[i]:12.3f} "
            f"{published_prof[i]:12.3f} "
            f"{published_offset[i]:+11.6f}"
        )

    out()

    # ==================================================================
    # TEST 6 — How much of raw published improvement disappears?
    # ==================================================================

    out("-" * 78)
    out("TEST 6 — Calibration contribution to published-spec signal")
    out("-" * 78)
    out()

    raw_pub_improvement = (
        published_raw[pub_raw_i]
        - chi2_l_raw
    )

    prof_pub_improvement = (
        published_prof[pub_prof_i]
        - chi2_l_prof
    )

    out(
        f"Published raw best Delta chi2      = "
        f"{raw_pub_improvement:+.12f}"
    )

    out(
        f"Published profiled best Delta chi2 = "
        f"{prof_pub_improvement:+.12f}"
    )

    out()

    removed_signal = (
        abs(raw_pub_improvement)
        - abs(prof_pub_improvement)
    )

    if abs(raw_pub_improvement) > 0:

        fraction_removed = (
            removed_signal
            / abs(raw_pub_improvement)
        )

    else:

        fraction_removed = np.nan

    out(
        f"|Delta chi2| removed after profiling = "
        f"{removed_signal:.12f}"
    )

    out(
        f"fraction of raw |Delta chi2| removed = "
        f"{fraction_removed:.6%}"
    )

    out()

    # ==================================================================
    # TEST 7 — H(0) semantics
    # ==================================================================

    out("-" * 78)
    out("TEST 7 — H(0) semantics at raw published best fit")
    out("-" * 78)
    out()

    eps_raw_best = float(
        EPS_VALUES[pub_raw_i]
    )

    H0_published_raw_best = (
        H0
        * (1.0 + eps_raw_best)
    )

    out(
        f"baseline H_Lambda(0) = "
        f"{H0:.6f}"
    )

    out(
        f"raw published best eps0 = "
        f"{eps_raw_best:+.6f}"
    )

    out(
        f"H_Psi(0) implied by published equation = "
        f"{H0_published_raw_best:.6f}"
    )

    out()

    # ==================================================================
    # Save scan
    # ==================================================================

    scan_file = (
        RESULTS_DIR
        / "04b_sn_calibration_scan.csv"
    )

    scan = np.column_stack(
        (
            EPS_VALUES,
            frozen_raw,
            frozen_prof,
            frozen_offset,
            published_raw,
            published_prof,
            published_offset,
            delta_fr_prof,
            delta_pub_prof,
        )
    )

    np.savetxt(
        scan_file,
        scan,
        delimiter=",",
        header=(
            "eps0,"
            "frozen_raw_chi2,"
            "frozen_profiled_chi2,"
            "frozen_best_offset_mag,"
            "published_raw_chi2,"
            "published_profiled_chi2,"
            "published_best_offset_mag,"
            "frozen_profiled_delta_chi2,"
            "published_profiled_delta_chi2"
        ),
        comments="",
        fmt="%.12f",
    )

    # ==================================================================
    # Summary
    # ==================================================================

    out("=" * 78)
    out("AUDIT 04b SUMMARY")
    out("=" * 78)
    out()

    out(
        "This audit separates an overall SN magnitude/calibration "
        "offset from redshift-dependent distance-shape information."
    )

    out()

    out(
        "Raw published-spec result:"
    )

    out(
        f"  eps0 = "
        f"{EPS_VALUES[pub_raw_i]:+.6f}"
    )

    out(
        f"  Delta chi2 = "
        f"{raw_pub_improvement:+.6f}"
    )

    out()

    out(
        "Offset-profiled published-spec result:"
    )

    out(
        f"  eps0 = "
        f"{EPS_VALUES[pub_prof_i]:+.6f}"
    )

    out(
        f"  Delta chi2 = "
        f"{prof_pub_improvement:+.6f}"
    )

    out()

    out(
        "Interpretation must be based on the profiled result "
        "when asking whether the enormous raw improvement is "
        "primarily an absolute-calibration effect."
    )

    out("=" * 78)

    report_file = (
        RESULTS_DIR
        / "04b_sn_calibration.txt"
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
