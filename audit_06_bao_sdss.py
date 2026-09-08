#!/usr/bin/env python3
"""
Audit 06B — SDSS/BOSS DR12 BAO statistical impact

Compares three interpretations of exactly the same DR12 data vector
and covariance:

1. FROZEN RAW
   Exact frozen-v2 implementation:
       model = [DM(z), H(z)]

2. FROZEN + CORRECT BOSS CONVENTION
       model = [
           DM(z) * rd_fid / rd,
           H(z)  * rd / rd_fid
       ]

3. PUBLISHED-v2 EQUATION + CORRECT BOSS CONVENTION
       H_psi(z) = H_LCDM(z) * (1 + eps0/(1+z))
   with DM obtained by independent numerical integration.

The frozen v2 repository is never modified.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


# ============================================================
# Paths
# ============================================================

AUDIT_DIR = Path(__file__).resolve().parent
FROZEN_ROOT = (AUDIT_DIR / "../psi-continuum-v2").resolve()

if not FROZEN_ROOT.exists():
    raise RuntimeError(f"Frozen v2 repository not found: {FROZEN_ROOT}")

sys.path.insert(0, str(FROZEN_ROOT))

DATA_DIR = FROZEN_ROOT / "data" / "bao"


# ============================================================
# Frozen-v2 imports
# ============================================================

from psi_continuum_v2.cosmology.data_loaders.bao_loader import load_bao_dr12
from psi_continuum_v2.cosmology.likelihoods.bao_likelihood import (
    bao_vector_model,
    chi2_bao,
)

from psi_continuum_v2.cosmology.models.lcdm_params import LCDMParams
from psi_continuum_v2.cosmology.models.psicdm_params import PsiCDMParams

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
    dL_lcdm,
)

from psi_continuum_v2.cosmology.background.psicdm import (
    H_psicdm,
    dL_psicdm,
)


# ============================================================
# Constants / fixed v2 assumptions
# ============================================================

C_LIGHT = 299792.458  # km/s

H0 = 70.0
OM0 = 0.3
N = 1.0

# Frozen v2 value
RD = 147.0  # Mpc

# BOSS DR12 fiducial sound horizon associated with consensus vector
RD_FID = 147.78  # Mpc

# Same epsilon grid as frozen joint scan
EPS_GRID = np.linspace(-0.10, 0.10, 201)


# ============================================================
# Published-v2 background equation
# ============================================================

def H_published(z, eps0):
    """
    Literal published-v2 defining equation:

        H_psi(z) = H_LCDM(z) * (1 + eps0/(1+z))
    """
    z = np.asarray(z, dtype=float)

    lcdm = LCDMParams(
        H0=H0,
        Om0=OM0,
        rd=RD,
    )

    return H_lcdm(z, lcdm) * (1.0 + eps0 / (1.0 + z))


def DM_published(z_eval, eps0):
    """
    Independent numerical integration:

        DM(z) = c * integral_0^z dz' / H_published(z')

    Flat background assumed, as in frozen v2.
    """
    z_eval = np.asarray(z_eval, dtype=float)

    zmax = float(np.max(z_eval))

    # Fine grid starting exactly at z=0.
    # Dense enough that integration error is negligible for this audit.
    grid = np.linspace(0.0, zmax, 200_001)

    H_grid = H_published(grid, eps0)

    chi = cumulative_trapezoid(
        C_LIGHT / H_grid,
        grid,
        initial=0.0,
    )

    return np.interp(z_eval, grid, chi)


# ============================================================
# Model vector builders
# ============================================================

def vector_frozen_raw(z, eps0):
    """
    Exact frozen-v2 SDSS implementation:

        [DM, H, DM, H, DM, H]

    No rd/rdfid scaling.
    """
    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=eps0,
        n=N,
        rd=RD,
    )

    dL = dL_psicdm(z, psi)
    DM = dL / (1.0 + z)
    H = H_psicdm(z, psi)

    return bao_vector_model(z, DM, H)


def vector_frozen_corrected(z, eps0):
    """
    Frozen background, but using correct BOSS DR12 observable convention:

        DM * rd_fid / rd
        H  * rd / rd_fid
    """
    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=eps0,
        n=N,
        rd=RD,
    )

    dL = dL_psicdm(z, psi)
    DM = dL / (1.0 + z)
    H = H_psicdm(z, psi)

    DM_obs = DM * (RD_FID / RD)
    H_obs = H * (RD / RD_FID)

    return bao_vector_model(z, DM_obs, H_obs)


def vector_published_corrected(z, eps0):
    """
    Literal published-v2 background equation + correct BOSS convention.
    """
    DM = DM_published(z, eps0)
    H = H_published(z, eps0)

    DM_obs = DM * (RD_FID / RD)
    H_obs = H * (RD / RD_FID)

    return bao_vector_model(z, DM_obs, H_obs)


# ============================================================
# LCDM cross-checks
# ============================================================

def vector_lcdm_frozen_raw(z):
    lcdm = LCDMParams(
        H0=H0,
        Om0=OM0,
        rd=RD,
    )

    dL = dL_lcdm(z, lcdm)
    DM = dL / (1.0 + z)
    H = H_lcdm(z, lcdm)

    return bao_vector_model(z, DM, H)


def vector_lcdm_corrected(z):
    lcdm = LCDMParams(
        H0=H0,
        Om0=OM0,
        rd=RD,
    )

    dL = dL_lcdm(z, lcdm)
    DM = dL / (1.0 + z)
    H = H_lcdm(z, lcdm)

    return bao_vector_model(
        z,
        DM * (RD_FID / RD),
        H * (RD / RD_FID),
    )


# ============================================================
# Scan helper
# ============================================================

def scan_eps(data_vec, cov, z, builder):
    vals = []

    for eps in EPS_GRID:
        vec = builder(z, float(eps))
        vals.append(chi2_bao(data_vec, cov, vec))

    return np.asarray(vals)


def best_result(curve):
    idx = int(np.argmin(curve))

    return {
        "index": idx,
        "eps": float(EPS_GRID[idx]),
        "chi2": float(curve[idx]),
        "boundary": idx == 0 or idx == len(EPS_GRID) - 1,
    }


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 78)
    print("AUDIT 06B — SDSS/BOSS DR12 BAO STATISTICAL IMPACT")
    print("=" * 78)

    print()
    print("Frozen repository:")
    print(FROZEN_ROOT)

    print()
    print("Fixed parameters:")
    print(f"H0      = {H0}")
    print(f"Om0     = {OM0}")
    print(f"n       = {N}")
    print(f"rd      = {RD} Mpc")
    print(f"rd_fid  = {RD_FID} Mpc")
    print(f"rd_fid / rd = {RD_FID / RD:.12f}")
    print(f"rd / rd_fid = {RD / RD_FID:.12f}")

    bao = load_bao_dr12(DATA_DIR)

    z = np.asarray(bao["z"], dtype=float)
    data_vec = np.asarray(bao["vec"], dtype=float)
    cov = np.asarray(bao["cov"], dtype=float)

    print()
    print("Data:")
    print(f"z = {z}")
    print("data vector:")
    print(data_vec)
    print(f"covariance shape = {cov.shape}")

    # --------------------------------------------------------
    # TEST 1 — exact frozen LCDM reproduction
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 1 — FROZEN LCDM REPRODUCTION")
    print("-" * 78)

    vec_lcdm_raw = vector_lcdm_frozen_raw(z)
    chi2_lcdm_raw = chi2_bao(data_vec, cov, vec_lcdm_raw)

    print("Frozen raw LCDM vector:")
    print(vec_lcdm_raw)
    print(f"chi2_LCDM_frozen_raw = {chi2_lcdm_raw:.12f}")

    reference = 11.14064185
    print(f"historical frozen reference = {reference:.12f}")
    print(
        "absolute difference       = "
        f"{abs(chi2_lcdm_raw - reference):.12e}"
    )

    pass_reproduce = abs(chi2_lcdm_raw - reference) < 1e-6
    print(f"RESULT: {'PASS' if pass_reproduce else 'FAIL'}")

    # --------------------------------------------------------
    # TEST 2 — convention correction at LCDM
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 2 — BOSS CONVENTION EFFECT AT LCDM")
    print("-" * 78)

    vec_lcdm_corr = vector_lcdm_corrected(z)
    chi2_lcdm_corr = chi2_bao(data_vec, cov, vec_lcdm_corr)

    print("Convention-corrected LCDM vector:")
    print(vec_lcdm_corr)

    print()
    print(f"chi2 frozen raw      = {chi2_lcdm_raw:.12f}")
    print(f"chi2 corrected       = {chi2_lcdm_corr:.12f}")
    print(
        "Delta chi2 convention = "
        f"{chi2_lcdm_corr - chi2_lcdm_raw:+.12f}"
    )

    # --------------------------------------------------------
    # TEST 3 — epsilon=0 identity
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 3 — EPSILON=0 IDENTITIES")
    print("-" * 78)

    vec_frozen_zero = vector_frozen_raw(z, 0.0)
    vec_frozen_corr_zero = vector_frozen_corrected(z, 0.0)
    vec_pub_corr_zero = vector_published_corrected(z, 0.0)

    diff_raw = np.max(np.abs(vec_frozen_zero - vec_lcdm_raw))
    diff_corr = np.max(np.abs(vec_frozen_corr_zero - vec_lcdm_corr))
    diff_pub = np.max(np.abs(vec_pub_corr_zero - vec_lcdm_corr))

    print(f"frozen Psi eps=0 vs frozen LCDM     : {diff_raw:.12e}")
    print(f"corrected frozen eps=0 vs corr LCDM : {diff_corr:.12e}")
    print(f"published eps=0 vs corr LCDM        : {diff_pub:.12e}")

    zero_pass = max(diff_raw, diff_corr, diff_pub) < 1e-5
    print(f"RESULT: {'PASS' if zero_pass else 'FAIL'}")

    # --------------------------------------------------------
    # TEST 4 — full epsilon scans
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 4 — EPSILON SCANS")
    print("-" * 78)

    curve_raw = scan_eps(
        data_vec,
        cov,
        z,
        vector_frozen_raw,
    )

    curve_frozen_corr = scan_eps(
        data_vec,
        cov,
        z,
        vector_frozen_corrected,
    )

    curve_pub_corr = scan_eps(
        data_vec,
        cov,
        z,
        vector_published_corrected,
    )

    best_raw = best_result(curve_raw)
    best_fc = best_result(curve_frozen_corr)
    best_pub = best_result(curve_pub_corr)

    idx_zero = int(np.argmin(np.abs(EPS_GRID)))

    base_raw = float(curve_raw[idx_zero])
    base_fc = float(curve_frozen_corr[idx_zero])
    base_pub = float(curve_pub_corr[idx_zero])

    print("FROZEN RAW:")
    print(f"  best eps       = {best_raw['eps']:+.6f}")
    print(f"  best chi2      = {best_raw['chi2']:.12f}")
    print(f"  chi2(eps=0)   = {base_raw:.12f}")
    print(f"  Delta chi2     = {best_raw['chi2'] - base_raw:+.12f}")
    print(f"  boundary       = {best_raw['boundary']}")

    print()
    print("FROZEN BACKGROUND + CORRECT BOSS CONVENTION:")
    print(f"  best eps       = {best_fc['eps']:+.6f}")
    print(f"  best chi2      = {best_fc['chi2']:.12f}")
    print(f"  chi2(eps=0)   = {base_fc:.12f}")
    print(f"  Delta chi2     = {best_fc['chi2'] - base_fc:+.12f}")
    print(f"  boundary       = {best_fc['boundary']}")

    print()
    print("PUBLISHED-v2 BACKGROUND + CORRECT BOSS CONVENTION:")
    print(f"  best eps       = {best_pub['eps']:+.6f}")
    print(f"  best chi2      = {best_pub['chi2']:.12f}")
    print(f"  chi2(eps=0)   = {base_pub:.12f}")
    print(f"  Delta chi2     = {best_pub['chi2'] - base_pub:+.12f}")
    print(f"  boundary       = {best_pub['boundary']}")

    # --------------------------------------------------------
    # TEST 5 — selected epsilon values
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 5 — SELECTED EPSILON VALUES")
    print("-" * 78)

    selected = [
        -0.10,
        -0.075,
        -0.05,
        -0.03,
        0.0,
        0.03,
        0.05,
        0.055,
        0.075,
        0.10,
    ]

    print(
        f"{'eps':>8} "
        f"{'frozen_raw':>16} "
        f"{'frozen_corr':>16} "
        f"{'published_corr':>18}"
    )

    for eps in selected:
        i = int(np.argmin(np.abs(EPS_GRID - eps)))

        print(
            f"{EPS_GRID[i]:+8.3f} "
            f"{curve_raw[i]:16.9f} "
            f"{curve_frozen_corr[i]:16.9f} "
            f"{curve_pub_corr[i]:18.9f}"
        )

    # --------------------------------------------------------
    # TEST 6 — isolate two effects at eps=0.05
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 6 — EFFECT DECOMPOSITION AT EPS=+0.05")
    print("-" * 78)

    eps_ref = 0.05
    i = int(np.argmin(np.abs(EPS_GRID - eps_ref)))

    raw = float(curve_raw[i])
    fc = float(curve_frozen_corr[i])
    pub = float(curve_pub_corr[i])

    print(f"frozen raw                       = {raw:.12f}")
    print(f"frozen + convention correction   = {fc:.12f}")
    print(f"published + convention correction= {pub:.12f}")

    print()
    print(
        "BAO convention contribution "
        f"(frozen corrected - raw) = {fc - raw:+.12f}"
    )

    print(
        "background-law contribution "
        f"(published corrected - frozen corrected) = {pub - fc:+.12f}"
    )

    print(
        "total raw->published corrected "
        f"= {pub - raw:+.12f}"
    )

    # Historical v2 value
    historical_eps05 = 10.63263995

    print()
    print(f"historical frozen BAO @ eps=.05 = {historical_eps05:.12f}")
    print(f"audit frozen raw @ eps=.05      = {raw:.12f}")
    print(
        "difference                         = "
        f"{raw - historical_eps05:+.12e}"
    )

    # --------------------------------------------------------
    # TEST 7 — vector-level comparison
    # --------------------------------------------------------

    print()
    print("-" * 78)
    print("TEST 7 — MODEL VECTOR COMPARISON AT EPS=+0.05")
    print("-" * 78)

    vraw = vector_frozen_raw(z, 0.05)
    vfc = vector_frozen_corrected(z, 0.05)
    vpub = vector_published_corrected(z, 0.05)

    print("data:")
    print(data_vec)

    print()
    print("frozen raw:")
    print(vraw)

    print()
    print("frozen convention-corrected:")
    print(vfc)

    print()
    print("published-v2 convention-corrected:")
    print(vpub)

    print()
    print("relative frozen convention correction:")
    print((vfc - vraw) / vraw)

    print()
    print("relative published-vs-frozen corrected:")
    print((vpub - vfc) / vfc)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("AUDIT 06B SUMMARY")
    print("=" * 78)

    print(
        "Frozen SDSS likelihood reproduction: "
        + ("PASS" if pass_reproduce else "FAIL")
    )

    print(
        "epsilon=0 consistency: "
        + ("PASS" if zero_pass else "FAIL")
    )

    same_best_fc = np.isclose(
        best_raw["eps"],
        best_fc["eps"],
        atol=1e-12,
    )

    same_best_pub = np.isclose(
        best_fc["eps"],
        best_pub["eps"],
        atol=1e-12,
    )

    print(
        "Raw frozen vs convention-corrected frozen "
        f"same best epsilon: {same_best_fc}"
    )

    print(
        "Corrected frozen vs corrected published "
        f"same best epsilon: {same_best_pub}"
    )

    print()
    print("Interpretation:")
    print(
        "1. The frozen SDSS DR12 numerical likelihood is first reproduced "
        "exactly before any correction is introduced."
    )
    print(
        "2. The BOSS rd_fid/rd convention effect is isolated from the "
        "PsiCDM background-law effect."
    )
    print(
        "3. The published-v2 curve changes only the background H(z) law "
        "and its derived DM(z); the data vector and covariance are unchanged."
    )
    print(
        "4. Any shift of the best epsilon between frozen and published "
        "curves therefore cannot be attributed to a different dataset "
        "or covariance."
    )


if __name__ == "__main__":
    main()
