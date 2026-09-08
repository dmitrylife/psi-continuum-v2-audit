#!/usr/bin/env python3
"""
Psi-Continuum v2 independent audit
===================================

AUDIT 08 — Joint likelihood reconstruction.

Purpose
-------
Reconstruct the complete joint epsilon scan from:

    Pantheon+SH0ES HF
    H(z) / Cosmic Chronometers
    SDSS/BOSS DR12 BAO
    DESI DR2 BAO

Two branches are evaluated.

A. FROZEN
   Exact reconstruction of the frozen-v2 numerical pipeline:
       - frozen PsiCDM background,
       - frozen SN likelihood,
       - frozen H(z) likelihood,
       - original SDSS raw [DM, H] convention,
       - frozen DESI observable construction.

B. PUBLISHED
   Keep the same datasets and likelihood definitions but replace the
   background relation by the literal published-v2 equation:

       H_Psi(z) = H_LCDM(z) * (1 + eps0/(1+z))

   For SDSS/BOSS DR12 use the corrected fiducial-rd convention established
   in Audit 06:

       DM * rd_fid / rd
       H  * rd / rd_fid

Important
---------
The PUBLISHED branch is a counterfactual audit reconstruction, not a new
cosmological inference.

In particular, its SN likelihood retains the raw frozen-v2 absolute
calibration treatment. Audit 04b showed that this treatment is highly
sensitive to the absolute SN normalization. Therefore the raw joint
published-spec minimum must not automatically be interpreted as a
statistically robust corrected-v2 result.

The frozen repository is never modified.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


# ======================================================================
# Paths
# ======================================================================

AUDIT_ROOT = Path(__file__).resolve().parent
V2_ROOT = (AUDIT_ROOT.parent / "psi-continuum-v2").resolve()
RESULTS_DIR = AUDIT_ROOT / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if not V2_ROOT.exists():
    raise RuntimeError(
        f"Frozen v2 repository not found:\n{V2_ROOT}"
    )

sys.path.insert(0, str(V2_ROOT))


# ======================================================================
# Frozen-v2 imports
# ======================================================================

from psi_continuum_v2.utils import get_data_path

from psi_continuum_v2.cosmology.data_loaders.pantheonplus_loader import (
    load_pantheonplus_hf,
)

from psi_continuum_v2.cosmology.data_loaders.hz_loader import (
    load_hz_compilation,
)

from psi_continuum_v2.cosmology.data_loaders.bao_loader import (
    load_bao_dr12,
)

from psi_continuum_v2.cosmology.data_loaders.desi_loader import (
    load_desi_dr2,
)

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
    mu_from_dL,
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

from psi_continuum_v2.cosmology.likelihoods.sn_likelihood import (
    chi2_sn_full_cov,
)

from psi_continuum_v2.cosmology.likelihoods.hz_likelihood import (
    chi2_hz,
)

from psi_continuum_v2.cosmology.likelihoods.bao_likelihood import (
    bao_vector_model,
    chi2_bao,
)

from psi_continuum_v2.cosmology.constants import (
    C_LIGHT,
)


# ======================================================================
# Configuration
# ======================================================================

H0 = 70.0
OM0 = 0.3
N_POWER = 1.0

RD = 147.0
RD_FID = 147.78

# Same epsilon grid as frozen joint scan.
EPS_GRID = np.linspace(-0.10, 0.10, 201)

# Independent integration grid for literal published equation.
SPEC_NZ = 200_001


# ======================================================================
# Historical frozen reference values
# ======================================================================

REFERENCE_SN_LCDM = 2609.19713024
REFERENCE_HZ_LCDM = 11.78688559
REFERENCE_SDSS_LCDM = 11.14064185
REFERENCE_DESI_LCDM = 40.29841997

REFERENCE_LCDM_TOTAL = 2672.42307766

REFERENCE_EPS005_TOTAL = 2666.25550047

REFERENCE_BEST_EPS = 0.075
REFERENCE_BEST_CHI2 = 2665.51086582


# ======================================================================
# Literal published-v2 background
# ======================================================================

def H_published(z, eps0):
    """
    Literal published-v2 defining equation:

        H_Psi(z) =
            H_LCDM(z) * (1 + eps0/(1+z))

    No additional z=0 normalization is introduced.
    """

    z = np.asarray(z, dtype=float)

    lcdm = LCDMParams(
        H0=H0,
        Om0=OM0,
        rd=RD,
    )

    H_lambda = np.asarray(
        H_lcdm(z, lcdm),
        dtype=float,
    )

    return H_lambda * (
        1.0 + float(eps0) / (1.0 + z)
    )


def DM_published(z_eval, eps0):
    """
    Independent flat-universe integration:

        DM(z) = c * integral_0^z dz' / H_published(z')

    The grid explicitly begins at z=0.
    """

    z_eval = np.asarray(z_eval, dtype=float)

    if z_eval.size == 0:
        return np.asarray([], dtype=float)

    zmax = float(np.max(z_eval))

    if zmax == 0.0:
        return np.zeros_like(z_eval)

    grid = np.linspace(
        0.0,
        zmax,
        SPEC_NZ,
    )

    H_grid = H_published(
        grid,
        eps0,
    )

    integral = cumulative_trapezoid(
        1.0 / H_grid,
        grid,
        initial=0.0,
    )

    DM_grid = C_LIGHT * integral

    return np.interp(
        z_eval,
        grid,
        DM_grid,
    )


def DH_published(z_eval, eps0):
    """
    Hubble distance:

        DH = c / H(z)
    """

    z_eval = np.asarray(z_eval, dtype=float)

    return (
        C_LIGHT
        / H_published(z_eval, eps0)
    )


def dL_published(z_eval, eps0):
    """
    Flat-universe luminosity distance:

        dL = (1+z) DM
    """

    z_eval = np.asarray(z_eval, dtype=float)

    return (
        (1.0 + z_eval)
        * DM_published(z_eval, eps0)
    )


# ======================================================================
# Generic covariance chi2
# ======================================================================

def chi2_cov(data, cov, model):
    """
    Generic Gaussian covariance chi2:

        chi2 = (d-m)^T C^{-1} (d-m)
    """

    data = np.asarray(data, dtype=float)
    cov = np.asarray(cov, dtype=float)
    model = np.asarray(model, dtype=float)

    diff = data - model

    return float(
        diff.T
        @ np.linalg.solve(cov, diff)
    )


# ======================================================================
# Pantheon+ SN
# ======================================================================

def sn_frozen(sn, eps0):
    """
    Exact frozen-v2 SN calculation.
    """

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=float(eps0),
        n=N_POWER,
        rd=RD,
    )

    dl = dL_psicdm(
        sn["z"],
        psi,
    )

    mu_model = mu_from_dL(dl)

    return float(
        chi2_sn_full_cov(
            sn["mu"],
            mu_model,
            sn["cov"],
        )
    )


def sn_published(sn, eps0):
    """
    Same frozen SN data/covariance/likelihood, but with the literal
    published-v2 background relation.
    """

    dl = dL_published(
        sn["z"],
        eps0,
    )

    mu_model = mu_from_dL(dl)

    return float(
        chi2_sn_full_cov(
            sn["mu"],
            mu_model,
            sn["cov"],
        )
    )


# ======================================================================
# H(z)
# ======================================================================

def hz_frozen(hz, eps0):
    """
    Exact frozen-v2 H(z) likelihood.

    Frozen API:
        chi2_hz(hzdata, H_model, params)
    """

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=float(eps0),
        n=N_POWER,
        rd=RD,
    )

    return float(
        chi2_hz(
            hz,
            H_psicdm,
            psi,
        )
    )


def hz_published(hz, eps0):
    """
    Frozen H(z) likelihood evaluated using the literal published equation.
    """

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=float(eps0),
        n=N_POWER,
        rd=RD,
    )

    def H_model(z, params):
        return H_published(
            z,
            params.eps0,
        )

    return float(
        chi2_hz(
            hz,
            H_model,
            psi,
        )
    )


# ======================================================================
# SDSS / BOSS DR12
# ======================================================================

def sdss_vector_frozen(z, eps0):
    """
    Exact frozen-v2 SDSS implementation:

        [DM, H, DM, H, DM, H]

    No rd/rdfid scaling.
    """

    z = np.asarray(z, dtype=float)

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=float(eps0),
        n=N_POWER,
        rd=RD,
    )

    dl = dL_psicdm(
        z,
        psi,
    )

    DM = dl / (1.0 + z)

    H = H_psicdm(
        z,
        psi,
    )

    return bao_vector_model(
        z,
        DM,
        H,
    )


def sdss_vector_published(z, eps0):
    """
    Literal published-v2 equation with the corrected BOSS DR12
    fiducial sound-horizon convention:

        DM * rd_fid / rd
        H  * rd / rd_fid
    """

    z = np.asarray(z, dtype=float)

    DM = DM_published(
        z,
        eps0,
    )

    H = H_published(
        z,
        eps0,
    )

    DM_scaled = (
        DM * RD_FID / RD
    )

    H_scaled = (
        H * RD / RD_FID
    )

    return bao_vector_model(
        z,
        DM_scaled,
        H_scaled,
    )


def sdss_frozen(bao, eps0):
    """
    Frozen SDSS likelihood.

    Frozen API:
        chi2_bao(data_vec, cov, model_vec)
    """

    model = sdss_vector_frozen(
        bao["z"],
        eps0,
    )

    return float(
        chi2_bao(
            bao["vec"],
            bao["cov"],
            model,
        )
    )


def sdss_published(bao, eps0):
    """
    Published equation + corrected BOSS convention.
    """

    model = sdss_vector_published(
        bao["z"],
        eps0,
    )

    return float(
        chi2_bao(
            bao["vec"],
            bao["cov"],
            model,
        )
    )


# ======================================================================
# DESI DR2
# ======================================================================

def desi_vector_frozen(z, labels, eps0):
    """
    Exact frozen DESI observable construction established in Audit 07.
    """

    z = np.asarray(z, dtype=float)

    psi = PsiCDMParams(
        H0=H0,
        Om0=OM0,
        eps0=float(eps0),
        n=N_POWER,
        rd=RD,
    )

    pred = []

    for zi, label in zip(z, labels):

        DM = float(
            DM_psicdm(
                zi,
                psi,
            )
        )

        DH = float(
            DH_psicdm(
                zi,
                psi,
            )
        )

        if label.startswith("DM"):

            value = DM / RD

        elif label.startswith("DH"):

            value = DH / RD

        elif label.startswith("DV"):

            DV = (
                DM
                * DM
                * zi
                * DH
            ) ** (1.0 / 3.0)

            value = DV / RD

        else:

            raise ValueError(
                f"Unknown DESI label: {label}"
            )

        pred.append(value)

    return np.asarray(
        pred,
        dtype=float,
    )


def desi_vector_published(z, labels, eps0):
    """
    DESI DR2 observables calculated from the literal published-v2
    background equation.
    """

    z = np.asarray(z, dtype=float)

    DM_all = DM_published(
        z,
        eps0,
    )

    DH_all = DH_published(
        z,
        eps0,
    )

    pred = []

    for i, (zi, label) in enumerate(
        zip(z, labels)
    ):

        DM = DM_all[i]
        DH = DH_all[i]

        if label.startswith("DM"):

            value = DM / RD

        elif label.startswith("DH"):

            value = DH / RD

        elif label.startswith("DV"):

            DV = (
                DM
                * DM
                * zi
                * DH
            ) ** (1.0 / 3.0)

            value = DV / RD

        else:

            raise ValueError(
                f"Unknown DESI label: {label}"
            )

        pred.append(value)

    return np.asarray(
        pred,
        dtype=float,
    )


def desi_frozen(desi, eps0):

    model = desi_vector_frozen(
        desi["z"],
        desi["labels"],
        eps0,
    )

    return chi2_cov(
        desi["vec"],
        desi["cov"],
        model,
    )


def desi_published(desi, eps0):

    model = desi_vector_published(
        desi["z"],
        desi["labels"],
        eps0,
    )

    return chi2_cov(
        desi["vec"],
        desi["cov"],
        model,
    )


# ======================================================================
# Data loading
# ======================================================================

def load_all_data():

    sn = load_pantheonplus_hf(
        get_data_path("pantheon_plus")
    )

    hz = load_hz_compilation(
        get_data_path("hz")
    )

    bao = load_bao_dr12(
        get_data_path("bao")
    )

    desi = load_desi_dr2(
        get_data_path("desi/dr2")
    )

    return (
        sn,
        hz,
        bao,
        desi,
    )


# ======================================================================
# Main
# ======================================================================

def main():

    print("=" * 86)
    print("AUDIT 08 — JOINT LIKELIHOOD RECONSTRUCTION")
    print("=" * 86)

    print()
    print("Frozen repository:")
    print(V2_ROOT)

    # ------------------------------------------------------------------
    # Load datasets
    # ------------------------------------------------------------------

    sn, hz, bao, desi = load_all_data()

    print()
    print("SETUP")
    print("-" * 86)

    print("SN N          :", len(sn["z"]))
    print("H(z) N        :", len(hz["z"]))
    print("SDSS vector N :", len(bao["vec"]))
    print("DESI vector N :", len(desi["vec"]))

    print("H0            :", H0)
    print("Om0           :", OM0)
    print("n             :", N_POWER)
    print("rd            :", RD)
    print("SDSS rd_fid   :", RD_FID)

    print(
        "eps grid       :",
        EPS_GRID[0],
        "...",
        EPS_GRID[-1],
    )

    print(
        "grid N         :",
        len(EPS_GRID),
    )

    # ------------------------------------------------------------------
    # Scan
    # ------------------------------------------------------------------

    rows = []

    print()
    print(
        f"Scanning {len(EPS_GRID)} epsilon values ..."
    )

    for k, eps0 in enumerate(EPS_GRID):

        # ==============================================================
        # A. Frozen branch
        # ==============================================================

        sn_f = sn_frozen(
            sn,
            eps0,
        )

        hz_f = hz_frozen(
            hz,
            eps0,
        )

        sdss_f = sdss_frozen(
            bao,
            eps0,
        )

        desi_f = desi_frozen(
            desi,
            eps0,
        )

        total_f = (
            sn_f
            + hz_f
            + sdss_f
            + desi_f
        )

        # ==============================================================
        # B. Published-spec branch
        # ==============================================================

        sn_p = sn_published(
            sn,
            eps0,
        )

        hz_p = hz_published(
            hz,
            eps0,
        )

        sdss_p = sdss_published(
            bao,
            eps0,
        )

        desi_p = desi_published(
            desi,
            eps0,
        )

        total_p = (
            sn_p
            + hz_p
            + sdss_p
            + desi_p
        )

        rows.append(
            [
                eps0,

                sn_f,
                hz_f,
                sdss_f,
                desi_f,
                total_f,

                sn_p,
                hz_p,
                sdss_p,
                desi_p,
                total_p,
            ]
        )

        if (
            k == 0
            or (k + 1) % 20 == 1
            or k == len(EPS_GRID) - 1
        ):
            print(
                f"{k+1:5d}/{len(EPS_GRID)} "
                f"eps={eps0:+.3f} "
                f"frozen={total_f:.6f} "
                f"published={total_p:.6f}"
            )

    arr = np.asarray(
        rows,
        dtype=float,
    )

    eps = arr[:, 0]

    total_f = arr[:, 5]
    total_p = arr[:, 10]

    # ------------------------------------------------------------------
    # Important indices
    # ------------------------------------------------------------------

    i0 = int(
        np.argmin(
            np.abs(eps)
        )
    )

    i003 = int(
        np.argmin(
            np.abs(eps - 0.03)
        )
    )

    i0031 = int(
        np.argmin(
            np.abs(eps - 0.031)
        )
    )

    i005 = int(
        np.argmin(
            np.abs(eps - 0.05)
        )
    )

    i0075 = int(
        np.argmin(
            np.abs(eps - 0.075)
        )
    )

    i_f = int(
        np.argmin(total_f)
    )

    i_p = int(
        np.argmin(total_p)
    )

    # ------------------------------------------------------------------
    # TEST 1
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 1 — FROZEN LCDM REPRODUCTION")
    print("=" * 86)

    print(
        "SN    :",
        f"{arr[i0,1]:.12f}",
    )

    print(
        "H(z)  :",
        f"{arr[i0,2]:.12f}",
    )

    print(
        "SDSS  :",
        f"{arr[i0,3]:.12f}",
    )

    print(
        "DESI  :",
        f"{arr[i0,4]:.12f}",
    )

    print(
        "TOTAL :",
        f"{arr[i0,5]:.12f}",
    )

    print()

    print(
        "reference SN    :",
        f"{REFERENCE_SN_LCDM:.12f}",
    )

    print(
        "reference H(z)  :",
        f"{REFERENCE_HZ_LCDM:.12f}",
    )

    print(
        "reference SDSS  :",
        f"{REFERENCE_SDSS_LCDM:.12f}",
    )

    print(
        "reference DESI  :",
        f"{REFERENCE_DESI_LCDM:.12f}",
    )

    print(
        "reference TOTAL :",
        f"{REFERENCE_LCDM_TOTAL:.12f}",
    )

    print()

    print(
        "SN difference   :",
        f"{arr[i0,1]-REFERENCE_SN_LCDM:+.12e}",
    )

    print(
        "H(z) difference :",
        f"{arr[i0,2]-REFERENCE_HZ_LCDM:+.12e}",
    )

    print(
        "SDSS difference :",
        f"{arr[i0,3]-REFERENCE_SDSS_LCDM:+.12e}",
    )

    print(
        "DESI difference :",
        f"{arr[i0,4]-REFERENCE_DESI_LCDM:+.12e}",
    )

    print(
        "TOTAL difference:",
        f"{arr[i0,5]-REFERENCE_LCDM_TOTAL:+.12e}",
    )

    # ------------------------------------------------------------------
    # TEST 2
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 2 — EPSILON = 0 NULL LIMIT")
    print("=" * 86)

    print(
        "Frozen total             :",
        f"{arr[i0,5]:.12f}",
    )

    print(
        "Published total          :",
        f"{arr[i0,10]:.12f}",
    )

    print()

    print(
        "Frozen SN                :",
        f"{arr[i0,1]:.12f}",
    )

    print(
        "Published SN             :",
        f"{arr[i0,6]:.12f}",
    )

    print(
        "SN difference            :",
        f"{arr[i0,6]-arr[i0,1]:+.12e}",
    )

    print()

    print(
        "Frozen H(z)              :",
        f"{arr[i0,2]:.12f}",
    )

    print(
        "Published H(z)           :",
        f"{arr[i0,7]:.12f}",
    )

    print(
        "H(z) difference          :",
        f"{arr[i0,7]-arr[i0,2]:+.12e}",
    )

    print()

    print(
        "Frozen SDSS              :",
        f"{arr[i0,3]:.12f}",
    )

    print(
        "Published corrected SDSS :",
        f"{arr[i0,8]:.12f}",
    )

    print(
        "SDSS convention shift    :",
        f"{arr[i0,8]-arr[i0,3]:+.12f}",
    )

    print()

    print(
        "Frozen DESI              :",
        f"{arr[i0,4]:.12f}",
    )

    print(
        "Published DESI           :",
        f"{arr[i0,9]:.12f}",
    )

    print(
        "DESI difference          :",
        f"{arr[i0,9]-arr[i0,4]:+.12e}",
    )

    print()
    print("NOTE:")
    print(
        "Published total is not expected to equal frozen total at eps=0,"
    )
    print(
        "because the published branch deliberately applies the corrected"
    )
    print(
        "SDSS/BOSS rd_fid convention."
    )

    # ------------------------------------------------------------------
    # TEST 3
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 3 — FROZEN EPSILON = +0.05 REPRODUCTION")
    print("=" * 86)

    print(
        "SN    :",
        f"{arr[i005,1]:.12f}",
    )

    print(
        "H(z)  :",
        f"{arr[i005,2]:.12f}",
    )

    print(
        "SDSS  :",
        f"{arr[i005,3]:.12f}",
    )

    print(
        "DESI  :",
        f"{arr[i005,4]:.12f}",
    )

    print(
        "TOTAL :",
        f"{arr[i005,5]:.12f}",
    )

    print()

    print(
        "reference total :",
        f"{REFERENCE_EPS005_TOTAL:.12f}",
    )

    print(
        "difference      :",
        f"{arr[i005,5]-REFERENCE_EPS005_TOTAL:+.12e}",
    )

    # ------------------------------------------------------------------
    # TEST 4
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 4 — FROZEN JOINT MINIMUM")
    print("=" * 86)

    frozen_delta = (
        total_f[i_f]
        - total_f[i0]
    )

    print(
        "best eps  :",
        f"{eps[i_f]:+.6f}",
    )

    print(
        "best chi2 :",
        f"{total_f[i_f]:.12f}",
    )

    print(
        "delta chi2:",
        f"{frozen_delta:+.12f}",
    )

    print(
        "boundary  :",
        bool(
            i_f == 0
            or i_f == len(eps) - 1
        ),
    )

    print()

    print(
        "reference best eps  :",
        f"{REFERENCE_BEST_EPS:+.6f}",
    )

    print(
        "reference best chi2 :",
        f"{REFERENCE_BEST_CHI2:.12f}",
    )

    print()

    print(
        "SN    :",
        f"{arr[i_f,1]:.12f}",
    )

    print(
        "H(z)  :",
        f"{arr[i_f,2]:.12f}",
    )

    print(
        "SDSS  :",
        f"{arr[i_f,3]:.12f}",
    )

    print(
        "DESI  :",
        f"{arr[i_f,4]:.12f}",
    )

    # ------------------------------------------------------------------
    # TEST 5
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 5 — PUBLISHED-v2 JOINT MINIMUM")
    print("=" * 86)

    published_delta = (
        total_p[i_p]
        - total_p[i0]
    )

    print(
        "best eps  :",
        f"{eps[i_p]:+.6f}",
    )

    print(
        "best chi2 :",
        f"{total_p[i_p]:.12f}",
    )

    print(
        "delta chi2:",
        f"{published_delta:+.12f}",
    )

    print(
        "boundary  :",
        bool(
            i_p == 0
            or i_p == len(eps) - 1
        ),
    )

    print()

    print(
        "SN    :",
        f"{arr[i_p,6]:.12f}",
    )

    print(
        "H(z)  :",
        f"{arr[i_p,7]:.12f}",
    )

    print(
        "SDSS  :",
        f"{arr[i_p,8]:.12f}",
    )

    print(
        "DESI  :",
        f"{arr[i_p,9]:.12f}",
    )

    # ------------------------------------------------------------------
    # TEST 6
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 6 — SELECTED EPSILON VALUES")
    print("=" * 86)

    print(
        f"{'eps':>8}"
        f"{'frozen total':>18}"
        f"{'dchi2 F':>16}"
        f"{'published total':>20}"
        f"{'dchi2 P':>16}"
    )

    selected = [
        -0.100,
        -0.050,
        -0.030,
        -0.022,
        0.000,
        0.030,
        0.031,
        0.050,
        0.075,
        0.100,
    ]

    for target in selected:

        j = int(
            np.argmin(
                np.abs(
                    eps - target
                )
            )
        )

        print(
            f"{eps[j]:+8.3f}"
            f"{total_f[j]:18.8f}"
            f"{total_f[j]-total_f[i0]:+16.8f}"
            f"{total_p[j]:20.8f}"
            f"{total_p[j]-total_p[i0]:+16.8f}"
        )

    # ------------------------------------------------------------------
    # TEST 7
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 7 — DATASET CONTRIBUTIONS TO JOINT DELTA CHI2")
    print("=" * 86)

    print()

    print(
        "FROZEN best epsilon:",
        f"{eps[i_f]:+.6f}",
    )

    frozen_names = [
        ("SN", 1),
        ("H(z)", 2),
        ("SDSS", 3),
        ("DESI", 4),
    ]

    frozen_sum = 0.0

    for name, col in frozen_names:

        delta = (
            arr[i_f, col]
            - arr[i0, col]
        )

        frozen_sum += delta

        print(
            f"{name:6s}: "
            f"{delta:+.12f}"
        )

    print(
        "SUM   :",
        f"{frozen_sum:+.12f}",
    )

    print(
        "TOTAL :",
        f"{frozen_delta:+.12f}",
    )

    print(
        "sum error:",
        f"{frozen_sum-frozen_delta:+.12e}",
    )

    print()

    print(
        "PUBLISHED best epsilon:",
        f"{eps[i_p]:+.6f}",
    )

    published_names = [
        ("SN", 6),
        ("H(z)", 7),
        ("SDSS", 8),
        ("DESI", 9),
    ]

    published_sum = 0.0

    for name, col in published_names:

        delta = (
            arr[i_p, col]
            - arr[i0, col]
        )

        published_sum += delta

        print(
            f"{name:6s}: "
            f"{delta:+.12f}"
        )

    print(
        "SUM   :",
        f"{published_sum:+.12f}",
    )

    print(
        "TOTAL :",
        f"{published_delta:+.12f}",
    )

    print(
        "sum error:",
        f"{published_sum-published_delta:+.12e}",
    )

    # ------------------------------------------------------------------
    # TEST 8
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 8 — DIRECT COMPARISON AT EPSILON = +0.05")
    print("=" * 86)

    print(
        f"{'dataset':<10}"
        f"{'frozen':>18}"
        f"{'published':>18}"
        f"{'P-F':>18}"
    )

    compare = [
        ("SN", 1, 6),
        ("H(z)", 2, 7),
        ("SDSS", 3, 8),
        ("DESI", 4, 9),
        ("TOTAL", 5, 10),
    ]

    for name, col_f, col_p in compare:

        fval = arr[i005, col_f]
        pval = arr[i005, col_p]

        print(
            f"{name:<10}"
            f"{fval:18.8f}"
            f"{pval:18.8f}"
            f"{pval-fval:+18.8f}"
        )

    # ------------------------------------------------------------------
    # TEST 9 — exact total checks
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("TEST 9 — JOINT TOTAL = SUM OF COMPONENTS")
    print("=" * 86)

    reconstructed_f = (
        arr[:, 1]
        + arr[:, 2]
        + arr[:, 3]
        + arr[:, 4]
    )

    reconstructed_p = (
        arr[:, 6]
        + arr[:, 7]
        + arr[:, 8]
        + arr[:, 9]
    )

    frozen_sum_error = float(
        np.max(
            np.abs(
                reconstructed_f
                - arr[:, 5]
            )
        )
    )

    published_sum_error = float(
        np.max(
            np.abs(
                reconstructed_p
                - arr[:, 10]
            )
        )
    )

    print(
        "Frozen max absolute sum error   :",
        f"{frozen_sum_error:.12e}",
    )

    print(
        "Published max absolute sum error:",
        f"{published_sum_error:.12e}",
    )

    # ------------------------------------------------------------------
    # Historical frozen control points
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("HISTORICAL FROZEN CONTROL POINTS")
    print("=" * 86)

    print(
        "eps = 0.000 :",
        f"{total_f[i0]:.12f}",
    )

    print(
        "eps = 0.030 :",
        f"{total_f[i003]:.12f}",
    )

    print(
        "eps = 0.031 :",
        f"{total_f[i0031]:.12f}",
    )

    print(
        "eps = 0.050 :",
        f"{total_f[i005]:.12f}",
    )

    print(
        "eps = 0.075 :",
        f"{total_f[i0075]:.12f}",
    )

    # ------------------------------------------------------------------
    # Save CSV
    # ------------------------------------------------------------------

    out_csv = (
        RESULTS_DIR
        / "08_joint_scan.csv"
    )

    header = (
        "eps,"
        "sn_frozen,"
        "hz_frozen,"
        "sdss_frozen,"
        "desi_frozen,"
        "total_frozen,"
        "sn_published,"
        "hz_published,"
        "sdss_published,"
        "desi_published,"
        "total_published"
    )

    np.savetxt(
        out_csv,
        arr,
        delimiter=",",
        header=header,
        comments="",
        fmt="%.12f",
    )

    print()
    print("CSV saved:")
    print(out_csv)

    # ------------------------------------------------------------------
    # Final audit checks
    # ------------------------------------------------------------------

    frozen_lcdm_ok = bool(
        abs(
            total_f[i0]
            - REFERENCE_LCDM_TOTAL
        ) < 1e-5
    )

    frozen_005_ok = bool(
        abs(
            total_f[i005]
            - REFERENCE_EPS005_TOTAL
        ) < 1e-5
    )

    frozen_best_eps_ok = bool(
        np.isclose(
            eps[i_f],
            REFERENCE_BEST_EPS,
            atol=1e-12,
            rtol=0.0,
        )
    )

    frozen_best_chi2_ok = bool(
        abs(
            total_f[i_f]
            - REFERENCE_BEST_CHI2
        ) < 1e-5
    )

    exact_sum_frozen = bool(
        frozen_sum_error < 1e-10
    )

    exact_sum_published = bool(
        published_sum_error < 1e-10
    )

    frozen_curve = (
        total_f
        - total_f[i0]
    )

    published_curve = (
        total_p
        - total_p[i0]
    )

    same_curve = bool(
        np.allclose(
            frozen_curve,
            published_curve,
            rtol=1e-7,
            atol=1e-7,
        )
    )

    same_best = bool(
        np.isclose(
            eps[i_f],
            eps[i_p],
            atol=1e-12,
            rtol=0.0,
        )
    )

    frozen_pass = bool(
        frozen_lcdm_ok
        and frozen_005_ok
        and frozen_best_eps_ok
        and frozen_best_chi2_ok
        and exact_sum_frozen
    )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print()
    print("=" * 86)
    print("AUDIT 08 SUMMARY")
    print("=" * 86)

    print()

    print(
        "Frozen eps=0 total       :",
        f"{total_f[i0]:.12f}",
    )

    print(
        "Published eps=0 total    :",
        f"{total_p[i0]:.12f}",
    )

    print()

    print(
        "Frozen best epsilon      :",
        f"{eps[i_f]:+.6f}",
    )

    print(
        "Frozen best chi2         :",
        f"{total_f[i_f]:.12f}",
    )

    print(
        "Frozen delta chi2        :",
        f"{frozen_delta:+.12f}",
    )

    print(
        "Frozen boundary          :",
        bool(
            i_f == 0
            or i_f == len(eps) - 1
        ),
    )

    print()

    print(
        "Published best epsilon   :",
        f"{eps[i_p]:+.6f}",
    )

    print(
        "Published best chi2      :",
        f"{total_p[i_p]:.12f}",
    )

    print(
        "Published delta chi2     :",
        f"{published_delta:+.12f}",
    )

    print(
        "Published boundary       :",
        bool(
            i_p == 0
            or i_p == len(eps) - 1
        ),
    )

    print()

    print(
        "Frozen LCDM reproduction :",
        frozen_lcdm_ok,
    )

    print(
        "Frozen eps=.05 reproduce :",
        frozen_005_ok,
    )

    print(
        "Frozen best eps reproduce:",
        frozen_best_eps_ok,
    )

    print(
        "Frozen best chi2 reproduce:",
        frozen_best_chi2_ok,
    )

    print(
        "Frozen exact-sum check   :",
        exact_sum_frozen,
    )

    print(
        "Published exact-sum check:",
        exact_sum_published,
    )

    print()

    print(
        "Same best epsilon        :",
        same_best,
    )

    print(
        "Same delta-chi2 curve    :",
        same_curve,
    )

    print()

    if frozen_pass:
        frozen_status = "PASS"
    else:
        frozen_status = "FAIL"

    if same_best and same_curve:
        comparison_status = "PASS"
    else:
        comparison_status = "MISMATCH"

    print(
        "Frozen reproduction      :",
        frozen_status,
    )

    print(
        "Published-spec comparison:",
        comparison_status,
    )

    print()

    print(
        "FINAL STATUS: "
        f"FROZEN REPRODUCTION {frozen_status} / "
        f"PUBLISHED-SPEC JOINT COMPARISON {comparison_status}"
    )


if __name__ == "__main__":
    main()
