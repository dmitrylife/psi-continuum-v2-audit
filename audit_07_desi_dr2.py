from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


# ============================================================
# Frozen repository
# ============================================================

AUDIT_ROOT = Path(__file__).resolve().parent
FROZEN_ROOT = (AUDIT_ROOT / "../psi-continuum-v2").resolve()

sys.path.insert(0, str(FROZEN_ROOT))

from psi_continuum_v2.cosmology.data_loaders.desi_loader import load_desi_dr2
from psi_continuum_v2.cosmology.models.lcdm_params import LCDMParams
from psi_continuum_v2.cosmology.models.psicdm_params import PsiCDMParams

from psi_continuum_v2.cosmology.background.lcdm import (
    H_lcdm,
    DM_lcdm,
    DH_lcdm,
)

from psi_continuum_v2.cosmology.background.psicdm import (
    H_psicdm,
    DM_psicdm,
    DH_psicdm,
)

from psi_continuum_v2.cosmology.constants import C_LIGHT


# ============================================================
# Configuration
# ============================================================

H0 = 70.0
OM0 = 0.3
N = 1.0
RD = 147.0

EPS_GRID = np.linspace(-0.10, 0.10, 201)

DESI_DIR = FROZEN_ROOT / "data/desi/dr2"


# ============================================================
# Generic chi2
# ============================================================

def chi2(data, cov, model):
    diff = np.asarray(data) - np.asarray(model)
    return float(diff.T @ np.linalg.solve(cov, diff))


# ============================================================
# Frozen DESI prediction
# ============================================================

def frozen_vector(z, labels, eps=None):
    if eps is None:
        params = LCDMParams(H0=H0, Om0=OM0)
        model = "lcdm"
        rd = params.rd
    else:
        params = PsiCDMParams(
            H0=H0,
            Om0=OM0,
            eps0=float(eps),
            n=N,
            rd=RD,
        )
        model = "psicdm"
        rd = params.rd

    pred = []

    for zi, label in zip(z, labels):
        if model == "lcdm":
            DM = float(DM_lcdm(zi, params))
            DH = float(DH_lcdm(zi, params))
        else:
            DM = float(DM_psicdm(zi, params))
            DH = float(DH_psicdm(zi, params))

        if label.startswith("DM"):
            value = DM / rd

        elif label.startswith("DH"):
            value = DH / rd

        elif label.startswith("DV"):
            DV = (DM * DM * zi * DH) ** (1.0 / 3.0)
            value = DV / rd

        else:
            raise ValueError(label)

        pred.append(value)

    return np.asarray(pred)


# ============================================================
# Published-v2 defining equation
#
# H_Psi(z) = H_LCDM(z) * (1 + eps/(1+z))
#
# IMPORTANT:
# no extra z=0 renormalisation is introduced here.
# ============================================================

def H_lcdm_independent(z):
    z = np.asarray(z, dtype=float)
    return H0 * np.sqrt(
        OM0 * (1.0 + z) ** 3 +
        (1.0 - OM0)
    )


def H_published(z, eps):
    z = np.asarray(z, dtype=float)
    return H_lcdm_independent(z) * (
        1.0 + float(eps) / (1.0 + z)
    )


def DM_from_H(z_target, Hfunc):
    """
    Independent flat-universe integration:
        DM(z) = c * integral_0^z dz'/H(z')
    """

    z_target = np.asarray(z_target, dtype=float)

    zmax = float(np.max(z_target))

    # Dense integration grid including z=0 exactly.
    grid = np.linspace(0.0, zmax, 100001)

    H_grid = Hfunc(grid)

    integral = cumulative_trapezoid(
        1.0 / H_grid,
        grid,
        initial=0.0,
    )

    DM_grid = C_LIGHT * integral

    return np.interp(z_target, grid, DM_grid)


def published_vector(z, labels, eps):
    z = np.asarray(z, dtype=float)

    H = H_published(z, eps)

    DM_all = DM_from_H(
        z,
        lambda zz: H_published(zz, eps),
    )

    DH_all = C_LIGHT / H

    pred = []

    for i, (zi, label) in enumerate(zip(z, labels)):
        DM = DM_all[i]
        DH = DH_all[i]

        if label.startswith("DM"):
            value = DM / RD

        elif label.startswith("DH"):
            value = DH / RD

        elif label.startswith("DV"):
            DV = (DM * DM * zi * DH) ** (1.0 / 3.0)
            value = DV / RD

        else:
            raise ValueError(label)

        pred.append(value)

    return np.asarray(pred)


# ============================================================
# Main audit
# ============================================================

def main():

    print("=" * 78)
    print("AUDIT 07 — DESI DR2 BAO")
    print("=" * 78)

    print()
    print("Frozen repository:")
    print(FROZEN_ROOT)

    desi = load_desi_dr2(DESI_DIR)

    z = np.asarray(desi["z"])
    labels = list(desi["labels"])
    obs = np.asarray(desi["vec"])
    cov = np.asarray(desi["cov"])

    print()
    print("SETUP")
    print("-" * 78)
    print("N observations :", len(obs))
    print("cov shape      :", cov.shape)
    print("H0             :", H0)
    print("Om0            :", OM0)
    print("n              :", N)
    print("rd             :", RD)
    print("eps grid       :", EPS_GRID[0], "...", EPS_GRID[-1])
    print("grid N         :", len(EPS_GRID))

    # --------------------------------------------------------
    # Test 1: covariance
    # --------------------------------------------------------

    print()
    print("TEST 1 — covariance diagnostics")
    print("-" * 78)

    sym_err = np.max(np.abs(cov - cov.T))
    eig = np.linalg.eigvalsh(cov)

    print("max |C-C^T|    :", f"{sym_err:.12e}")
    print("min eigenvalue :", f"{eig.min():.12e}")
    print("max eigenvalue :", f"{eig.max():.12e}")
    print("positive def.  :", bool(np.all(eig > 0)))

    cov_pass = sym_err < 1e-12 and np.all(eig > 0)

    print("STATUS          :", "PASS" if cov_pass else "FAIL")

    # --------------------------------------------------------
    # Test 2: frozen LCDM reproduction
    # --------------------------------------------------------

    print()
    print("TEST 2 — frozen LCDM")
    print("-" * 78)

    vec_lcdm = frozen_vector(z, labels, eps=None)
    chi_lcdm = chi2(obs, cov, vec_lcdm)

    reference = 40.29841997

    print("chi2 LCDM      :", f"{chi_lcdm:.12f}")
    print("reference      :", f"{reference:.12f}")
    print("difference     :", f"{chi_lcdm-reference:+.12e}")

    lcdm_pass = abs(chi_lcdm - reference) < 1e-6

    print("STATUS         :", "PASS" if lcdm_pass else "FAIL")

    # --------------------------------------------------------
    # Test 3: null limit
    # --------------------------------------------------------

    print()
    print("TEST 3 — eps=0 null limit")
    print("-" * 78)

    vec_frozen_zero = frozen_vector(z, labels, eps=0.0)
    vec_pub_zero = published_vector(z, labels, eps=0.0)

    chi_frozen_zero = chi2(obs, cov, vec_frozen_zero)
    chi_pub_zero = chi2(obs, cov, vec_pub_zero)

    print("LCDM chi2          :", f"{chi_lcdm:.12f}")
    print("frozen eps=0 chi2  :", f"{chi_frozen_zero:.12f}")
    print("published eps=0    :", f"{chi_pub_zero:.12f}")

    print(
        "max vector diff frozen/LCDM:",
        f"{np.max(np.abs(vec_frozen_zero-vec_lcdm)):.12e}",
    )

    print(
        "max vector diff pub/LCDM   :",
        f"{np.max(np.abs(vec_pub_zero-vec_lcdm)):.12e}",
    )

    null_frozen = np.allclose(
        vec_frozen_zero, vec_lcdm,
        rtol=1e-7, atol=1e-8
    )

    null_pub = np.allclose(
        vec_pub_zero, vec_lcdm,
        rtol=1e-6, atol=1e-6
    )

    print("frozen null STATUS :", "PASS" if null_frozen else "FAIL")
    print("pub null STATUS    :", "PASS" if null_pub else "FAIL")

    # --------------------------------------------------------
    # Test 4: scans
    # --------------------------------------------------------

    print()
    print("TEST 4 — DESI-only epsilon scans")
    print("-" * 78)

    chi_frozen = []
    chi_pub = []

    for eps in EPS_GRID:
        vf = frozen_vector(z, labels, eps)
        vp = published_vector(z, labels, eps)

        chi_frozen.append(chi2(obs, cov, vf))
        chi_pub.append(chi2(obs, cov, vp))

    chi_frozen = np.asarray(chi_frozen)
    chi_pub = np.asarray(chi_pub)

    i_f = int(np.argmin(chi_frozen))
    i_p = int(np.argmin(chi_pub))

    eps_f = EPS_GRID[i_f]
    eps_p = EPS_GRID[i_p]

    print("FROZEN:")
    print("  best eps       :", f"{eps_f:+.6f}")
    print("  best chi2      :", f"{chi_frozen[i_f]:.12f}")
    print(
        "  delta chi2     :",
        f"{chi_frozen[i_f]-chi_lcdm:+.12f}",
    )
    print(
        "  boundary       :",
        bool(i_f == 0 or i_f == len(EPS_GRID)-1),
    )

    print()
    print("PUBLISHED EQUATION:")
    print("  best eps       :", f"{eps_p:+.6f}")
    print("  best chi2      :", f"{chi_pub[i_p]:.12f}")
    print(
        "  delta chi2     :",
        f"{chi_pub[i_p]-chi_lcdm:+.12f}",
    )
    print(
        "  boundary       :",
        bool(i_p == 0 or i_p == len(EPS_GRID)-1),
    )

    # --------------------------------------------------------
    # Test 5: selected eps
    # --------------------------------------------------------

    print()
    print("TEST 5 — selected epsilon values")
    print("-" * 78)

    selected = [-0.10, -0.05, 0.0, 0.03, 0.05, 0.075, 0.10]

    print(
        f"{'eps':>8} "
        f"{'chi2 frozen':>16} "
        f"{'delta frozen':>16} "
        f"{'chi2 published':>18} "
        f"{'delta published':>18}"
    )

    for eps in selected:
        vf = frozen_vector(z, labels, eps)
        vp = published_vector(z, labels, eps)

        cf = chi2(obs, cov, vf)
        cp = chi2(obs, cov, vp)

        print(
            f"{eps:+8.3f} "
            f"{cf:16.8f} "
            f"{cf-chi_lcdm:+16.8f} "
            f"{cp:18.8f} "
            f"{cp-chi_lcdm:+18.8f}"
        )

    # --------------------------------------------------------
    # Test 6: model-vector discrepancy
    # --------------------------------------------------------

    print()
    print("TEST 6 — frozen vs published observable vectors")
    print("-" * 78)

    for eps in [0.03, 0.05, 0.075, -0.05]:
        vf = frozen_vector(z, labels, eps)
        vp = published_vector(z, labels, eps)

        rel = np.abs(vf - vp) / np.maximum(np.abs(vp), 1e-30)

        print(
            f"eps={eps:+.3f}  "
            f"max_rel={np.max(rel):.8e}  "
            f"rms_rel={np.sqrt(np.mean(rel**2)):.8e}"
        )

    # --------------------------------------------------------
    # Test 7: verify old DESI test interpretation
    # --------------------------------------------------------

    print()
    print("TEST 7 — legacy DESI test interpretation")
    print("-" * 78)

    print("Legacy test used eps0 = 0.0")
    print("LCDM chi2         :", f"{chi_lcdm:.12f}")
    print("Psi eps=0 chi2    :", f"{chi_frozen_zero:.12f}")
    print(
        "Delta            :",
        f"{chi_frozen_zero-chi_lcdm:+.12e}",
    )

    legacy_pass = abs(chi_frozen_zero - chi_lcdm) < 1e-8

    print(
        "Interpretation    :",
        "null-limit equality confirmed"
        if legacy_pass else
        "unexpected discrepancy"
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 78)
    print("AUDIT 07 SUMMARY")
    print("=" * 78)

    print(
        "DESI data/covariance integrity :",
        "PASS" if cov_pass else "FAIL",
    )

    print(
        "Frozen LCDM reproduction       :",
        "PASS" if lcdm_pass else "FAIL",
    )

    print(
        "Frozen eps=0 null limit         :",
        "PASS" if null_frozen else "FAIL",
    )

    print(
        "Published eps=0 null limit      :",
        "PASS" if null_pub else "FAIL",
    )

    same_best = abs(eps_f - eps_p) < 1e-12
    same_curve = np.allclose(
        chi_frozen,
        chi_pub,
        rtol=1e-7,
        atol=1e-7,
    )

    print()
    print("Frozen best eps                 :", f"{eps_f:+.6f}")
    print("Published-equation best eps     :", f"{eps_p:+.6f}")
    print("Same best epsilon               :", same_best)
    print("Same chi2 curve                 :", same_curve)

    print()
    print(
        "Frozen DESI delta chi2          :",
        f"{chi_frozen[i_f]-chi_lcdm:+.12f}",
    )

    print(
        "Published DESI delta chi2       :",
        f"{chi_pub[i_p]-chi_lcdm:+.12f}",
    )

    print()
    print(
        "FINAL STATUS: FROZEN REPRODUCTION "
        + ("PASS" if lcdm_pass else "FAIL")
        + " / PUBLISHED-SPEC COMPARISON COMPLETED"
    )


if __name__ == "__main__":
    main()
