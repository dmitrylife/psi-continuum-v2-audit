# Ψ-Continuum v2 --- Independent Scientific Code Audit

Independent reproducibility and model-consistency audit of the frozen
**Ψ-Continuum v2 (v0.2.3)** numerical implementation.

> **Main finding:** the frozen v2 numerical pipeline is internally
> reproducible, but the non-zero-$\varepsilon_0$ ΨCDM
> background implemented in the frozen code is not mathematically
> equivalent to the defining equation documented for v2.

The original `psi-continuum-v2` repository is treated as a **frozen
scientific artifact** and is not modified by this audit.

## Scope

The audit checks the full numerical chain used by v2:

-   ΨCDM background implementation;
-   implementation/dependency trace;
-   cosmological distances;
-   Pantheon+ supernova likelihood;
-   supernova calibration sensitivity;
-   $H(z)$ / Cosmic Chronometers;
-   SDSS/BOSS DR12 BAO;
-   DESI DR2 BAO;
-   complete joint likelihood reconstruction.

The audit distinguishes three questions: **reproducibility** of frozen
results, **model consistency** between code and the documented equation,
and **likelihood consistency** after identified
implementation/convention differences.

## Defining equation under audit

The documented v2 background relation is


$$
H_\Psi(z)=H_\Lambda(z)\left(1+\frac{\varepsilon_0}{1+z}\right).
$$


The frozen implementation instead uses a normalized modified-$E^2$
prescription that enforces $H_\Psi(0)=H_0$ for every
$\varepsilon_0$. Both formulations recover the same ΛCDM limit
at $\varepsilon_0=0$, but they are not equivalent for non-zero
$\varepsilon_0$.

## Main results

### 1. Frozen numerical pipeline is reproducible

The individual frozen likelihoods and the complete joint scan are
reproduced to numerical precision. For the joint ΛCDM point,


$$
\chi^2_{\mathrm{total}}=2672.423077655228,
$$


compared with the historical frozen reference (2672.423077660000).

The reconstructed frozen joint minimum is


$$
\varepsilon_{0,\mathrm{best}}=+0.075,\qquad
\chi^2_{\mathrm{min}}=2665.510865815413,\qquad
\Delta\chi^2=-6.912211839815.
$$


Thus the audit does **not** find a general numerical reproducibility
failure in the frozen software.

### 2. Published equation and frozen implementation differ

For non-zero $\varepsilon_0$, the implemented background
expansion differs from the documented v2 equation. This discrepancy
propagates into cosmological distances and materially changes the
likelihood curves for Pantheon+ SN, $H(z)$, SDSS/BOSS DR12, DESI DR2,
and the joint analysis.

Consequently, non-zero-$\varepsilon_0$ preferences obtained
from the frozen pipeline describe the **implemented frozen model**, not
the documented v2 defining equation.

### 3. SDSS/BOSS DR12 has an additional convention mismatch

The frozen SDSS/BOSS model vector does not explicitly apply the standard
fiducial sound-horizon factors associated with the DR12 consensus
observables. The frozen code uses $r_d=147.0$ Mpc, whereas the fiducial
value associated with the BOSS DR12 vector is
$r_{d,\mathrm{fid}}=147.78$ Mpc.

This produces a smaller but genuine and numerically non-negligible BAO
convention mismatch.

### 4. Raw Pantheon+ comparison is calibration-sensitive

Substituting the literal documented equation into the frozen raw SN
likelihood produces a very large apparent improvement in
$\chi^2$. A diagnostic profiling of one global additive
distance-modulus offset removes approximately **98.46%** of this raw
improvement.

Therefore the large raw SN and raw published-spec joint improvements
must not be interpreted directly as evidence for the documented ΨCDM
equation.

## Audit status

  -----------------------------------------------------------------------
  Audit                   Test                    Result
  ----------------------- ----------------------- -----------------------
  01                      ΨCDM background         **FAIL ---
                                                  model/implementation
                                                  mismatch**

  02                      Implementation trace    **TRACE CONFIRMED**

  03                      Distances               **INTERNAL PASS /
                                                  PUBLISHED-SPEC FAIL**

  04                      Pantheon+ SN            **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  04b                     SN calibration          **CALIBRATION EFFECT
                          diagnostic              CONFIRMED / BEST-FIT
                                                  UNRESOLVED**

  05                      $H(z)$                  **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  06                      SDSS/BOSS DR12          **FROZEN PASS /
                                                  CONVENTION +
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  07                      DESI DR2                **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  08                      Joint likelihood        **FROZEN PASS /
                                                  PUBLISHED-SPEC JOINT
                                                  MISMATCH**
  -----------------------------------------------------------------------

## Repository structure

``` text
psi-continuum-v2-audit/
├── README.md
├── AUDIT_REPORT.md
├── audit_01_model.py
├── audit_02_trace.py
├── audit_03_distances.py
├── audit_04_sn.py
├── audit_04b_sn_calibration.py
├── audit_05_hz.py
├── audit_06_bao_sdss.py
├── audit_07_desi_dr2.py
├── audit_08_joint.py
└── results/
    ├── 01_model.txt
    ├── 02_trace.txt
    ├── 03_distances.txt
    ├── 04_sn.txt
    ├── 04_sn_scan.csv
    ├── 04b_sn_calibration.txt
    ├── 04b_sn_calibration_scan.csv
    ├── 05_hz.txt
    ├── 05_hz_scan.csv
    ├── audit_06_bao_sdss.txt
    ├── 07_desi_dr2.txt
    ├── 08_joint.txt
    └── 08_joint_scan.csv
```

`AUDIT_REPORT.md` contains the complete scientific interpretation and
numerical results. The individual scripts and files under `results/`
provide the reproducibility trail.

## Relationship to frozen v2

This repository does **not** replace or silently correct
`psi-continuum-v2`.

``` text
psi-continuum-v2/
    frozen original scientific implementation

psi-continuum-v2-audit/
    independent reproducibility and consistency audit
```

The frozen repository should remain unchanged so that the exact software
state associated with v2 remains inspectable.

## Interpretation

The central audit result is:

> **INTERNAL NUMERICAL REPRODUCIBILITY CONFIRMED / PUBLISHED
> MODEL-IMPLEMENTATION MISMATCH CONFIRMED**

The frozen joint minimum at (\varepsilon_0=+0.075) is
reproducible, but it is a result for the background model actually
implemented by the frozen code. It must not be presented as a likelihood
result for the documented v2 equation.

Conversely, the raw counterfactual likelihood obtained by directly
inserting the documented equation into the old likelihood pipeline is
**not** a corrected v2 inference, because the SN treatment is strongly
affected by absolute calibration and the SDSS convention also requires
correction.

A scientifically interpretable inference for the literal documented
equation would require a separately specified analysis with consistent
nuisance/calibration treatment and BAO conventions. That would
constitute new analysis rather than a modification of this frozen audit.

## Scientific policy

-   **No silent correction:** frozen v2 remains unchanged.
-   **Reproducibility first:** historical numerical results are
    independently reconstructed before interpretation.
-   **Model and implementation are tested separately:** internal code
    consistency does not establish correspondence with the documented
    equation.
-   **Audit and new research remain separate:** this repository
    documents v2; it is not a replacement cosmological analysis.

## Citation and scientific record

If v2 results are cited or discussed, the distinction between the
**documented v2 equation** and the **frozen implemented model** should
be stated explicitly.

If a formal correction to the scientific record is required, it should
be issued transparently as a validation/erratum note rather than by
rewriting the frozen repository.

## Status

**Audit completed:** Audits 01--08.

**Frozen v2:** unchanged.

**Overall result:** **INTERNAL NUMERICAL REPRODUCIBILITY CONFIRMED /
PUBLISHED MODEL-IMPLEMENTATION MISMATCH CONFIRMED.**

For detailed numerical evidence, methodology, and dataset-by-dataset
conclusions, see `AUDIT_REPORT.md`.
