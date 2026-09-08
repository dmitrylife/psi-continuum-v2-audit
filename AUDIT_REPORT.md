# Ψ-Continuum v2 --- Independent Scientific Code Audit

**Audit target:** frozen `psi-continuum-v2` / v0.2.3\
**Audit policy:** frozen v2 is not modified; all tests run externally in
`psi-continuum-v2-audit`.

## Executive status

The completed audits show that the frozen numerical pipeline is
internally reproducible in the tested components, but the implemented
non-zero-ε ΨCDM background is not the background defined by the
published v2 equation. This propagates into distances and changes the
SN, H(z), and SDSS/BOSS DR12 likelihood behaviour. A separate SDSS DR12
fiducial sound-horizon convention mismatch is also confirmed.

## Audit 01 --- ΨCDM background model

**Status: FAIL --- FROZEN IMPLEMENTATION INCONSISTENT WITH PUBLISHED
DEFINING EQUATION**

Published v2 defines \[
H\_`\Psi`{=tex}(z)=H\_`\Lambda`{=tex}(z)`\left`{=tex}(1+`\frac{\varepsilon_0}{1+z}`{=tex}`\right`{=tex}).
\]

Frozen code instead implements a normalized modified (E\^2(z)) law that
enforces (H\_`\Psi`{=tex}(0)=H_0) for every ε. The ΛCDM null limit
passes exactly, but non-zero ε does not reproduce the published
equation. At (H_0=70), ε=0.03, frozen code gives
(H\_`\Psi`{=tex}(0)=70), while the published equation gives 72.1 km s⁻¹
Mpc⁻¹.

**Conclusion:** the frozen ΨCDM background implementation is
inconsistent with the published v2 defining equation.

## Audit 02 --- Implementation trace

**Status: TRACE CONFIRMED WITH LOCAL WRAPPER**

The numerical dependency chain from the frozen ΨCDM background through
(H(z)), luminosity/angular distances and the analysis scripts was
traced. No direct implementation of the published
(`\varepsilon`{=tex}\_0/(1+z)) factor was found in the package
background path. A local DESI wrapper reconstructs (D_M=d_L/(1+z)) but
does not introduce an independent background law.

**Conclusion:** the discrepancy of Audit 01 feeds the actual frozen
analysis pipeline.

## Audit 03 --- Cosmological distances

**Status: INTERNAL PASS / PUBLISHED-SPEC FAIL**

Frozen distance identities are internally consistent. Frozen (d_L)
agrees with an independent integration of frozen (H(z)) at roughly
(2.8`\times10`{=tex}\^{-6}) relative accuracy, and the ΛCDM null limit
passes. For non-zero ε, distances derived from the published equation
differ at the percent level; for ε=0.03, (d_L) differs by about 3% at
low z and about 2% by z=3.

**Conclusion:** the distance layer correctly propagates the frozen
background, but therefore also propagates the model-definition mismatch.

## Audit 04 --- Pantheon+ SN likelihood

**Status: FROZEN REPRODUCTION PASS / PUBLISHED-SPEC STATISTICAL MISMATCH
CONFIRMED**

For 1701 Pantheon+ Hubble-flow entries with the frozen full covariance,
\[ `\chi`{=tex}\^2\_{`\Lambda`{=tex}`\mathrm{CDM}`{=tex}}=2609.19713024.
\] The frozen scan over
(-0.1`\le`{=tex}`\varepsilon`{=tex}\_0`\le0.1`{=tex}) reaches its lower
boundary at ε=-0.100 with \[
`\chi`{=tex}^2=2593.27825030,`\qquad `{=tex}`\Delta`{=tex}`\chi`{=tex}^2=-15.91888.
\] Using the literal published background in the same raw frozen
likelihood gives a radically different curve, with an unprofiled minimum
near ε=+0.055 and (`\Delta`{=tex}`\chi`{=tex}\^2`\simeq-831.60`{=tex}).
Audit 04b tests the calibration origin of this enormous raw change.

**Conclusion:** frozen SN is reproducible, but frozen and published
backgrounds do not generate the same SN likelihood curve.

## Audit 04b --- Pantheon+ SN calibration and (H_0) semantics

**Status: CALIBRATION EFFECT CONFIRMED / RESIDUAL SHAPE EFFECT DETECTED
/ BEST-FIT UNRESOLVED**

Audit 04b tested whether the very large improvement obtained when the
published v2 background equation is inserted into the frozen Pantheon+
likelihood is primarily caused by an overall supernova calibration shift
or by a redshift-dependent change in the luminosity-distance relation.

The frozen Pantheon+ data vector, full STAT+SYS covariance matrix,
distance-modulus conversion, and likelihood definition were preserved. A
single global additive distance-modulus offset was then analytically
profiled.

For the published v2 equation, the unprofiled calculation gives

\[ `\varepsilon`{=tex}\_0 = +0.055, `\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2 = -831.603. \]

After profiling over one global magnitude offset, the minimum within the
tested interval becomes

\[ `\varepsilon`{=tex}*0 = -0.100, `\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2*{`\rm prof`{=tex}} = -12.769. \]

Thus, profiling removes approximately (818.83) units of the apparent
(`\chi`{=tex}\^2) improvement, corresponding to (98.46%) of its original
magnitude.

This demonstrates that the extremely large unprofiled improvement is
dominated by an absolute SN calibration/normalization effect. In
particular, the raw published-spec minimum at
(`\varepsilon`{=tex}\_0=+0.055) implies

\[ H\_`\Psi`{=tex}(0) = H\_`\Lambda`{=tex}(0)(1+`\varepsilon`{=tex}\_0)
= 73.85 {`\rm km\,s^{-1}\,Mpc^{-1}`{=tex}} \]

for the adopted baseline
(H\_`\Lambda`{=tex}(0)=70 {`\rm km\,s^{-1}\,Mpc^{-1}`{=tex}}).

A residual redshift-dependent effect remains after removal of the global
offset. The published-spec model reaches
(`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm prof`{=tex}}`\simeq-12.77`{=tex}),
whereas the frozen implementation reaches only
(`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm prof`{=tex}}`\simeq-1.40`{=tex})
over the same (`\varepsilon`{=tex}\_0) interval. Therefore, the
discrepancy between the published equation and the frozen implementation
cannot be attributed solely to their different normalization at (z=0);
their predicted distance-redshift shapes also differ.

However, the profiled published-spec minimum occurs at the lower scan
boundary, (`\varepsilon`{=tex}\_0=-0.100). Consequently, this audit does
not determine a best-fit value of (`\varepsilon`{=tex}\_0) for the
published model.

**Conclusion:** the apparent
(`\Delta`{=tex}`\chi`{=tex}\^2`\simeq-831.6`{=tex}) improvement must not
be interpreted as evidence for the published ΨCDM dynamics. It is
predominantly an absolute-calibration effect. A smaller but non-zero
shape-dependent difference remains and must be treated separately. No
interior SN-only best-fit (`\varepsilon`{=tex}\_0) for the published
equation is established within the frozen v2 scan range.

## Audit 05 --- H(z) / Cosmic Chronometers

**Status: FROZEN REPRODUCTION PASS / PUBLISHED-SPEC STATISTICAL MISMATCH
CONFIRMED**

The frozen v2 H(z) analysis was independently reproduced using the
original 32-point H(z) compilation, the frozen diagonal likelihood, and
the original 81-point scan over (-0.2
`\leq `{=tex}`\varepsilon`{=tex}\_0 `\leq 0.2`{=tex}).

The frozen ΛCDM result is

\[ `\chi`{=tex}\^2\_{`\Lambda{\rm CDM}`{=tex}} = 11.78688559, \]

and the frozen ΨCDM implementation gives

\[
`\varepsilon`{=tex}*{0,`\rm best`{=tex}}\^{`\rm frozen`{=tex}}=+0.055,
`\qquad`{=tex} `\chi`{=tex}\^2*{`\rm best`{=tex}}=11.74585998,
`\qquad`{=tex} `\Delta`{=tex}`\chi`{=tex}\^2=-0.04103. \]

Both the ΛCDM calculation and the frozen ΨCDM scan were reproduced to
numerical precision.

The audit then repeated the same likelihood calculation with the same
data, uncertainties, (H_0), (`\Omega`{=tex}\_m), and
(`\varepsilon`{=tex}\_0) grid, changing only the ΨCDM background
relation to the published v2 defining equation,

\[ H\_`\Psi`{=tex}(z)=H\_`\Lambda`{=tex}(z)
`\left`{=tex}(1+`\frac{\varepsilon_0}{1+z}`{=tex}`\right`{=tex}). \]

Under the published equation, the minimum occurs at

\[ `\varepsilon`{=tex}*{0,`\rm best`{=tex}}\^{`\rm published`{=tex}}=0,
`\qquad`{=tex} `\chi`{=tex}\^2*{`\rm best`{=tex}}=11.78688559,
`\qquad`{=tex} `\Delta`{=tex}`\chi`{=tex}\^2=0. \]

Thus the shallow non-zero minimum at
(`\varepsilon`{=tex}\_0`\simeq`{=tex}+0.055) found by the frozen H(z)
analysis is not reproduced when the published defining equation is used.

At the frozen best-fit value (`\varepsilon`{=tex}\_0=+0.055), the
published equation instead gives

\[ `\chi`{=tex}\^2=13.69558711, `\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=+1.90870, \]

relative to ΛCDM.

Because H(z) is evaluated directly, this discrepancy does not arise from
distance integration or supernova magnitude calibration. It is a direct
statistical consequence of the different background expansion laws
implemented by the frozen code and specified by the published v2
equation.

**Conclusion:** the frozen H(z) numerical pipeline is internally
reproducible, but its non-zero (`\varepsilon`{=tex}\_0) preference
belongs to the implemented frozen background model and cannot be
attributed to the published v2 background equation.

## Audit 06A --- SDSS/BOSS DR12 BAO observable convention

**Status: CONVENTION MISMATCH CONFIRMED**

The frozen v2 SDSS DR12 data file contains the six BOSS DR12 consensus
measurements at z = 0.38, 0.51, and 0.61.

The corresponding published BOSS observables are

\[ D_M(z)`\frac{r_{s,\rm fid}}{r_s(z_d)}`{=tex} \]

and

\[ H(z)`\frac{r_s(z_d)}{r_{s,\rm fid}}`{=tex}, \]

with

\[ r\_{s,`\rm fid`{=tex}}=147.78 {`\rm Mpc`{=tex}}. \]

The frozen v2 loader and diagnostic code instead describe these
quantities as (D_M/r_s) and (H(z)r_s), which is not the convention
represented by the numerical data vector.

More importantly, the frozen joint likelihood constructs its SDSS model
vector directly from

\[ D_M(z),`\qquad `{=tex}H(z), \]

without applying the BOSS fiducial sound-horizon factors

\[ `\frac{r_{s,\rm fid}}{r_s}`{=tex}
`\quad`{=tex}`\text{and}`{=tex}`\quad`{=tex}
`\frac{r_s}{r_{s,\rm fid}}`{=tex}. \]

Therefore the frozen SDSS DR12 likelihood is strictly equivalent to the
published BOSS likelihood only in the special case

\[ r_s=r\_{s,`\rm fid`{=tex}}. \]

The frozen parameter classes instead adopt (r_d=147.0) Mpc, whereas the
BOSS DR12 fiducial value associated with this data vector is (147.78)
Mpc.

**Conclusion:** the SDSS DR12 data vector is correctly recognizable as
the standard BOSS consensus vector, but its observable convention is
misdocumented in the frozen v2 code and the required fiducial
sound-horizon scaling is not explicitly implemented. With the frozen
choice (r_d=147.0) Mpc, this produces a small but genuine convention
mismatch in the BAO model vector.

## Audit 06B --- SDSS/BOSS DR12 BAO statistical impact

**Status: FROZEN REPRODUCTION PASS / BAO CONVENTION MISMATCH CONFIRMED /
PUBLISHED-SPEC STATISTICAL MISMATCH CONFIRMED**

The frozen ΛCDM BAO likelihood is reproduced: \[
`\chi`{=tex}\^2=11.140641849565, \] agreeing with the historical frozen
reference to (4.35`\times10`{=tex}\^{-10}).

Applying the BOSS fiducial sound-horizon convention while retaining the
frozen background gives \[ `\chi`{=tex}\^2=7.951527983437,
`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm convention`{=tex}}=-3.189113866128.
\]

For the ε scan, frozen raw and convention-corrected frozen backgrounds
both decrease toward the positive scan boundary: \[
`\varepsilon`{=tex}*{0,`\rm best`{=tex}}\^{`\rm frozen`{=tex}}=+0.100,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-0.93081, \] and \[
`\varepsilon`{=tex}*{0,`\rm best`{=tex}}\^{`\rm frozen,corr`{=tex}}=+0.100,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-0.57272. \]

With the literal published-v2 background and the same corrected BOSS
convention, the behaviour changes qualitatively: \[
`\varepsilon`{=tex}*{0,`\rm best`{=tex}}\^{`\rm published`{=tex}}`\simeq-0.023`{=tex},`\quad`{=tex}
`\chi`{=tex}\^2*{`\rm best`{=tex}}=3.01799833,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-4.93353. \]

At ε=+0.05, the convention correction changes χ² by about -3.00, while
changing the frozen background to the published law changes it by about
+41.74.

**Conclusion:** the frozen SDSS likelihood is reproducible, but it
contains a genuine fiducial-scaling mismatch for the adopted (r_d=147.0)
Mpc. More importantly, the published background produces qualitatively
different ε dependence, so the frozen SDSS non-zero-ε behaviour cannot
be attributed to the published v2 equation.

## Audit 07 --- DESI DR2 BAO

**Status: FROZEN REPRODUCTION PASS / PUBLISHED-SPEC STATISTICAL MISMATCH
CONFIRMED**

The DESI DR2 Gaussian BAO likelihood was independently reconstructed
using the frozen 13-element data vector and its full
(13`\times13`{=tex}) covariance matrix. The covariance is symmetric and
positive definite. The frozen ΛCDM value is reproduced to numerical
precision:

\[ `\chi`{=tex}\^2\_{`\Lambda{\rm CDM}`{=tex}}=40.298419969637. \]

At (`\varepsilon`{=tex}\_0=0), the frozen ΨCDM branch coincides exactly
with ΛCDM. Across the tested interval
(-0.1`\leq`{=tex}`\varepsilon`{=tex}\_0`\leq0.1`{=tex}), however, the
two non-zero-(`\varepsilon`{=tex}\_0) likelihood curves are
qualitatively different. The frozen implementation decreases toward the
positive scan boundary,

\[
`\varepsilon`{=tex}\_{0,`\rm best`{=tex}}\^{`\rm frozen`{=tex}}=+0.100,`\qquad`{=tex}
`\chi`{=tex}\^2=20.62428459,`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-19.67414, \]

so no interior frozen DESI-only best fit is established in this
interval.

Using the literal published-v2 equation instead gives an interior
minimum near

\[
`\varepsilon`{=tex}\_{0,`\rm best`{=tex}}\^{`\rm published`{=tex}}`\simeq-0.022`{=tex},`\qquad`{=tex}
`\chi`{=tex}\^2=12.32747374,`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-27.97095. \]

At (`\varepsilon`{=tex}\_0=+0.05), the frozen implementation gives
(`\chi`{=tex}\^2=27.30368388) and
(`\Delta`{=tex}`\chi`{=tex}\^2=-12.99474), whereas the published
equation gives approximately (`\chi`{=tex}\^2=285.16974) and
(`\Delta`{=tex}`\chi`{=tex}\^2=+244.87132).

The legacy DESI diagnostic that returned identical ΛCDM and ΨCDM values
used (`\varepsilon`{=tex}\_0=0); its equality is therefore an expected
null-limit result and is not evidence of a DESI likelihood failure.

**Conclusion:** the DESI DR2 data vector, covariance treatment, and
ratio observables are internally sound in the tested frozen pipeline.
The positive-(`\varepsilon`{=tex}\_0) DESI preference of the frozen
implementation belongs to the normalized frozen (E\^2) model and is not
a result for the published v2 defining equation.

## Audit 08 --- Joint likelihood reconstruction

**Status: FROZEN REPRODUCTION PASS / PUBLISHED-SPEC JOINT COMPARISON
MISMATCH**

The complete joint likelihood was reconstructed on a common 201-point
grid, (-0.1`\leq`{=tex}`\varepsilon`{=tex}\_0`\leq0.1`{=tex}), using
Pantheon+ SN, 32 H(z) measurements, SDSS/BOSS DR12 BAO, and DESI DR2
BAO.

Two branches were evaluated:

1.  **Frozen branch:** the exact frozen-v2 background and likelihood
    conventions, including the historical raw SDSS vector.
2.  **Published-spec comparison branch:** the literal published-v2
    background equation with the corrected SDSS/BOSS fiducial-(r_d)
    convention. The original raw SN likelihood was retained so that the
    effect of changing the background law could be isolated;
    consequently this branch is a diagnostic comparison and not a
    corrected cosmological fit.

### Frozen reconstruction

At (`\varepsilon`{=tex}\_0=0), the independently reconstructed
components are

\[ `\chi`{=tex}\^2\_{`\rm SN`{=tex}}=2609.197130241862,`\quad`{=tex}
`\chi`{=tex}\^2\_{H(z)}=11.786885594164,`\quad`{=tex}
`\chi`{=tex}\^2\_{`\rm SDSS`{=tex}}=11.140641849565,`\quad`{=tex}
`\chi`{=tex}\^2\_{`\rm DESI`{=tex}}=40.298419969637, \]

giving

\[ `\chi`{=tex}\^2\_{`\rm total`{=tex}}=2672.423077655228. \]

This differs from the historical frozen reference (2672.423077660000) by
only (-4.77`\times10`{=tex}\^{-9}). The historical
(`\varepsilon`{=tex}\_0=+0.05) point is also reproduced:

\[ `\chi`{=tex}\^2\_{`\rm total`{=tex}}(0.05)=2666.255500468532, \]

against the reference (2666.255500470000).

The complete frozen scan has its minimum at

\[ `\boxed{\varepsilon_{0,\rm best}^{\rm frozen}=+0.075}`{=tex} \]

with

\[ `\chi`{=tex}\^2\_{`\rm min`{=tex}}=2665.510865815413,`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-6.912211839815. \]

The component changes relative to (`\varepsilon`{=tex}\_0=0) at this
minimum are

\[
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm SN`{=tex}}=+10.93035,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{H(z)}=-0.03468,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm SDSS`{=tex}}=-0.72923,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm DESI`{=tex}}=-17.07866. \]

Their sum agrees with the joint (`\Delta`{=tex}`\chi`{=tex}\^2) to
approximately (7`\times10`{=tex}\^{-13}). Thus the frozen joint
improvement is driven overwhelmingly by DESI DR2, partly offset by the
SN penalty.

### Literal published-spec comparison

At (`\varepsilon`{=tex}\_0=0), the published-spec comparison branch
gives

\[ `\chi`{=tex}\^2\_{`\rm total`{=tex}}=2669.238605650585. \]

It is not expected to equal the frozen ΛCDM total because this
comparison branch deliberately applies the corrected SDSS/BOSS
fiducial-(r_d) convention. The SDSS convention change alone shifts the
(`\varepsilon`{=tex}\_0=0) contribution from (11.14064185) to
approximately (7.95152728).

The raw published-spec joint comparison reaches a minimum at

\[ `\varepsilon`{=tex}\_0`\simeq`{=tex}+0.041,`\qquad`{=tex}
`\chi`{=tex}\^2`\simeq2100.05910543`{=tex},`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2`\simeq-569.17950`{=tex}. \]

This number must **not** be interpreted as evidence for the published
ΨCDM equation. At this point the component changes are approximately

\[
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm SN`{=tex}}=-785.163,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{H(z)}=+1.036,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm SDSS`{=tex}}=+31.091,`\quad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2\_{`\rm DESI`{=tex}}=+183.856. \]

The enormous negative total is therefore produced by the raw SN
normalization/calibration response and is opposed by both BAO datasets.
Audit 04b independently showed that (98.46%) of the very large raw SN
improvement is removed by profiling one global magnitude offset. The
published-spec joint branch is consequently a model-definition
diagnostic, not a statistically valid corrected-v2 parameter inference.

The frozen and published-spec branches have different best-fit
(`\varepsilon`{=tex}\_0) values and different
(`\Delta`{=tex}`\chi`{=tex}\^2(`\varepsilon`{=tex}\_0)) curves. In both
branches, the joint total equals the sum of the four independently
reconstructed components to machine precision.

**Conclusion:** the complete frozen joint pipeline is numerically
reproducible. Its minimum at (`\varepsilon`{=tex}\_0=+0.075) is a result
of the frozen implemented background model and frozen likelihood
conventions. It cannot be attributed to the published v2 defining
equation.

## Final audit matrix

  -----------------------------------------------------------------------
  Audit                   Component               Status
  ----------------------- ----------------------- -----------------------
  01                      ΨCDM background         **FAIL** --- frozen
                          definition              implementation differs
                                                  from published equation

  02                      Implementation trace    **TRACE CONFIRMED** ---
                                                  discrepancy reaches the
                                                  numerical pipeline

  03                      Distances               **INTERNAL PASS /
                                                  PUBLISHED-SPEC FAIL**

  04                      Pantheon+ SN            **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  04b                     SN calibration          **CALIBRATION EFFECT
                          diagnostic              CONFIRMED / BEST-FIT
                                                  UNRESOLVED**

  05                      H(z)                    **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  06A                     SDSS/BOSS observable    **CONVENTION MISMATCH
                          convention              CONFIRMED**

  06B                     SDSS/BOSS statistical   **FROZEN PASS /
                          impact                  CONVENTION MISMATCH /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  07                      DESI DR2 BAO            **FROZEN PASS /
                                                  PUBLISHED-SPEC
                                                  MISMATCH**

  08                      Joint likelihood        **FROZEN PASS /
                                                  PUBLISHED-SPEC JOINT
                                                  MISMATCH**
  -----------------------------------------------------------------------

## Final scientific conclusion

The audit establishes two distinct facts that must be kept separate.

First, the **frozen v2 software is substantially internally
reproducible**. The tested ΛCDM values, individual frozen likelihood
calculations, historical control points, and complete frozen joint scan
can be reconstructed to numerical precision. In particular, the
independent joint reconstruction recovers the frozen minimum

\[ `\varepsilon`{=tex}\_0=+0.075,`\qquad`{=tex}
`\Delta`{=tex}`\chi`{=tex}\^2=-6.91221. \]

Therefore the principal problem is not an inability to reproduce the
frozen software output.

Second, the \*\*non-zero-(`\varepsilon`{=tex}\_0) model implemented by
the frozen software is not the model defined by the published v2
equation\*\*,

\[
H\_`\Psi`{=tex}(z)=H\_`\Lambda`{=tex}(z)`\left`{=tex}(1+`\frac{\varepsilon_0}{1+z}`{=tex}`\right`{=tex}).
\]

The frozen code instead uses a normalized modified-(E\^2) prescription
that enforces (H\_`\Psi`{=tex}(0)=H_0). The two formulations share the
exact (`\varepsilon`{=tex}\_0=0) ΛCDM limit but diverge for non-zero
(`\varepsilon`{=tex}\_0). The discrepancy propagates into distances and
materially changes the SN, H(z), SDSS/BOSS, DESI DR2, and joint
likelihood curves.

A second, independent issue is present in the frozen SDSS/BOSS DR12
treatment: the standard BOSS fiducial sound-horizon scaling is not
explicitly applied even though the frozen parameter class uses
(r_d=147.0) Mpc rather than the (147.78) Mpc fiducial value associated
with the data vector. This is a smaller effect than the background-law
discrepancy but is numerically non-negligible.

The Pantheon+ analysis also requires a separate statistical
qualification. When the literal published equation is inserted into the
raw frozen SN likelihood, most of the very large apparent
(`\chi`{=tex}\^2) improvement is an absolute calibration/normalization
effect. Profiling one global magnitude offset removes approximately
(98.46%) of that raw improvement. A residual shape-dependent effect
remains, but its minimum lies at the tested scan boundary; Audit 04b
therefore does not establish a published-model SN-only best fit.

Accordingly:

-   numerical preferences obtained from the frozen
    non-zero-(`\varepsilon`{=tex}\_0) pipeline are valid descriptions of
    the **implemented frozen model**, not of the published v2 defining
    equation;
-   the frozen joint minimum at (`\varepsilon`{=tex}\_0=+0.075) must not
    be reported as a likelihood result for the published equation;
-   the raw published-spec joint minimum near
    (`\varepsilon`{=tex}\_0=+0.041) and
    (`\Delta`{=tex}`\chi`{=tex}\^2`\simeq-569`{=tex}) must likewise not
    be reported as evidence for the published equation, because it is
    dominated by the known SN calibration response;
-   a scientifically interpretable fit of the literal published equation
    requires a separately specified inference analysis with explicit
    nuisance/calibration treatment and consistent BAO conventions. Such
    an analysis would be new work and is outside this frozen-v2 audit.

## Audit disposition

**Frozen v2 repository:** remains unchanged and archived.

**Audit repository:** contains the independent validation scripts and
numerical outputs for Audits 01--08.

**Overall audit result:** **INTERNAL NUMERICAL REPRODUCIBILITY CONFIRMED
/ PUBLISHED MODEL-IMPLEMENTATION MISMATCH CONFIRMED.**

The appropriate scientific response is transparent documentation of the
discrepancy and, if the published v2 result is part of a public
scientific record, a validation/erratum note that clearly distinguishes
the published equation from the frozen implemented model. No silent
modification of the frozen v2 repository is warranted.
