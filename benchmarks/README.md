# Ground-truth benchmark: AR(2) |λ| against established rhythm-detection methods

This suite answers two questions that reanalysis of real datasets cannot, because real
data has no known answer:

1. **Does |λ| recover the persistence it claims to measure**, and how precise is it at
   the sampling densities circadian experiments actually use?
2. **What does |λ| detect that JTK_CYCLE, ARSER, Lomb–Scargle and RAIN miss, and where
   do they beat it?**

Comparators are the reference implementations, not reimplementations: JTK_CYCLE, ARSER
and Lomb–Scargle via `MetaCycle::meta2d`, and RAIN from Bioconductor. `dryR` and
`LimoRhyde` are deliberately excluded — both test *differential* rhythmicity between
conditions, which is a different task from single-series detection.

## Reproducing

```bash
python3 simulate.py            # 2,000 series x {12,18,24,48} timepoints, 6 scenario families
Rscript run_incumbents.R       # JTK / ARS / LS via MetaCycle, plus RAIN
python3 analyse.py             # AR(2) fits, bootstrap CIs, AUC, |λ| recovery
python3 calibrated.py          # power at a matched 5% false-positive rate
python3 figures.py
```

R dependencies: `MetaCycle` (CRAN), `rain` (Bioconductor). Sampling is 2 h; the
period search window given to every method is 20–28 h.

## Scenarios

| scenario | process | truth |
|---|---|---|
| `ar2_<λ>` | AR(2) with complex roots of modulus λ ∈ {0.30, 0.50, 0.618, 0.75, 0.90}, 24 h eigenperiod, noise-driven | rhythmic, known \|λ\| |
| `cosinor` | sustained 24 h cosine + white noise | rhythmic |
| `damped` | 24 h cosine with exponentially decaying amplitude | rhythmic |
| `ultradian` | sustained 12 h cosine | rhythmic, non-24 h |
| `white` | white noise | arrhythmic (calibration null) |
| `red` | AR(1), φ = 0.8 | arrhythmic but autocorrelated |

Because `|λ|` is a model fit rather than a p-value, methods are compared at a **matched
null**: each score is thresholded so its false-positive rate on white noise is exactly
5%, and the call rate in every other scenario is then directly comparable.
`PAR2_support` is the fraction of residual-bootstrap refits that yield complex roots
with an eigenperiod inside the 20–28 h window; `PAR2_lambda` is the raw point estimate.

## What the benchmark shows

**1. |λ| is the most powerful method on noise-driven oscillators — the regime real
transcript dynamics occupy.** At n = 24 and true |λ| = 0.618, `PAR2_support` detects 92%
against 56% (JTK), 62% (LS), 51% (ARSER) and 40% (RAIN). The advantage is largest for
weakly persistent oscillators (48% vs 25–34% at |λ| = 0.30), which is where circadian
screens lose most genes.

**2. |λ| is period-agnostic, and the others are not.** On 12 h ultradian series the
p-value methods searching 20–28 h detect ≤0.5% (RAIN 34% at n = 48), while the raw |λ|
detects 57% at n = 24 and 97% at n = 48. This is the clearest justification for the metric: it measures
persistence without being told the period.

**3. The incumbents win on deterministic waveforms.** For a clean sustained cosine, JTK,
LS, ARSER and RAIN all reach 100% at n = 24, `PAR2_support` 83%, raw |λ| 12%. For
amplitude-decaying cosines the gap is worse (0.98 vs 0.56). A method built for a fixed
sinusoid beats an AR(2) fit whenever the signal really is a fixed sinusoid.

**4. The important limitation: autocorrelation is not oscillation.** Calibrated to 5%
FPR on white noise, `PAR2_support` calls plain AR(1) red noise "circadian" 85% of the
time at n = 24, against 38% (JTK) and 29% (RAIN). |λ| measures persistence and does so
correctly here — red noise *is* persistent — but a high |λ| on its own is not evidence
of rhythmicity, and the eigenperiod window does not fix that. Any claim that a gene is
an oscillator needs a rhythm test alongside |λ|, not instead of it.

**5. |λ| is biased upward, and the bias concentrates exactly where the biology is
interpreted.** At n = 24 a true |λ| of 0.30 estimates 0.51 (bias +0.21, RMSE 0.25);
at 0.90 the bias is −0.02. Bootstrap 95% intervals under-cover at low persistence
(77% at |λ| = 0.30) and are well calibrated above ~0.6. Mean interval width at n = 24
is 0.46 at |λ| = 0.618, so per-gene point estimates cannot be compared at two decimal
places.

Consequence for hierarchy claims: because low-|λ| series are inflated more than
high-|λ| series, an observed clock-minus-target gap is *compressed* relative to the
truth — the direction of the bias does not manufacture a gap, it shrinks one. That is
reassuring for the clock > target result and fatal for any narrow per-gene claim about
a specific value such as 1/φ ≈ 0.618, which the CI width already rules out.

Full tables in `results/` (`detection.csv`, `calibrated_power.csv`, `recovery.csv`,
`sensitivity_by_scenario.csv`); figures in `figures/`.
