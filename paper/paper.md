---
title: 'par2-circadian: temporal persistence of gene expression from second-order autoregressive eigenvalues'
tags:
  - Python
  - chronobiology
  - circadian rhythms
  - time series
  - gene expression
  - autoregressive models
authors:
  - name: Michael Whiteside
    orcid: 0009-0000-0643-5791
    affiliation: 1
affiliations:
  - name: Independent researcher, United Kingdom
    index: 1
date: 5 September 2026
bibliography: paper.bib
---

# Summary

`par2-circadian` fits a second-order autoregressive model, AR(2), to a gene
expression time series and reports the modulus of the dominant root of its
characteristic polynomial, $|\lambda|$, as a scalar measure of *temporal
persistence*: how strongly a transcript's recent past constrains its next
value. From the two fitted coefficients the package derives, per gene, the root
type (complex roots indicate an oscillatory decay mode, real roots a monotone
one), the persistence half-life in sampling intervals, the intrinsic
eigenperiod where roots are complex, and an optional residual-bootstrap
confidence interval on $|\lambda|$ obtained by resampling residuals and
refitting, so that no biological replicates are required. It provides a Python
API, a command-line interface for CSV expression matrices, and a helper that
stratifies fitted genes into user-supplied gene sets and reports their
layer-wise medians. It depends only on NumPy [@harris2020numpy], SciPy
[@virtanen2020scipy] and statsmodels [@seabold2010statsmodels], installs from
PyPI, and is tested against a bundled example dataset.

# Statement of need

Circadian transcriptomics is dominated by *rhythm detection*: JTK_CYCLE
[@hughes2010jtk], ARSER [@yang2010arser], RAIN [@thaben2014rain] and the
harmonised `MetaCycle` wrapper [@wu2016metacycle] each test a series against an
assumed periodic waveform and return a p-value for the null of arrhythmicity.
This answers whether a gene cycles, which is the right question for building
atlases [@zhang2014atlas] and is well served by existing tools
[@hughes2017guidelines; @hutchison2015limits]. It does not answer a distinct
question that arises once a dataset is in hand: *how long does this
transcript's dynamics remember?* Persistence is a continuous property, it is
defined whether or not the series is periodic, and it does not require the
period to be specified in advance.

`par2-circadian` exists to make that quantity routine to compute and to compare
across genes, tissues and conditions. An AR(2) process is the smallest linear
model that admits a complex root pair, and therefore the smallest one that can
represent damped oscillation rather than mere exponential decay; its eigenvalue
modulus is a single interpretable number bounded by the stationarity triangle.
The package targets researchers who already have a time-series matrix and want
a per-gene dynamical summary — as a covariate, a ranking, or a descriptor of how
a perturbation changes dynamics — without committing to a 24 h template.

# Comparison with existing software

The comparison is one of task, not accuracy. `dryR` [@weger2021dryr] and
`LimoRhyde` [@singer2019limorhyde] test *differential* rhythmicity between
conditions; `CosinorPy` [@moskon2020cosinorpy] and `DiscoRhythm`
[@carlucci2020discorhythm] fit cosinor models and screen for rhythmicity;
`pyBOAT` [@monke2020pyboat] performs wavelet analysis of instantaneous period
and phase. None returns a persistence statistic, and all but the wavelet
approach require a period or a period window.

The repository includes a ground-truth benchmark (`benchmarks/`) against
reference implementations of JTK_CYCLE, ARSER, Lomb–Scargle and RAIN on
simulated series with known dynamics, thresholded so that every score has a 5%
false-positive rate on white noise. Two results characterise where the package
is useful. On noise-driven AR(2) oscillators with $|\lambda| = 0.618$ at 24
timepoints, the bootstrap-supported $|\lambda|$ criterion detects 92% against
40–62% for the p-value methods. On sustained 12 h series, methods searching a
20–28 h window detect $\leq$0.5% while raw $|\lambda|$ detects 57%, because it
is period-agnostic. Conversely, on a clean sustained 24 h cosine the incumbents
reach 100% and raw $|\lambda|$ 12%: where the signal really is a fixed
sinusoid, a sinusoidal test is better.

# Limitations

Two limitations are measured in the benchmark and should govern use.
First, **autocorrelation is not oscillation**: calibrated on white noise,
$|\lambda|$ classifies plain AR(1) red noise as a positive 85% of the time, and
correctly so, since red noise is genuinely persistent. A high $|\lambda|$ is
therefore not evidence of rhythmicity, and any claim that a gene oscillates
needs a rhythm test alongside $|\lambda|$, not instead of it. Second, the
estimator is **biased upward at low persistence** — at 24 timepoints a true
$|\lambda|$ of 0.30 estimates 0.51, while at 0.90 the bias is negligible — and
bootstrap intervals under-cover in that regime. Point estimates from short
series should not be compared at two decimal places.

# Acknowledgements

The method is described in a preprint [@whiteside2026ar2]. The author declares
a patent application (GB2518973.9) relating to applications of the method;
the software itself is released under Apache-2.0.

# References
