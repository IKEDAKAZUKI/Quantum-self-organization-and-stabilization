# Equations and public functions

Equation and section references refer to the associated study, *Quantum entanglement as a resource for self-organization and stabilization*.

| Function or implementation | Article | Meaning |
|---|---|---|
| `qtp.bounds.comparison_constants` | Appendix B.2, Eqs. (38)–(42) | alpha, c0, c1, leakage steady state, Ue and the asymptotic bound |
| `qtp.bounds.formation_bound` | Eq. (28) | Formation bound for ell(0)=0 and the supplied E0 |
| `qtp.bounds.maintenance_bound` | Eq. (30) | Bound conditional on the measured initial leakage ell_f |
| `qtp.bounds.decay_convolution` | Equal-rate limits following Eq. (30) | Stable difference-of-exponentials convolution |
| `qtp.bounds.source_moments` | Eqs. (12), (35)–(36) | Constant-rate coherence and leakage expectations |
| `qtp.observables.occupation_contrast` | Eq. (3), centered storage convention | Decode the integer/leakage arrays |
| `qtp.observables.spectrum` | Eq. (18), Appendix D | Fixed-band Fourier intensity and spectral statistics |
| `qtp.geometry.measure` | Sec. 3.1 and `ANALYSIS_DETAILS.md` | Fixed morphology and spatial-extent criteria |
| `critical_jump.cpp`, field rates | Eqs. (8)–(11) | Exact continuous-time number/phase jump sampling |
| `critical_jump.cpp`, mode 1 | Appendix C | Exact phase average of nonlinear projected rates |

The maintenance bound uses the measured leakage fraction of the starting configuration, including its deviation from the stationary fraction. Both history bounds use the same stable convolution. The default d must be positive; invalid coefficient regimes raise an error rather than returning a misleading bound.

For alpha != k the convolution is evaluated as

`exp(-min(alpha,k)*t) * [-expm1(-abs(alpha-k)*t)] / abs(alpha-k)`.

For alpha == k the function returns `t*exp(-alpha*t)`. This handles equal and nearly equal rates without cancellation. Unit tests compare it with high-precision decimal evaluation and with direct integration of the contrast/leakage equations.

The source threshold check evaluates the homogeneous deterministic linear description with rapid phase relaxation. It is not an eigenanalysis of the full noisy many-body generator. The figures' growth curve uses the continuous-wavevector maximum for the reported parameters. The discrete reciprocal-grid checks are reported separately.

Geometry indicators, exact thresholds, endpoint inference and figure-selection rules are specified in [ANALYSIS_DETAILS.md](ANALYSIS_DETAILS.md). Endpoint statistics include all realizations. Figure inputs are numerical configurations, ensemble summaries and analytical curves. Image interpolation is only a rendering choice and is not part of the classifier.

## Admissible nonselective controls

Section 3.1 fixes the physical modes, finite support and generator channels. An
additional Lindblad generator M_t must have a nonpositive contribution to
mean squared contrast and zero contribution to mean leakage for every state on
that support. A discrete CPTP intervention obeys the corresponding adjoint
inequality/equality for these two observables. The separable comparison
requires the unconditional state to be mode-separable at every time, including
the continuous-projection limit. No postselection is performed.

Number projection preserves both diagonal observables. Appendix C derives its
limiting generator and its nonlinear rate average. Tests in
`tests/test_admissibility.py` check these moment identities on finite matrices
and give a first-moment-preserving map that increases the second moment. They
illustrate the need for the second-moment condition; the proof is in the article.
