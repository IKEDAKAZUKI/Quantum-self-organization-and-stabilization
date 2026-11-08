# Equations and public functions

Equation and section references follow the manuscript with Eqs. (1)–(56), *Quantum entanglement as a resource for self-organization and stabilization*, with the environmental model in Section 5 and numerical methods in Appendix D.

| Function or implementation | Article | Meaning |
|---|---|---|
| `qtp.bounds.comparison_constants` | Appendix B.2, Eqs. (43)–(47) | alpha, c0, c1, leakage steady state, Ue and the asymptotic bound |
| `qtp.bounds.formation_bound` | Eq. (28) | Formation bound for ell(0)=0 and the supplied E0 |
| `qtp.bounds.maintenance_bound` | Eq. (30) | Bound conditional on the measured initial leakage ell_f |
| `qtp.bounds.decay_convolution` | Equal-rate limits following Eq. (30) | Stable difference-of-exponentials convolution |
| `qtp.bounds.source_moments` | Eq. (12), Appendix A, Eqs. (40)–(41) | Constant-rate coherence and leakage expectations |
| `qtp.observables.occupation_contrast` | Eq. (3), centered storage convention | Decode the integer/leakage arrays |
| `qtp.observables.spectrum` | Eq. (18), Appendix D | Fixed-band Fourier intensity and spectral statistics |
| `qtp.geometry.measure` | Sec. 3.1 and `ANALYSIS_DETAILS.md` | Fixed morphology and spatial-extent criteria |
| `critical_jump.cpp`, field rates | Eqs. (8)–(11) | Exact continuous-time number/phase jump sampling |
| `critical_jump.cpp`, mode 1 | Eq. (32), Appendix C, Eqs. (53)–(55) | Exact phase average of nonlinear projected rates |

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

## Environmental extension

The [environmental model description](../classicalization/docs/MODEL.md) uses the following equations from Section 5 and Appendix D.

| Quantity or implementation | Article | Meaning |
|---|---|---|
| `Q_x` in the enlarged drift | Eq. (33) | Coherence-dependent transport before conversion and classical transport afterward |
| `A_x`, `Lambda_x` and environmental jumps | Eqs. (34)–(36) | Local activity, environmental response rate and number-state recording |
| Bounds on mean local negativity | Eq. (37) | Lower and upper bounds for the unconditional two-mode reduced states |
| Environmental parameters | Eq. (56) | Response rates, activity scale and classical transport gain |

The original formation and maintenance bounds apply to the fixed channels of Eq. (10). Their stated scope does not include the additional transport activated by the environmental extension.
