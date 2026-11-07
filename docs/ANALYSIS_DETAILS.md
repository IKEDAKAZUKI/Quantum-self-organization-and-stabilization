# Numerical analysis and figure construction

This document specifies the fixed analysis of the simulated configurations. The
article defines the dynamics, the physical-mode partition, the observables and the
separable bounds. The numerical definitions below are implemented in
`qtp/observables.py`, `qtp/morphology.py`, `qtp/geometry.py`, `qtp/analysis.py`
and `qtp/verification.py`. All configuration arrays and model coefficients retain
their stored precision.

## Initial preparation and population sampling

The initial occupation field is a Gaussian random field, centered and scaled to
root-mean-square amplitude 0.02 before rounding to allowed contrast values.
The formation bound uses the actual second moment of the rounded configuration,

$$
 E_0=M^{-1}\sum_x u_x(0)^2.
$$

At initial coherence zero and no leakage, the equal mixture of
$|n,+\rangle$ and $|n,-\rangle$ is the number-diagonal mixture of
$|n,L\rangle$ and $|n,R\rangle$. The initial state is thus separable
across the physical oscillator modes. Each physical oscillator has
$2U/\varepsilon+2$ available Fock levels, which gives 202 for the common
parameters. The local invariant space has $3(2U/\varepsilon+1)$ states.

The invariant density matrix is diagonal in the joint $|n,s\rangle$ basis.
Its populations evolve with transitions $-\to+$ at rate $r+\eta$,
$+\to-$ at rate $\eta$, $\pm\to0$ at rate
$\kappa_{\mathrm l}$, and $0\to\pm$ at rate
$\kappa_{\mathrm r}/2$, together with the field jumps. Continuous-time
Gillespie sampling generates the recorded configurations. The source and
initialization hashes and `phase_draw_order` specify the random initialization
convention. See `DATA_DICTIONARY.md` and `SIMULATION.md` for the stored labels
and executable inputs.

Formation and maintenance use 48 by 56 cells, periodic boundaries and unit
lattice spacing. Formation runs to $t_f=600$, followed by quantum and
number-projected continuations of duration $\tau=120$. The dephasing data
use 24 by 28 cells and four independent realizations per parameter set.
Numerical times are in units of $D_4^{-1}$, with $D_4=1$. The bound
coefficient is evaluated at $\theta=0.15$.

## Spectral observables

For each parameter set the band is fixed by the homogeneous growth rate at
$e=1$, evaluated on the reciprocal grid of that finite array,

$$
 \mathcal B=\{\boldsymbol k\ne0:\sigma(\boldsymbol k;1)>0\},\qquad
 P_{\mathcal B}(u)=\sum_{\boldsymbol k\in\mathcal B}|A_{\boldsymbol k}(u)|^2.
$$

Let $v_x=u_x-M^{-1}\sum_yu_y$ and
$\sigma_v=(M^{-1}\sum_xv_x^2)^{1/2}$. A configuration is resolved if
$P_{\mathcal B}\ge0.005$ and
$P_{\mathcal B}/\sigma_v^2\ge0.40$. A zero-variance configuration is
unresolved. For a resolved field, the spectral nematic order is

$$
 S_2=\frac{\left|\sum_{\boldsymbol k\in\mathcal B}
 |A_{\boldsymbol k}(u)|^2 e^{2i\phi_{\boldsymbol k}}\right|}
 {P_{\mathcal B}(u)},\qquad
 \phi_{\boldsymbol k}=\arg(k_x+ik_y).
$$

## Morphology and spatial extent

Positive and negative domains satisfy $v_x>0.2\sigma_v$ and
$-v_x>0.2\sigma_v$, respectively. A signed component is a connected
component of either domain. Connectivity uses eight neighbors with periodic
identifications. Its area $A$ is the number of cells. Its perimeter $p$
counts exposed horizontal and vertical unit edges, separately from the
connectivity rule.

Components with fewer than five cells are not compact. For each remaining
component, the center is the circular mean of its coordinates in each periodic
direction. The sample covariance of shortest periodic displacements from that
center uses normalization $1/(A-1)$. With eigenvalues
$\lambda_{\max},\lambda_{\min}$, compact components satisfy all of

$$
 A\ge5,\qquad \frac{4\pi A}{\max(p^2,1)}\ge0.38,\qquad
 \sqrt{\frac{\lambda_{\max}+0.25}{\lambda_{\min}+0.25}}\le2.3.
$$

The additive 0.25 is a fixed regularization of the covariance eigenvalues. The
compact-area fraction is compact positive-component area divided by total
positive-domain area, or zero for an empty positive domain.

The morphology rules are applied in this order:

1. A resolved field with $S_2\ge0.85$ is a stripe.
2. A resolved field failing the stripe condition is a spot pattern if it has
   at least four compact positive components, compact-area fraction at least
   0.65, and largest positive-component area divided by $M$ strictly below
   0.12.
3. A resolved field with $S_2<0.85$ and compact-area fraction strictly below
   0.40 is a labyrinth.
4. Remaining resolved configurations are mixed.

These conditions were fixed with separate pilot configurations. The geometry
indicator adds the following spatial-extent requirements. A stripe must contain
a noncontractible cycle in at least one signed component. A labyrinth must
contain a winding signed component with area at least $0.15M$. Compact
positive-component centers of a spot pattern must occupy at least six tiles of
a fixed 3 by 3 partition, with its origin at the coordinate origin. Tile
occupancy counts centers, not overlap of component areas with tiles.

A winding component contains a noncontractible cycle. The coordinate lift
around a closed loop must have a nonzero system-length displacement. Merely
crossing the displayed image seam need not give a noncontractible component.
The same rules apply to quantum dynamics, projected dynamics and the
independent-preparation ensembles.

Each target geometry has indicator $\chi_g$. Its configuration intensity is
$G_g=\chi_gP_{\mathcal B}$. The ensemble includes all configurations,
with zero intensity for outcomes outside the target. Geometry classification
uses arrays before any rendering interpolation.

## Sampling units and endpoint inference

The main comparison uses sixteen independent stochastic realizations for each
parameter set, with one initial field and one disorder realization fixed.
All formation endpoints initialize both a quantum and a projected continuation.
The two evolution types have the same initial occupations and leakage, and
continuations are independent across starting configurations within each type.
The primary data consist of 48 formation and 96 continuation histories.

Formation satisfies the target criteria in 16 of 16 stripe configurations and
15 of 16 configurations for each of spot and labyrinth. Every configuration is
included in the means. Projected continuation endpoints have zero target
intensity while their residual band power can be nonzero.

For endpoint realization $i$, write $G_i=G_g(u_i)$. The formation bound
is $B_i=B_{\rm sep}(t_f)$, with that preparation's $E_0$. The maintenance
bound is $B_i=B_f(\tau;\ell_{f,i})$, with the leakage measured in the
starting configuration. Descriptive ratios divide mean intensity by mean bound.
They are distinct from the individual ratios used by the test below.

Under the separable-history hypothesis and a starting configuration $Z_i$,

$$
 \mathbb E[G_i\mid Z_i,H_{\rm sep}]\le B_i(Z_i)
 \quad\Rightarrow\quad
 \mathbb E[G_i/B_i\mid H_{\rm sep}]\le1.
$$

Set $X_i=G_i/B_i$, $Y_i=(1+X_i)/2$, and
$\mathcal E=\prod_{i=1}^{16}Y_i$. Independence gives
$\mathbb E[\mathcal E\mid H_{\rm sep}]\le1$. A zero intensity contributes
one half to the product. Markov's inequality and a Bonferroni correction across
the six prespecified formation and quantum-maintenance endpoints give adjusted
$p=\min(1,6/\mathcal E)$. The largest value is approximately
$1.9\times10^{-3}$, for spots. This tests a hypothesis about the
unconditional expectation, conditional on the simulation model and the
independence assumptions. It is not a simultaneous-in-time confidence bound.

Standard errors are calculated across the sixteen independent realizations for
primary means and across the four realizations for dephasing means. Shaded
bands and error bars display one standard error. Endpoint error bars in
Figure 3(d) divide the standard error of the mean intensity by the mean bound.

The additional ensembles contain four independent initial-field/disorder pairs,
with two stochastic realizations for each pair and parameter set. Their 24
formation and 48 continuation histories are a separate descriptive analysis.
Between-block standard errors use the four block means. They are not pooled
into the six primary endpoint tests. Spot counts change from 7 of 8 after
formation to 6 of 8 after quantum continuation. One configuration's
compact-area fraction crosses the fixed threshold while band power and spatial
coverage remain nearly unchanged, giving zero target intensity at that
continuation endpoint. Detailed results are in `data/reference/preparation_endpoints.csv`
and `data/reference/preparation_summary.csv`.

## Spatial-extent sensitivity

With morphology rules fixed, the compact-spot-center requirement was checked
at five, six and seven occupied tiles, and the winding-labyrinth area
requirement at $0.12M$, $0.15M$ and $0.18M$. The main-ensemble
endpoint classifications and pattern intensities are identical over these
choices. The defaults remain six tiles and $0.15M$.

## Representative configurations and plotted quantities

For each parameter set, Figure 1 uses the first configuration in the
predetermined realization order satisfying the target indicator at $t_f=600$.
The earlier columns use that same realization at $t=0,20,60$. Figures 3(a)
and 3(b) show the selected labyrinth before and after quantum continuation.
This selection determines only the illustration. Ensemble averages retain all
outcomes and both continuation types.

Heatmaps display occupation contrast decoded from the physical occupations and
leakage, with coordinates divided by system length. Figure 2(a) is the exact
analytical source expectation. Its transport and pattern-intensity panels use
sixteen realizations per parameter set. Figure 4(a) evaluates the homogeneous
linear description; Figure 4(b) averages band power over four dephasing cycles
per parameter set. Figure generation records its inputs under the output
folder. The reference figure PDFs are in `figures/`.

## Independent numerical validation

`python verify_results.py small --output results/small.json` compares the
Lindblad generator with the invariant population process in an independently
constructed small density-matrix calculation. A separate two-dimer calculation
applies physical-number projection to the quantum generator and compares it
with a number-and-leakage Markov generator assembled from the phase-averaged
rates. The comparisons include nonuniform physical L/R orientation
probabilities. The transition rates and propagated probabilities agree within
floating-point precision on the tested inputs.

`check` reconstructs endpoint observables and statistics and matches each
continuation's initial occupations and leakage to its formation endpoint.
`smoke` compiles the C++ generator and checks short quantum/projected
trajectories against recorded arrays. The main formation and maintenance data
remain within the occupation cutoff. Per-run event diagnostics are retained
in the corresponding `_engine.json` files. These checks and their numerical
tolerances are implemented in the verification functions and tests.
