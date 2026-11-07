# Data dictionary

## Ensembles and independence

`data/primary/` contains 48 formation trajectories (16 for each of stripe, spot and labyrinth), each followed by one quantum and one number-projected continuation. There are 144 primary histories in total. The 16 stochastic realizations share a fixed number preparation and a fixed disorder realization within each family. Every endpoint enters its ensemble, including a zero target intensity for an unsuccessful geometry.

`data/preparations/` contains four independent initial-field/disorder blocks. Each block has two stochastic realizations per family: 24 formed parents and 48 continuations, 72 histories total. Between-block standard errors use the four block means, not eight independent preparations. These ensembles are separate from the six primary endpoint tests.

`data/dephasing/` contains 12 cycles, four per family, used in Figure 4. These cycles are separate from the 216 formation/maintenance histories. Their initial occupation arrays and recorded parameters are included.

## Raw files

For each run prefix `ID`:

* `ID.npz`: arrays `t`, `n`, `s`, scalar `eps`.
* `ID_params.json`: full-precision physical parameters and run metadata.
* `ID.cfg`: generator configuration. `init` is resolved relative to this file.
* `ID_init.txt`: flattened centered integer and initial label, two columns, x varying fastest.
* `ID.csv`: event-log observables and cumulative counts for primary/preparation runs.
* `ID_engine.json`: cumulative event counts, `boundary_hits`, and `stencil_error` for primary/preparation runs.

The dephasing dataset includes trajectories, parameters and the original initialization files. Short-prefix reruns verify these initializations. New executions also generate the log and engine diagnostics.

NPZ shapes are `t: (nt,)`, `n,s: (nt, ny, nx)`. `n` is int32, `s` is int8, `t` is float64. Data are ordered by increasing saved time.

| Ensemble | nx × ny | Saved times | Statistical role |
|---|---|---|---|
| Primary formation | 48 × 56 | 0,10,...,600 | 16 independent dynamical seeds per family, fixed preparation/disorder |
| Primary maintenance | 48 × 56 | 0,5,...,120 | Quantum/projected branches paired by their common parent |
| Independent-preparation formation | 48 × 56 | 0,10,...,600 | Two seeds in each of four preparation blocks |
| Independent-preparation maintenance | 48 × 56 | 0,5,...,120 | Paired branches from every formed parent |
| Dephasing cycle | 24 × 28 | 0,2,...,650 | Four dynamical realizations per family |

Formation summaries in the article use the stored 20-unit subset. Maintenance uses all saved times. Cycle curves use all saved times. Interpolation is not added to the stored trajectories.

## Integer labels and physical occupations

Let `j` denote the saved value `n`, and let `N=j+U/eps` be the physical number-sector index in the article. For nonleaked quantum cells, `s=+1` or `-1` represents

`(|N,N+1> + s |N+1,N>)/sqrt(2)`.

For `s=0`, the cell is in the leakage state `|N,N>`. The field is

`u = eps * (j - 0.5 * (s == 0))`.

In `mode=0`, the unconditional local coherence is obtained by averaging the signed label over the appropriate ensemble. Its spatial average within a realization is `s.mean()`.

In `mode=1`, the program samples the closed number-and-leakage marginal of the continuously number-projected dynamics. `s=1` means nonleaked, and `s=0` means leaked. It does not label `|N,+>`. Physical coherence is exactly zero in this control. `qtp.observables.mean_coherence(s, mode)` implements the distinction.

The value `s=9` is allowed only in initialization text files: it instructs the generator to sample a phase label using `m0`. Saved NPZ labels are never 9. The `phase_draw_order` configuration flag specifies the uniform-to-label map: 0 uses plus-first (primary and preparations), 1 uses minus-first (cycles). Both have probability `(1+m0)/2` for +1. The difference matters for deterministic seed reproduction, not transition rates. `phaseinit=1` overrides phase preparation to +1; none of the reported main ensembles uses that override.

## Parameters and metadata

`mode`: 0 for quantum, 1 for exact number projection. `family` is the intended parameter family, not a successful classifier result. `id` is the trajectory identifier and matches the output filename prefix and execution report. For a copied input with a different identifier, the runner records that identifier in `source_id`. `group` distinguishes `confirmation`, `maintenance`, `fresh_formation`, and `fresh_maintenance`. `parent` points to the corresponding formation ID in the same directory. Continuation time starts at zero even though the parent ends at t=600.

`seed` initializes the dynamical random-number stream, including sentinel-9 initialization. `disseed` initializes static onsite and coupling disorder. `initseed`, when present, records the separate initial-number preparation seed. Exact initialization is supplied as a file; the public runner does not infer it from a seed. `rep` distinguishes stochastic realizations. `block` identifies an independently prepared number-field/disorder pair, and is the independent unit in the descriptive preparation tests.

The physics keys match Table 1: `Db=D_bg`, `smooth=s_0`, `leak=kappa_l`, `repair=kappa_r`, `disorder=w`, with all remaining coefficient names as written. `eps` and `U` define the finite occupation range. `T` and `save` specify terminal time and save interval. `E0` is the actual initial mean squared contrast. For cycles, `t1=250`, `t2=320`, `eta_mid=30` define the dephasing window.

`engine_sha256` and `initial_sha256` identify the source and initialization used to generate each recorded trajectory. New runs record the hashes of the inputs actually used. `provenance/generator.json` records the generator hashes associated with the distributed datasets.

## Event CSV and engine diagnostics

`time` is the saved time. `energy` is `mean(u**2)`, the mean squared contrast, not a Hamiltonian energy. `mean` is `mean(u)`. `coherence` is `mean(s)` for mode 0 and zero for mode 1. `leak_fraction` is `mean(s==0)`. `transport_power` is the weighted contrast-growth contribution T in Eq. (22). `weighted_leak` is the spatial mean of P0,x f_x. The local transport witness subtracts `2*J*weighted_leak + 2*J*eps**2*leak_fraction`, where `J=D2*(1+disorder)`.

`field_events`, `phase_events`, `leak_events`, and `repair_events` are cumulative counts. `max_abs` is the maximum absolute contrast. `boundary_hits` counts events that reach an occupation endpoint, and `stencil_error` is the maximum discrepancy between incrementally updated and directly recomputed linear stencil values. These are retained numerical diagnostics.

## Derived summaries

`G` is a single configuration's target pattern intensity and `chi` is its fixed target indicator. `power` is band power. `B` is the relevant separable bound. `mean_G`, `SE_G`, `mean_B` and `ratio=mean_G/mean_B` are descriptive summaries. Failure contributes G=0 and remains in n. `hits` counts target geometries. The endpoint inference uses individual `G_i/B_i` and the fixed factors `(1+G_i/B_i)/2`; it does not substitute the descriptive ratio. `SE_between_blocks` refers to the four independent preparation-block means.

`data/reference/` contains saved summaries for regression checks. Fresh analyses and figures are written under the selected output directory, not over these reference files.

## Input snapshots for new simulations

New output parameter files record `configuration_source_format` (`cfg` or
`json`) and `configuration_source_sha256` for the actual selected input file.
`engine_sha256` identifies the compiled source, and `initial_sha256` identifies
the chosen initialization. The output `init` path is relative to the output
parameter file; `E0` is recomputed from that initialization and `eps`.

These fields describe one executed configuration. Editing a sibling file does
not alter the recorded source or any stored trajectory.

Both CFG and JSON inputs are validated against the same parameter schema before defaults are applied. Unknown or duplicate names are rejected. An omitted documented optional parameter uses its engine default; a misspelled name never selects a default silently.

`data/reference/preparation_endpoints.csv` and `preparation_summary.csv` contain the endpoint and block-level reference statistics for the independent-preparation ensembles. Event diagnostics are recorded once per trajectory in `ID_engine.json`.
