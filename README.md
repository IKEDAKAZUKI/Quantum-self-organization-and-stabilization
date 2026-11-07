# Quantum self-organization and stabilization

**Local quantum entanglement as a resource for macroscopic spatial order.**

[![Stripe, spot and labyrinth patterns in a bosonic lattice](figures/research_overview.png)](figures/README.md)

*Representative occupation-contrast configurations at $t=600$, shown on the same color scale. [Explore the figures](figures/README.md).*

Simulation code, numerical data and analysis tools for **Quantum entanglement as a resource for self-organization and stabilization**, by **Kazuki Ikeda**, University of Massachusetts Boston.

## Research overview

Can entanglement be necessary for a spatial pattern to form and persist? This study addresses that question in a noisy bosonic lattice. Entanglement between the two modes within each cell supports transport that amplifies spatial variations, producing stripes, spots and labyrinths across the array.

The calculations connect three aspects of this process:

- **Formation.** Quantum dynamics produces mean pattern intensities above rigorous bounds for separable evolution under the specified physical constraints.
- **Stabilization.** Quantum and continuously number-projected dynamics start from the same occupation and leakage configurations. Quantum evolution maintains the pattern intensity, while the projected evolution loses the target geometry.
- **Control.** Increasing local dephasing suppresses pattern amplification. Restoring the original dephasing rate permits recovery.

The simulation samples the exact invariant population jump process of the Lindblad model. The entanglement criteria refer to the unconditional physical-mode state, not to a single sampled configuration. The separable bounds include arbitrary classical correlations and the control conditions described in [Methods](docs/METHODS.md).

## Getting started

Clone the repository and run the commands from its root directory. The tested environment is Python **3.13.5** with the dependencies pinned in `requirements-tested.txt`. A C++17 compiler is needed only for new simulations and generator tests.

```bash
git clone https://github.com/IKEDAKAZUKI/Quantum-self-organization-and-stabilization.git
cd Quantum-self-organization-and-stabilization
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-tested.txt
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in PowerShell. The ready-to-view [figures](figures/README.md) and [reference tables](data/reference/) do not require running a simulation.

## Reproduce the analysis and figures

```bash
python analyze_trajectories.py --output results/analysis
python make_figures.py --analysis results/analysis --output results/figures
```

These commands read the recorded trajectories, recompute the observables and ensemble statistics, and produce the four study figures plus the overview image in PDF and PNG formats. Plotted-value tables and input hashes are written alongside the figures. All realizations enter the ensemble averages, including configurations outside the target geometry. Generated outputs are kept under `results/`, separate from the distributed data.

The complete analysis uses the saved trajectories. **No new stochastic simulation is required to reproduce the figures.** To render only the figures from the included summary tables:

```bash
python make_figures.py --analysis data/summary --output results/figures
```

## Data

All **228 recorded trajectories**, their parameters, initial configurations and event diagnostics are included.

| Dataset | Contents | Trajectories |
| --- | --- | ---: |
| [`data/primary/`](data/primary/) | 48 formation runs and their matched quantum and number-projected continuations | 144 |
| [`data/preparations/`](data/preparations/) | Independent initial-field and disorder blocks, with matched continuations | 72 |
| [`data/dephasing/`](data/dephasing/) | Four suppression-and-recovery cycles for each pattern family | 12 |

The primary formation and maintenance ensembles use $48\times56$ cells and 16 realizations per family. Dephasing cycles use $24\times28$ cells. The independent-preparation analysis uses four preparation/disorder blocks as its statistical units.

The [data dictionary](docs/DATA_DICTIONARY.md) explains array axes, state labels, run identifiers and metadata. In particular, stored `n` is a centered integer label, and occupation contrast is `u = eps * (n - 0.5 * (s == 0))`. In projected mode, `s=1` denotes a nonleaked cell, not nonzero coherence. The CSV field `energy` denotes mean squared occupation contrast, not Hamiltonian energy.

## Validation

```bash
python verify_results.py audit --output results/integrity.json
python -m unittest discover -s tests -v
python verify_results.py check --output results/endpoints.json
python verify_results.py small --output results/generators.json
python verify_results.py smoke --output results/short_trajectories.json
```

These checks cover file integrity, parameter handling, all 216 formation/maintenance endpoints, 144 matched continuation starts, independent finite-system generator calculations, and short quantum/projected simulations. A failed check returns a nonzero exit status. The `thresholds` and `dephasing` actions provide additional checks. See [validation results](validation/README.md) for the recorded reports.

## Run a simulation

```bash
python run_simulation.py confirm_stripe_00 --until 10 --output results/quantum_run
python run_simulation.py maint_stripe_00_1 --until 10 --output results/projected_run
```

Omit `--until` to run the full duration. Explicit CFG and JSON inputs are supported, and existing outputs are not overwritten. Set `CXX` to select a C++17 compiler. See [Simulation guide](docs/SIMULATION.md) for custom parameters, direct C++ execution and binary conversion. Bitwise agreement of new stochastic trajectories depends on the compiler and random-number implementation.

## Repository guide

| Location | Purpose |
| --- | --- |
| [`qtp/`](qtp/) | Observables, separable bounds, geometry classification, analysis and plotting |
| [`src/critical_jump.cpp`](src/critical_jump.cpp) | Continuous-time quantum and number-projected jump dynamics |
| [`figures/`](figures/) | Study figures, overview image and captions |
| [`docs/`](docs/) | Methods, data formats, analysis definitions and simulation instructions |
| [`tests/`](tests/) | Automated numerical and software tests |
| [`provenance/`](provenance/) | Generator and initialization hashes |

## Citation and permissions

Please cite **Kazuki Ikeda, “Quantum entanglement as a resource for self-organization and stabilization”**, and identify the repository version or commit used. Machine-readable software citation metadata is provided in [CITATION.cff](CITATION.cff).

Reuse and redistribution permissions are described in [LICENSE.md](LICENSE.md). Questions and reproducible bug reports can be submitted through the repository's issue tracker.
