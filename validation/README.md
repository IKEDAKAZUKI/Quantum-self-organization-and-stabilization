# Numerical and software validation

The reports in this directory were generated from the distributed code and recorded data. All listed checks passed in the environment recorded in [environment.json](environment.json).

| Check | Scope | Report |
| --- | --- | --- |
| Unit tests | 71 tests covering bounds, configuration handling, initial conditions, generators and figure labels | [unit_tests.json](unit_tests.json) |
| Full analysis | 3,888 primary configurations, 72 independent-preparation endpoints and 3,912 dephasing configurations | [analysis.json](analysis.json) |
| Endpoints and matched starts | 216 endpoints, 144 matched continuation starts and 72 projected endpoints | [check.json](check.json) |
| Finite-system generators | Independent 81-dimensional density-matrix calculations and 200 five-cell projected-rate cases | [small.json](small.json) |
| Short trajectory reproduction | Exact occupation, phase-label and time arrays for quantum and projected prefixes through time 10 | [smoke.json](smoke.json) |
| Instability thresholds | Continuous-wavevector and finite-grid thresholds for the three parameter families | [thresholds.json](thresholds.json) |
| Dephasing recovery | Recovery ratios recomputed from all 12 recorded cycles | [dephasing.json](dephasing.json) |

The full-analysis results agree with the reference observables to floating-point precision. The independent density-matrix comparisons have maximum discrepancies below $10^{-14}$. Short-trajectory equality applies to the recorded toolchain; cross-platform bitwise equality is not assumed.

The recorded trajectories are used directly for the ensemble analysis. Validation does not require regenerating all 228 simulations.

## Run the checks

From the repository root:

```bash
python verify_results.py audit
python -m unittest discover -s tests -v
python verify_results.py check
python verify_results.py small
python verify_results.py smoke
python verify_results.py thresholds
python verify_results.py dephasing
```

Add `--output results/<name>.json` to a verification command to save its report. Each failed check returns a nonzero exit status. The `audit` command checks the distributed files against `MANIFEST_SHA256.json`. Generated files under `results/` are not part of that manifest.

Figure input hashes and panel-boundary checks are recorded in [figure_provenance.json](../figures/figure_provenance.json).
