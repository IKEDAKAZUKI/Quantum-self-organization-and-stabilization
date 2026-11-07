# Figure input summaries

These tables are recomputed from the recorded trajectories by `analyze_trajectories.py`.

- `primary_summary.csv` contains formation and maintenance means, standard errors, separable bounds and transport diagnostics for each family and evolution mode.
- `dephasing_summary.csv` contains the mean fixed-band power and its standard error across the four dephasing cycles in each family.

To render the figures directly from these summaries, run from the repository root:

```bash
python make_figures.py --analysis data/summary --output results/figures
```

The complete analysis command regenerates both tables under `results/analysis`.
See [Analysis details](../../docs/ANALYSIS_DETAILS.md) for the statistical definitions.
