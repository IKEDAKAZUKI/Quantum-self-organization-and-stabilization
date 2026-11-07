# Running the continuous-time generator

The root-level `run_simulation.py` command combines compilation, configuration, execution and conversion.

## Input selection

- An explicit `.cfg` is parsed directly. Its settings take precedence over any
  differently valued sibling JSON because the sibling is not read.
- An explicit JSON path selects exactly that file. Recorded run IDs,
  extensionless prefixes and NPZ paths select `ID_params.json`.
- `init` is resolved relative to the selected file. CFG requires it. JSON can
  omit it only to use the conventional colocated `ID_init.txt`.
- Optional engine settings take the defaults in `qtp.io.ENGINE_DEFAULTS`, which
  correspond to `src/critical_jump.cpp`. There is no CFG/JSON merge.
- CFG parsing accepts blank lines and full-line `#` comments. Numeric values
  must be finite. Seeds are exact unsigned 64-bit integers. Duplicate and
  unknown parameter keys are rejected in both CFG and JSON before defaults are applied. Inline comments after numeric values are not
  part of the syntax. Path values, including spaces, extend to the line end.

Use a copy of a CFG together with its initialization when changing parameters.
Analyze the newly generated `ID_params.json`/NPZ, rather than pairing an edited
CFG with an old trajectory. Recorded-trajectory loading rejects CFG input to
prevent that accidental relabelling of scientific data.

## Commands
All distributed configurations use a colocated `ID_init.txt` and can also be run directly. Paths inside a configuration are resolved relative to that configuration, independently of the shell's current directory.

From the project root:

```bash
mkdir -p results/build
c++ -O3 -std=c++17 src/critical_jump.cpp -o results/build/critical_jump
results/build/critical_jump data/primary/confirm_stripe_00.cfg results/direct/confirm_stripe_00
python convert_trajectory.py results/direct/confirm_stripe_00.bin --eps 0.01 --output results/direct/confirm_stripe_00.npz
```

The direct example runs the complete recorded duration. For a shorter check use

```bash
python run_simulation.py confirm_stripe_00 --until 10 --output results/short
```

The output directory must not already contain files with the same run ID. The wrapper copies the initialization and writes parameters matching the executed prefix. It leaves the packaged input untouched. Generated artifacts may be placed under `results/` or in an external directory, not inside the packaged data or source directories.

The binary header is three little-endian int32 values: nx, ny, M=nx*ny. Each subsequent record is float64 time, M int32 centered occupations, and M int8 state labels. Records have no padding. Conversion rejects truncated files and invalid labels. NPZ contains t, n, s and eps; CSV and JSON remain alongside it.

The C++ executable checks the same dynamical keys and integer fields as the Python CFG reader. It rejects unknown or duplicate keys, nonintegral integer fields, seeds outside the unsigned 64-bit range, nonfinite numbers, invalid rates, unsupported modes, and inconsistent dephasing schedules before creating trajectory files. Integer-valued decimal and scientific forms such as `phase_draw_order 1.0` and `seed 1e3` are accepted without rounding. Initialization paths are reported in errors for missing, truncated, oversized or invalid initial configurations. Dynamics use the Gillespie algorithm and a Fenwick tree for cumulative event-rate selection. The `phase_draw_order` setting specifies the initialization convention used by each dataset. Data-generator hashes and this convention are documented in `provenance/`.

## Output identifiers

The selected input filename determines the output prefix. The saved parameters
and execution report use that prefix as `id`. If a copied input contains a
different `id`, it is saved as `source_id`. The numerical parent relationship
remains in `parent` and is independent of the output identifier.
