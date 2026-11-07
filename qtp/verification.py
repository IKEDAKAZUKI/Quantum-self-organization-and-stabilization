"""Check published inputs, ensemble statistics and independent finite-system calculations.

Each action writes an optional JSON report and exits nonzero on a failed check.
Input trajectories and reference summaries are read-only.
"""
from __future__ import annotations
from qtp import __version__
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any
ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True

def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)

def digest(p: Path) -> str:
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def audit() -> dict[str, Any]:
    m = json.loads((ROOT / 'MANIFEST_SHA256.json').read_text())
    differences = []
    for name, sha in m['sha256'].items():
        p = ROOT / name
        require(p.resolve().is_relative_to(ROOT.resolve()), 'Invalid manifest path: ' + name)
        if not p.is_file():
            differences.append({'path': name, 'reason': 'missing'})
        elif digest(p) != sha:
            differences.append({'path': name, 'reason': 'hash mismatch'})
    return dict(action='audit', bundle=m['package'], files_checked=len(m['sha256']), differences=differences, errors=[] if not differences else ['Archive integrity failed'])

def prepare_science(tmp: Path) -> None:
    os.environ['NUMBA_CACHE_DIR'] = str(tmp / 'numba')
    os.environ['MPLCONFIGDIR'] = str(tmp / 'mpl')
    for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
        os.environ[key] = '1'
    (tmp / 'validation').mkdir(exist_ok=True)

def check() -> dict[str, Any]:
    import numpy as np
    import pandas as pd
    from .geometry import measure
    from .bounds import formation_bound, maintenance_bound
    primary = pd.read_csv(ROOT / 'data/reference/geometry_confirmation_timeseries.csv')
    expected_primary = primary.sort_values('t').groupby('id', sort=False).tail(1).set_index('id')
    expected_fresh = pd.read_csv(ROOT / 'data/reference/preparation_endpoints.csv').set_index('id')
    records = []
    matched = []
    err = []
    fft_err = []
    for subdir, expected, count in [('primary', expected_primary, 144), ('preparations', expected_fresh, 72)]:
        folder = ROOT / 'data' / subdir
        paths = sorted(folder.glob('*_params.json'))
        require(len(paths) == count, f'{subdir}: expected {count} histories, found {len(paths)}')
        require(len(expected) == count, f'{subdir}: reference table count mismatch')
        for ppath in paths:
            p = json.loads(ppath.read_text())
            name = p['id']
            ref = expected.loc[name]
            with np.load(folder / (name + '.npz'), allow_pickle=False) as z:
                require(len(z['t']) >= 2, name + ': fewer than two saved times')
                t = float(z['t'][-1])
                require(abs(t - p['T']) < 1e-10, name + ': endpoint time mismatch')
                require(tuple(z['n'].shape[1:]) == (p['ny'], p['nx']), name + ': array shape mismatch')
                u = p['eps'] * (z['n'][-1] - 0.5 * (z['s'][-1] == 0))
                initial_ell = float((z['s'][0] == 0).mean())
                is_maintenance = p['group'] in ('maintenance', 'fresh_maintenance')
                if is_maintenance:
                    with np.load(folder / (p['parent'] + '.npz'), allow_pickle=False) as parent:
                        n_ok = np.array_equal(z['n'][0], parent['n'][-1])
                        l_ok = np.array_equal(z['s'][0] == 0, parent['s'][-1] == 0)
                    require(n_ok and l_ok, name + ': unmatched parent number/leakage start')
                    matched.append(name)
            f = measure(u, p, p['family'])
            B = float(maintenance_bound([t], p, initial_ell)[0] if is_maintenance else formation_bound([t], p)[0])
            require(f['morphology'] == ref['morphology'], name + ': changed morphology')
            require(f['chi'] == int(ref['chi']), name + ': changed fixed gate')
            for key in ('G', 'power', 'nematic', 'compact_area_share', 'compact_tiles', 'largest_wrapping_fraction'):
                delta = abs(float(f[key]) - float(ref[key]))
                require(delta < 2e-12, f'{name}: {key} difference {delta}')
                err.append(delta)
            delta = abs(B - float(ref['B']))
            require(delta < 2e-12, name + ': envelope mismatch')
            err.append(delta)
            kx = 2 * np.pi * np.arange(p['nx']) / p['nx']
            ky = 2 * np.pi * np.arange(p['ny']) / p['ny']
            wx, wy = (2 - 2 * np.cos(kx), 2 - 2 * np.cos(ky))
            omega = wx[None, :] + wy[:, None]
            band = -p['delta'] + (p['D2'] + p['Dcl']) * omega - p['D4'] * omega ** 2 - p['zeta'] * wy[:, None] > 0
            band[0, 0] = False
            fft_power = float((abs(np.fft.fft2(u) / u.size) ** 2)[band].sum())
            fft_err.append(abs(fft_power - f['power']))
            require(fft_err[-1] < 2e-12, name + ': independent FFT disagreement')
            row = dict(id=name, family=p['family'], group=p['group'], mode=p['mode'], t=t, B=B, G=f['G'], power=f['power'], chi=f['chi'])
            if 'block' in p:
                row['block'] = p['block']
            records.append(row)
    df = pd.DataFrame(records)
    require(len(matched) == 144, 'Expected 96 primary and 48 follow-up children')
    primary_table = pd.read_csv(ROOT / 'data/reference/geometry_primary_endpoints.csv')
    groups = []
    summary_errors = []
    for (group, family, mode), q in df[df.group.isin(['confirmation', 'maintenance'])].groupby(['group', 'family', 'mode']):
        ref = primary_table[(primary_table.group == group) & (primary_table.family == family) & (primary_table['mode'] == mode)].iloc[0]
        n = len(q)
        g, B = (float(q.G.mean()), float(q.B.mean()))
        se = float(q.G.std(ddof=1) / np.sqrt(n))
        le = float(np.log((1 + q.G.to_numpy() / q.B.to_numpy()) / 2).sum())
        values = dict(mean_G=g, mean_B=B, SE_G=se, ratio=g / B, mean_power=float(q.power.mean()), log_evidence=le, familywise_markov_bound=float(min(1, 6 * np.exp(-le))))
        require(n == int(ref['n']) and q.chi.sum() == ref['hits'], 'Primary membership changed')
        for key, val in values.items():
            delta = abs(val - float(ref[key]))
            summary_errors.append(delta)
            require(delta < 2e-12, f'{group}/{family}/{mode}: summary {key}')
        groups.append(dict(group=group, family=family, mode=int(mode), n=n, accepted=int(q.chi.sum()), **values))
    fresh_table = pd.read_csv(ROOT / 'data/reference/preparation_summary.csv')
    followups = []
    for (group, family, mode), q in df[df.group.str.startswith('fresh_')].groupby(['group', 'family', 'mode']):
        ref = fresh_table[(fresh_table.group == group) & (fresh_table.family == family) & (fresh_table['mode'] == mode)].iloc[0]
        bm = q.groupby('block').G.mean()
        g, B = (float(q.G.mean()), float(q.B.mean()))
        se = float(bm.std(ddof=1) / np.sqrt(len(bm)))
        vals = dict(mean_G=g, ratio=g / B, SE_between_blocks=se)
        for key, val in vals.items():
            delta = abs(val - float(ref[key]))
            summary_errors.append(delta)
            require(delta < 2e-12, f'{group}/{family}/{mode}: block summary {key}')
        followups.append(dict(group=group, family=family, mode=int(mode), n=len(q), blocks=len(bm), accepted=int(q.chi.sum()), **vals))
    proj = df[(df['mode'] == 1) & df.group.isin(['maintenance', 'fresh_maintenance'])]
    require(len(proj) == 72 and (proj.G == 0).all(), 'Projected endpoint geometry mismatch')
    return dict(action='check', raw_endpoints_reclassified=216, primary_endpoints=144, followup_endpoints=72, matched_starts=len(matched), projected_endpoints=len(proj), all_projected_target_scores_zero=True, max_endpoint_difference=float(max(err)), max_independent_fft_difference=float(max(fft_err)), max_summary_difference=float(max(summary_errors)), primary_summaries=groups, followup_block_summaries=followups, errors=[])

def small(tmp: Path) -> dict[str, Any]:
    from .validation import density as validate
    from .validation import projected as validate_projected
    validate.ROOT = tmp
    validate_projected.ROOT = tmp
    density = validate.full_check()
    ppt = validate.witness_check()
    projected = validate_projected.full_projected()
    rates = validate_projected.local_five()
    require(density['max_density_difference'] < 1e-10, 'Full-density evolution disagreement')
    require(ppt > -1e-10, 'Weighted-PPT finite-matrix check failed')
    require(projected['max_marginal_error'] < 1e-10, 'Projected evolution disagreement')
    require(projected['lumpability_residual'] < 1e-10, 'Projected lumpability disagreement')
    require(rates['max_rate_error'] < 1e-10, 'Five-cell rate disagreement')
    return dict(action='small', full_density=density, weighted_ppt=ppt, projected_density=projected, five_cell_rates=rates, errors=[])

def smoke(tmp: Path) -> dict[str, Any]:
    """Recompile the generator and reproduce quantum and projected prefixes."""
    import numpy as np
    from .simulation import compile_engine, simulate
    from .io import load_trajectory
    engine = compile_engine(tmp / 'build')
    cases = []
    for name in ['confirm_stripe_00', 'maint_stripe_00_1']:
        report = simulate(name, tmp / name, until=10.0, engine=engine)
        p, reference = load_trajectory(name)
        with np.load(tmp / name / (name + '.npz'), allow_pickle=False) as z:
            selected = reference.t <= 10 + 1e-09
            equal = dict(times_equal=bool(np.array_equal(z['t'], reference.t[selected])), integers_equal=bool(np.array_equal(z['n'], reference.n[selected])), phase_labels_equal=bool(np.array_equal(z['s'], reference.s[selected])))
            require(all(equal.values()), name + ': short-prefix mismatch')
            cases.append(dict(id=name, compared_times=z['t'].tolist(), **equal))
    return dict(action='smoke', cases=cases, tolerance=0, generator_sha256=digest(ROOT / 'src/critical_jump.cpp'), reproducibility_scope='Exact arrays on this toolchain; cross-platform bitwise identity is not assumed.', errors=[])

def thresholds() -> dict[str, Any]:
    """Evaluate the homogeneous deterministic growth relation, Eqs. (15)--(17)."""
    import numpy as np
    from .io import load_parameters
    rows = []
    for name in ['confirm_stripe_00', 'confirm_spot_00', 'confirm_labyrinth_01']:
        _, p = load_parameters(name)
        ell = p['leak'] / (p['leak'] + p['repair'])
        ec = (2 * np.sqrt(p['delta'] * p['D4']) - p['Dcl']) / p['D2']
        e = p['r'] * (1 - ell) / (p['r'] + 2 * p['eta'] + p['leak'])
        eh = p['r'] * (1 - ell) / (p['r'] + 60 + p['leak'])
        wx = 4 * np.sin(np.pi * np.arange(p['nx']) / p['nx']) ** 2
        wy = 4 * np.sin(np.pi * np.arange(p['ny']) / p['ny']) ** 2
        w = wx[None, :] + wy[:, None]

        def growth(x):
            return float((-p['delta'] + (p['D2'] * x + p['Dcl']) * w - p['D4'] * w * w - p['zeta'] * wy[:, None]).max())
        nonzero = w > 0
        ec_grid = np.full_like(w, np.inf)
        ec_grid[nonzero] = ((p['delta'] + p['D4']*w*w + p['zeta']*wy[:, None]
                            - p['Dcl']*w)[nonzero] / (p['D2']*w[nonzero]))
        discrete_ec = float(ec_grid.min())
        require(2 * ell < discrete_ec < e, 'Discrete-grid coherence threshold ordering failed')
        require(2 * ell < ec < e, 'Coherence threshold ordering failed')
        require(growth(2 * ell) < 0 < growth(e) and growth(eh) < 0, 'Growth-sign mismatch')
        rows.append(dict(family=p['family'], separable_coherence_limit=2 * ell, coherence_threshold=float(ec), discrete_coherence_threshold=discrete_ec, quantum_stationary_coherence=e, high_dephasing_coherence=eh, separable_max_growth=growth(2 * ell), quantum_max_growth=growth(e), high_dephasing_max_growth=growth(eh)))
    return dict(action='thresholds', approximation='Homogeneous deterministic limit with rapid phase relaxation', rows=rows, errors=[])

def dephasing() -> dict[str, Any]:
    """Recompute recovery ratios from all twelve dephasing-cycle arrays."""
    import numpy as np
    from .io import load_trajectory
    from .observables import spectrum
    families = []
    for family in ['stripe', 'spot', 'labyrinth']:
        ratios = []
        for rep in range(4):
            p, z = load_trajectory(f'cycle_{family}_{rep}')
            idx = [np.flatnonzero(np.isclose(z.t, t, rtol=0, atol=1e-09)) for t in (240, 650)]
            require(all((len(j) == 1 for j in idx)), 'Missing cycle evaluation time')
            v = spectrum(z.u[[int(j[0]) for j in idx]], p)['power']
            ratios.append(float(v[1] / v[0]))
        families.append(dict(family=family, n=4, individual_ratios=ratios, mean_individual_recovery_ratio=float(np.mean(ratios))))
    return dict(action='dephasing', reference_time=240, evaluation_time=650, families=families, errors=[])

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', choices=['audit', 'check', 'small', 'smoke', 'thresholds', 'dephasing'])
    parser.add_argument('--output', type=Path, help='Explicit destination for the JSON report (otherwise stdout only).')
    args = parser.parse_args()
    try:
        if args.action == 'audit':
            report = audit()
        else:
            with tempfile.TemporaryDirectory(prefix='qtp-validation-') as td:
                tmp = Path(td)
                prepare_science(tmp)
                report = {'check': check, 'small': lambda: small(tmp), 'smoke': lambda: smoke(tmp), 'thresholds': thresholds, 'dephasing': dephasing}[args.action]()
        report.update(software_version=__version__, status='failed' if report.get('errors') else 'passed')
    except Exception as exc:
        report = dict(action=args.action, software_version=__version__, status='failed', errors=[f'{type(exc).__name__}: {exc}'])
    encoded = json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n'
    if args.output:
        from .safety import output_directory
        output_directory(args.output.parent)
        args.output.write_text(encoded)
    print(encoded)
    return 1 if report.get('errors') else 0
if __name__ == '__main__':
    raise SystemExit(main())
