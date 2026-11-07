"""Recompute ensemble observables from the distributed integer trajectories."""
from __future__ import annotations
from pathlib import Path
from typing import Iterable
import json
import hashlib
import numpy as np
import pandas as pd
from .io import ROOT, load_trajectory
from .bounds import formation_bound, maintenance_bound
from .geometry import measure
from .observables import spectrum, mean_coherence
from .safety import output_directory


def _summary(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows = []
    for key, q in frame.groupby(columns, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        n = len(q)
        g = float(q.G.mean())
        B = float(q.B.mean())
        loge = float(np.log((1 + q.G.to_numpy() / q.B.to_numpy()) / 2).sum())
        row = dict(zip(columns, key))
        row.update(n=n, hits=int(q.chi.sum()), mean_G=g, SE_G=float(q.G.std(ddof=1) / np.sqrt(n)), mean_B=B, ratio=g / B, mean_power=float(q.power.mean()), log_evidence=loge, familywise_markov_bound=float(min(1, 6 * np.exp(-loge))))
        for name in ['coherence', 'leak_fraction', 'transport_power', 'transport_witness']:
            if name in q:
                row[name] = float(q[name].mean())
                row['SE_' + name] = float(q[name].std(ddof=1) / np.sqrt(n))
        rows.append(row)
    return pd.DataFrame(rows)


def analyze_primary(output: Path) -> dict:
    """Reclassify every saved main-ensemble configuration used in Figures 2--3.

    Formation is evaluated on the reported 20-unit grid, maintenance on every
    saved time. Failed or non-target geometries remain zero in each average.
    """
    output = output_directory(output)
    rows = []
    cache = output / 'configuration_cache'
    cache.mkdir(exist_ok=True)
    code_hash = hashlib.sha256(b''.join(((ROOT / 'qtp' / name).read_bytes() for name in ['analysis.py', 'geometry.py', 'morphology.py', 'bounds.py', 'observables.py']))).hexdigest()
    for file in sorted((ROOT / 'data/primary').glob('*_params.json')):
        name = file.name[:-12]
        stamp = hashlib.sha256(file.read_bytes() + file.with_name(name + '.npz').read_bytes() + file.with_name(name + '.csv').read_bytes() + code_hash.encode()).hexdigest()
        cache_csv = cache / (name + '.csv')
        cache_json = cache / (name + '.json')
        if cache_csv.is_file() and cache_json.is_file() and (json.loads(cache_json.read_text()).get('sha256') == stamp):
            rows.extend(pd.read_csv(cache_csv).to_dict('records'))
            continue
        before = len(rows)
        p, z = load_trajectory(file)
        logs = pd.read_csv(file.with_name(p['id'] + '.csv'))
        if len(logs) != len(z.t) or not np.allclose(logs.time, z.t, rtol=0, atol=1e-10):
            raise ValueError(f"Log/snapshot time mismatch: {p['id']}")
        u = z.u
        formation = p['group'] == 'confirmation'
        leakage = float((z.s[0] == 0).mean())
        B = formation_bound(z.t, p) if formation else maintenance_bound(z.t, p, leakage)
        coh = mean_coherence(z.s, int(p['mode']))
        if not np.allclose(coh, logs.coherence, rtol=0, atol=2e-12):
            raise ValueError(f"Coherence semantics mismatch: {p['id']}")
        J = p['D2'] * (1 + p['disorder'])
        for j, t in enumerate(z.t):
            if formation and (not np.isclose(t % 20, 0, rtol=0, atol=1e-09)):
                continue
            f = measure(u[j], p, p['family'])
            lr = logs.iloc[j]
            rows.append(dict(id=p['id'], group=p['group'], family=p['family'], mode=p['mode'], rep=p['rep'], t=float(t), B=float(B[j]), **f, coherence=float(coh[j]), leak_fraction=float(lr.leak_fraction), transport_power=float(lr.transport_power), transport_witness=float(lr.transport_power - 2 * J * lr.weighted_leak - 2 * J * p['eps'] ** 2 * lr.leak_fraction)))
        pd.DataFrame(rows[before:]).to_csv(cache_csv, index=False)
        cache_json.write_text(json.dumps({'sha256': stamp, 'id': name}) + '\n')
    raw = pd.DataFrame(rows)
    summary = _summary(raw, ['group', 'family', 'mode', 't'])
    end = summary[(summary.group == 'confirmation') & (summary.t == 600) | (summary.group == 'maintenance') & (summary.t == 120)]
    raw.to_csv(output / 'primary_timeseries.csv', index=False)
    summary.to_csv(output / 'primary_summary.csv', index=False)
    end.to_csv(output / 'primary_endpoints.csv', index=False)
    ref = pd.read_csv(ROOT / 'data/reference/geometry_confirmation_timeseries.csv')
    joined = raw.merge(ref, on=['id', 't'], suffixes=('', '_reference'), validate='one_to_one')
    if len(joined) != len(ref) or len(raw) != len(ref):
        raise ValueError('Primary time-series membership changed')
    errors = {key: float(np.max(np.abs(joined[key] - joined[key + '_reference']))) for key in ['G', 'B', 'power', 'nematic', 'compact_area_share', 'transport_witness']}
    if max(errors.values()) > 2e-12:
        raise ValueError(f'Time-series regression failed: {errors}')
    return dict(dataset='primary', histories=144, configurations=len(raw), max_differences=errors, status='passed')


def analyze_preparations(output: Path) -> dict:
    """Recompute the 72 preparation/disorder endpoints; use four blocks for SE."""
    output = output_directory(output)
    rows = []
    for file in sorted((ROOT / 'data/preparations').glob('*_params.json')):
        p, z = load_trajectory(file)
        t = float(z.t[-1])
        f = measure(z.u[-1], p, p['family'])
        B = maintenance_bound(t, p, float((z.s[0] == 0).mean())) if p['group'] == 'fresh_maintenance' else formation_bound(t, p)
        rows.append(dict(id=p['id'], family=p['family'], group=p['group'], mode=p['mode'], block=p['block'], rep=p['rep'], t=t, B=float(B), **f))
    raw = pd.DataFrame(rows)
    summary = []
    for key, q in raw.groupby(['group', 'family', 'mode']):
        blocks = q.groupby('block').G.mean()
        summary.append(dict(zip(['group', 'family', 'mode'], key), n=len(q), blocks=len(blocks), hits=int(q.chi.sum()), mean_G=float(q.G.mean()), mean_B=float(q.B.mean()), ratio=float(q.G.mean() / q.B.mean()), SE_between_blocks=float(blocks.std(ddof=1) / np.sqrt(len(blocks)))))
    raw.to_csv(output / 'preparation_endpoints.csv', index=False)
    pd.DataFrame(summary).to_csv(output / 'preparation_summary.csv', index=False)
    return dict(dataset='preparations', histories=len(raw), independent_blocks=4, status='passed')


def analyze_dephasing(output: Path) -> dict:
    """Recompute cycle band-power curves and individual pre/post recovery ratios."""
    output = output_directory(output)
    rows = []
    ratios = []
    for fam in ['stripe', 'spot', 'labyrinth']:
        for rep in range(4):
            name = f'cycle_{fam}_{rep}'
            p, z = load_trajectory(name)
            s = spectrum(z.u, p)
            for j, t in enumerate(z.t):
                rows.append(dict(id=name, family=fam, rep=rep, t=float(t), power=float(s['power'][j])))
            i0 = np.flatnonzero(z.t == 240)
            i1 = np.flatnonzero(z.t == 650)
            if len(i0) != 1 or len(i1) != 1:
                raise ValueError('Missing cycle normalization times')
            ratios.append(dict(id=name, family=fam, rep=rep, recovery=float(s['power'][i1[0]] / s['power'][i0[0]])))
    raw = pd.DataFrame(rows)
    summary = raw.groupby(['family', 't']).power.agg(n='size', mean='mean', SE='sem').reset_index()
    raw.to_csv(output / 'dephasing_timeseries.csv', index=False)
    summary.to_csv(output / 'dephasing_summary.csv', index=False)
    pd.DataFrame(ratios).to_csv(output / 'dephasing_recovery.csv', index=False)
    return dict(dataset='dephasing', histories=12, configurations=len(raw), status='passed')


def analyze_all(output: Path) -> dict:
    output = output_directory(output)
    results = [func(output) for func in (analyze_primary, analyze_preparations, analyze_dephasing)]
    report = dict(status='passed', analyses=results)
    (output / 'analysis_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
