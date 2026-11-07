"""Determine a common recorded time after all baseline environmental conversions."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np


def completed_reference_time(root: Path) -> tuple[float, dict]:
    """Return the first common saved time at which all baseline records have h=1.

    Completion is assessed for the recorded finite ensemble, not asserted for
    every branch of the unconditional density operator at a fixed finite time.
    """
    rows = []
    common = None
    for path in sorted((root / 'trajectories').glob('main_autonomous_*.npz')):
        metadata = json.loads(path.with_name(path.stem + '_engine.json').read_text())
        with np.load(path, allow_pickle=False) as z:
            times = z['t']
            complete = np.flatnonzero(np.all(z['h'] == 1, axis=(1, 2)))
            if len(complete) == 0:
                raise ValueError(f'{path.stem}: environmental conversion is incomplete.')
            first = float(times[complete[0]])
            saved = set(float(x) for x in times)
        exact = float(metadata['all_activated_time'])
        if not 0 <= exact <= first:
            raise ValueError(f'{path.stem}: inconsistent conversion event time.')
        common = saved if common is None else common & saved
        rows.append({'id': path.stem, 'conversion_event_time': exact,
                     'first_complete_saved_time': first})
    if not rows or not common:
        raise ValueError('No baseline autonomous ensemble or common observation grid.')
    latest = max(row['first_complete_saved_time'] for row in rows)
    candidates = sorted(t for t in common if t >= latest)
    if not candidates:
        raise ValueError('No common saved time after the final environmental conversion.')
    reference = candidates[0]
    return reference, {
        'symbol': 't_*', 'value': reference,
        'definition': 'First common saved time at which every baseline autonomous environmental record has h_x=1 at every cell.',
        'baseline_runs': len(rows),
        'latest_conversion_event_time': max(row['conversion_event_time'] for row in rows),
        'runs': rows,
    }
