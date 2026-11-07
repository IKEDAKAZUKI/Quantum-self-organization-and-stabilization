#!/usr/bin/env python3
"""Check the completion-time normalization and the plotted publication data."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from reference_time import completed_reference_time
ROOT=Path(__file__).resolve().parent
reference, definition=completed_reference_time(ROOT)
summary=pd.read_csv(ROOT/'analysis/summary.csv')
series=pd.read_csv(ROOT/'analysis/timeseries.csv')
values=pd.read_csv(ROOT/'figures/figure5_values.csv')
info=json.loads((ROOT/'figures/figure5_provenance.json').read_text())
assert definition['baseline_runs']==48
assert reference==info['reference_time']['value']
assert max(r['canvas_overflow_points'] for r in info['panel_layouts'])<=.5
max_error=0.
for family in ('stripe','spot','labyrinth'):
    runs=series[(series.tag=='main')&(series.protocol=='autonomous')&(series.family==family)]
    assert runs.id.nunique()==16
    for run,g in runs.groupby('id'):
        base=g[g.t==reference]
        assert len(base)==1 and base.activated_fraction.iloc[0]==1
        expected=g.power.to_numpy()/float(base.power.iloc[0])
        err=float(np.max(abs(expected-g.power_over_t_star.to_numpy())))
        assert err<1e-12
        max_error=max(max_error,err)
    direct=runs.groupby('t').power_over_t_star.agg(['mean','sem']).reset_index()
    plotted=values[(values.panel=='c')&(values.family==family)]
    joined=plotted.merge(direct,on='t',suffixes=('_plot','_direct'))
    assert np.allclose(joined.mean_plot,joined.mean_direct,rtol=0,atol=1e-12)
    assert np.allclose(joined.standard_error,joined['sem'],rtol=0,atol=1e-12)
    at_reference=plotted[plotted.t==reference]
    assert np.isclose(at_reference['mean'].iloc[0],1.)
    assert np.isclose(at_reference.standard_error.iloc[0],0.)
    lo=plotted['mean']-plotted.standard_error
    hi=plotted['mean']+plotted.standard_error
    assert lo.min()>info['panel_c_ylim'][0] and hi.max()<info['panel_c_ylim'][1]
assert all(x['h_mean']==1 for x in info['maps'] if x['t']>=reference)
report={'status':'passed','baseline_runs':48,'reference_time':reference,
        'latest_conversion_event_time':definition['latest_conversion_event_time'],
        'maximum_ratio_error':max_error,'all_reference_normalizations_equal_one':True,
        'all_error_bands_inside_axes':True,'all_postreference_snapshots_converted':True,
        'all_panels_inside_canvas':True}
(ROOT/'validation/publication_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
