#!/usr/bin/env python3
"""Validate trajectories, source provenance, and irreversible classical sectors."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
sources={hashlib.sha256(p.read_bytes()).hexdigest():str(p.relative_to(ROOT)) for p in (ROOT/'src').glob('*.cpp')}
records=[]
for f in sorted((ROOT/'trajectories').glob('*.npz')):
    p=json.loads(f.with_name(f.stem+'_params.json').read_text())
    e=json.loads(f.with_name(f.stem+'_engine.json').read_text())
    log=pd.read_csv(f.with_suffix('.csv'));z=np.load(f)
    n,s,h,t=z['n'],z['s'],z['h'],z['t']
    assert n.shape==s.shape==h.shape==(len(t),p['ny'],p['nx'])
    assert np.all(np.diff(t)>0) and t[0]==0 and t[-1]==p['T']
    assert np.all(np.isin(s,[-1,0,1])) and np.all(np.isin(h,[0,1]))
    assert np.all(np.diff(h.astype(int),axis=0)>=0)
    assert np.max(abs(n))<=round(p['U']/p['eps'])
    assert not np.any((h==1)&(s==-1))
    assert p['source_sha256'] in sources
    init=ROOT/'configs'/p['init'];assert hashlib.sha256(init.read_bytes()).hexdigest()==p['initial_sha256']
    coherence=np.mean(np.where(h==0,s,0),axis=(1,2)) if p['quantum_enabled'] else np.zeros(len(t))
    assert np.allclose(coherence,log.coherence,atol=1e-14,rtol=0)
    assert np.allclose((s==0).mean((1,2)),log.leakage,atol=1e-14,rtol=0)
    assert np.allclose(h.mean((1,2)),log.activated_fraction,atol=1e-14,rtol=0)
    events=pd.read_csv(f.with_name(f.stem+'_events.csv'))
    assert len(events)==e['activation_events']
    if len(events):assert events.site.nunique()==len(events) and events.t.min()>=0 and events.t.max()<=p['T']
    if e['all_activated_time']>=0:
        j=np.searchsorted(t,e['all_activated_time']);assert np.all(h[j:]==1) and np.all(coherence[j:]==0)
    assert e['max_stencil_error']<1e-10
    records.append(dict(id=f.stem,source=sources[p['source_sha256']],events=e['field_events']+e['phase_events']+e['leak_events']+e['repair_events']+e['activation_events'],boundary_hits=e['boundary_hits'],all_activated_time=e['all_activated_time'],status='passed'))
result=dict(status='passed',completed_runs=len(records),total_events=sum(x['events'] for x in records),boundary_hits=sum(x['boundary_hits'] for x in records),runs=records)
(ROOT/'validation'/'data_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:v for k,v in result.items() if k!='runs'})
