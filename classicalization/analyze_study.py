#!/usr/bin/env python3
"""Analyze every completed run without rejecting morphology outcomes."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from qtp.observables import spectrum
from qtp.geometry import measure
from reference_time import completed_reference_time
ROOT=Path(__file__).resolve().parent


def correlation(a,b):
    a=a.ravel()-a.mean();b=b.ravel()-b.mean();d=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/d) if d else 0.

def main():
    t_star, reference = completed_reference_time(ROOT)
    (ROOT / "analysis").mkdir(exist_ok=True)
    (ROOT / "analysis" / "reference_time.json").write_text(json.dumps(reference, indent=2) + "\n")
    series=[];ends=[];morph=[]
    for file in sorted((ROOT/'trajectories').glob('*.npz')):
        name=file.stem;p=json.loads(file.with_name(name+'_params.json').read_text())
        z=np.load(file);u=p['eps']*(z['n']-.5*(z['s']==0));t=z['t']
        log=pd.read_csv(file.with_suffix('.csv'))
        assert np.allclose(t,log.t)
        spec=spectrum(u,p);log['power']=spec['power']
        first=json.loads(file.with_name(name+'_engine.json').read_text())
        for k in ['family','rep','protocol','tag']:log[k]=p[k]
        log['id']=name
        isq=p['quantum_enabled']==1
        log['negativity_lower']=np.maximum(log.coherence-2*log.leakage,0)/4
        log['negativity_upper']=log.inactive_nonleaked_fraction/2 if isq else 0.
        j40=int(np.argmin(abs(t-40)));jstar=int(np.argmin(abs(t-t_star)))
        log['power_over_t40']=log.power/log.power.iloc[j40]
        log['power_over_t_star']=log.power/log.power.iloc[jstar]
        log['correlation_t40']=[correlation(a,u[j40]) for a in u]
        log['correlation_t_star']=[correlation(a,u[jstar]) for a in u]
        # A trajectory is classical in its environment-record sector after the
        # final activation. The unconditional density retains unsampled tails.
        tc=first['all_activated_time']
        ja=int(np.searchsorted(t,tc)) if tc>=0 else -1
        if ja>=0:
            log['power_over_absorption']=log.power/log.power.iloc[ja]
        else:log['power_over_absorption']=np.nan
        series.append(log)
        # Geometry is checked at these recorded times, not at every SSA event.
        chosen=[0,20,40,60,100,150,200,300,t_star,600,900,1200]
        indices={int(np.argmin(abs(t-v))) for v in chosen if v<=t[-1]}
        if ja>=0:indices.add(ja)
        mcache={}
        for j in sorted(indices):
            m=measure(u[j],p,p['family']);mcache[j]=m
            morph.append(dict(id=name,family=p['family'],rep=p['rep'],protocol=p['protocol'],tag=p['tag'],t=float(t[j]),activated=float(z['h'][j].mean()),morphology=m['morphology'],chi=m['chi'],G=m['G'],power=m['power'],any_wrap=m['any_wrap'],largest_wrapping_fraction=m['largest_wrapping_fraction'],compact_tiles=m['compact_tiles']))
        row=dict(id=name,family=p['family'],rep=p['rep'],protocol=p['protocol'],tag=p['tag'],T=float(t[-1]),first_activation=first['first_activation_time'],absorption_time=tc,activated_final=float(z['h'][-1].mean()),power_final=float(log.power.iloc[-1]),power_t40=float(log.power.iloc[j40]),power_t_star=float(log.power.iloc[jstar]),power_ratio_t40=float(log.power_over_t40.iloc[-1]),power_ratio_t_star=float(log.power_over_t_star.iloc[-1]),correlation_t40=float(log.correlation_t40.iloc[-1]),correlation_t_star=float(log.correlation_t_star.iloc[-1]),coherence_t_point2=float(log.coherence.iloc[np.argmin(abs(t-.2))]),morphology_final=mcache[len(t)-1]['morphology'],chi_final=mcache[len(t)-1]['chi'],G_final=mcache[len(t)-1]['G'],chi_t_star=mcache[jstar]['chi'],G_t_star=mcache[jstar]['G'],boundary_hits=first['boundary_hits'],field_events=first['field_events'],phase_events=first['phase_events'],max_stencil_error=first['max_stencil_error'])
        if ja>=0:
            row.update(power_at_absorption=float(log.power.iloc[ja]),power_ratio_absorption=float(log.power_over_absorption.iloc[-1]),min_power_ratio_after_absorption=float(log.power_over_absorption.iloc[ja:].min()),correlation_absorption=correlation(u[-1],u[ja]),morphology_at_absorption=mcache[ja]['morphology'],chi_at_absorption=mcache[ja]['chi'],chi_all_checked_after_absorption=int(all(x['chi'] for j,x in mcache.items() if j>=ja)))
        ends.append(row)
    S=pd.concat(series,ignore_index=True);E=pd.DataFrame(ends);G=pd.DataFrame(morph)
    for name,data in [('timeseries',S),('endpoints',E),('morphology',G)]:data.to_csv(ROOT/'analysis'/f'{name}.csv',index=False)
    keys=['tag','protocol','family','t']
    cols=['power','coherence','leakage','activated_fraction','negativity_lower','negativity_upper','quantum_transport','classical_transport','power_over_t40','power_over_t_star','correlation_t40','correlation_t_star']
    ag=S.groupby(keys)[cols].agg(['mean','sem','count']);ag.columns=['_'.join(x) for x in ag.columns];ag.reset_index().to_csv(ROOT/'analysis'/'summary.csv',index=False)
    # Apply the witness after ensemble averaging, retaining its physical meaning.
    summary=pd.read_csv(ROOT/'analysis'/'summary.csv')
    summary['ensemble_negativity_lower']=np.maximum(summary.coherence_mean-2*summary.leakage_mean,0)/4
    summary.to_csv(ROOT/'analysis'/'summary.csv',index=False)
    out=[]
    for key,g in E.groupby(['tag','protocol','family']):
        d=dict(zip(['tag','protocol','family'],key));d['runs']=len(g)
        for c in ['absorption_time','power_final','power_ratio_t_star','correlation_t_star','power_ratio_absorption','min_power_ratio_after_absorption','coherence_t_point2']:
            if c in g:d[c+'_mean']=g[c].mean();d[c+'_sem']=g[c].sem()
        d['morphologies']=g.morphology_final.value_counts().to_dict();d['target_hits']=int(g.chi_final.sum());d['fully_activated_runs']=int((g.activated_final==1).sum());d['boundary_hits']=int(g.boundary_hits.sum());out.append(d)
    def finite_json(value):
        if isinstance(value, dict): return {k: finite_json(v) for k,v in value.items()}
        if isinstance(value, list): return [finite_json(v) for v in value]
        if isinstance(value, (float, np.floating)) and not np.isfinite(value): return None
        return value
    (ROOT/'analysis'/'results.json').write_text(json.dumps(finite_json(out),indent=2,allow_nan=False)+'\n')
    main_rows=[]
    for family in ('stripe', 'spot', 'labyrinth'):
        g=E[(E.tag=='main')&(E.protocol=='autonomous')&(E.family==family)]
        main_rows.append(dict(family=family,realizations=len(g),reference_time=t_star,
            conversion_time_mean=g.absorption_time.mean(),conversion_time_sem=g.absorption_time.sem(),
            power_retention_percent=100*g.power_ratio_t_star.mean(),
            power_retention_percent_sem=100*g.power_ratio_t_star.sem(),
            configuration_correlation=g.correlation_t_star.mean(),
            configuration_correlation_sem=g.correlation_t_star.sem(),
            final_target_morphologies=int(g.chi_final.sum())))
    pd.DataFrame(main_rows).to_csv(ROOT/'analysis'/'main_results.csv',index=False)
    print(E.groupby(['tag','protocol','family'])[['power_final','activated_final','chi_final','power_ratio_t_star','correlation_t_star']].mean().to_string())
if __name__=='__main__':main()
