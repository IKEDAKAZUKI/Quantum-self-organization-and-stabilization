#!/usr/bin/env python3
"""Independent field-moment and number-projected-rate checks."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(285194)
projection_error=0.;drift_error=0.;variance_error=0.
for _ in range(8192):
    u=rng.uniform(-1,1,5);h=rng.integers(0,2,5);leak=rng.integers(0,2,5)
    edge=rng.uniform(1.2,1.9,4);f0=rng.normal();eff=.91;eps=.01;smooth=.002;Db=1e-5
    w=edge*(u[0]-u[1:])/2
    f=f0+np.sum(w*eff*(h[0]+h[1:]))
    roots=[]
    for mask in range(32):
        signs=np.array([1 if (mask>>k)&1 else -1 for k in range(5)])*(1-h)*(1-leak)
        z=f0+np.sum(w*(signs[0]+signs[1:]+eff*(h[0]+h[1:])))
        roots.append(np.hypot(z,smooth))
    original=float(np.mean(roots))
    aa=[float(w[k]) for k in range(4) if h[k+1]==0 and leak[k+1]==0 and w[k]!=0]
    if h[0]==0 and leak[0]==0 and w.sum()!=0:aa.append(float(w.sum()))
    fast=np.mean([np.sqrt((f+sum(a*(1 if (mask>>j)&1 else -1) for j,a in enumerate(aa)))**2+smooth**2) for mask in range(1<<len(aa))])
    projection_error=max(projection_error,abs(original-fast))
    qp=(fast+f)/(2*eps)+Db/eps**2;qm=(fast-f)/(2*eps)+Db/eps**2
    assert qp>=0 and qm>=0
    drift_error=max(drift_error,abs(eps*(qp-qm)-f))
    variance_error=max(variance_error,abs(eps**2*(qp+qm)-(eps*original+2*Db)))
report=dict(cases=8192,max_projected_mean_root_error=projection_error,max_drift_error=drift_error,max_variance_error=variance_error)
assert projection_error<1e-12 and drift_error<1e-12 and variance_error<1e-12
(ROOT/'validation'/'rate_checks.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
