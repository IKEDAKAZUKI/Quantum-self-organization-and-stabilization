#!/usr/bin/env python3
"""Parameter checks and quantum-only reference for the autonomous model."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import json
from run_study import execute,FAMILIES,ROOT
cases=[]
for f in FAMILIES:
    cases.append((f,0,'no_environment',1200,48,56,'reference',{}))
    for rep in range(2):cases.append((f,rep,'autonomous',1200,48,56,'initial_coherent',{'m0':.91}))
    for tag,p in [('slow_activation',{'gamma':.001}),('slow_propagation',{'beta':.04}),('lower_classical_gain',{'classical_efficiency':.85}),('higher_classical_gain',{'classical_efficiency':1.0}),('higher_activity_scale',{'activity_scale':.15})]:
        cases.append((f,0,'autonomous',1200,48,56,tag,p))
records=[]
with ThreadPoolExecutor(max_workers=2) as pool:
    for future in as_completed([pool.submit(execute,*args) for args in cases]):
        r=future.result();records.append(r);print(json.dumps(r),flush=True)
        (ROOT/'validation'/'extension_execution.json').write_text(json.dumps(records,indent=2)+'\n')
