#!/usr/bin/env python3
"""Run a fixed-generator environmental classicalization study from weak fields."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import hashlib, json, subprocess, time
import numpy as np

ROOT=Path(__file__).resolve().parent
FAMILIES=('stripe','spot','labyrinth')
BASE=dict(nx=48,ny=56,eps=.01,U=1.,delta=.3,D4=1.,b=1.,c=.1,r=20.,eta=1.,leak=.0002,repair=2.,Db=1e-5,smooth=.002,Dcl=.001,disorder=.05,gamma=.002,beta=.08,activity_scale=.10,classical_efficiency=.91,quantum_enabled=1,environment_initial=0,m0=0.,verify_every=1000000)
FAMILY=dict(stripe=dict(D2=1.6,zeta=.8,nu=0.),spot=dict(D2=1.4,zeta=0.,nu=.55),labyrinth=dict(D2=1.8,zeta=0.,nu=0.))

def convert(src:Path,dest:Path,eps:float):
    with src.open('rb') as f:
        nx,ny,m=map(int,np.fromfile(f,'<i4',3))
        dtype=np.dtype([('t','<f8'),('n','<i4',(m,)),('s','i1',(m,)),('h','i1',(m,))])
        if nx*ny!=m or (src.stat().st_size-12)%dtype.itemsize: raise ValueError('Malformed trajectory')
        a=np.fromfile(f,dtype)
    np.savez_compressed(dest,t=a['t'],n=a['n'].reshape(-1,ny,nx),s=a['s'].reshape(-1,ny,nx),h=a['h'].reshape(-1,ny,nx),eps=eps)
    src.unlink()

def execute(family:str,rep:int,protocol:str,T:float,nx:int,ny:int,tag:str,overrides:dict)->dict:
    p=BASE|FAMILY[family]|dict(nx=nx,ny=ny,T=T,save=5.,seed=931000+1000*FAMILIES.index(family)+rep,disseed=73000+rep)|overrides
    if protocol=='no_quantum':p['quantum_enabled']=0;p['r']=0
    elif protocol=='no_classical_gain':p['classical_efficiency']=0
    elif protocol=='no_environment':p['gamma']=p['beta']=0
    elif protocol=='projected_control':p['quantum_enabled']=0
    elif protocol=='classical_ready':p['quantum_enabled']=0;p['r']=0;p['environment_initial']=1
    elif protocol!='autonomous':raise ValueError(protocol)
    name=f'{tag}_{protocol}_{family}_{rep:02d}'
    out=ROOT/'trajectories'/name
    if out.with_suffix('.npz').is_file():
        saved=json.loads(Path(str(out)+'_params.json').read_text())
        different=[key for key,value in p.items() if saved.get(key)!=value]
        if different:
            raise ValueError(f'{name}: existing run differs in {different}; choose a new --tag.')
        return dict(id=name,status='existing')
    initial=ROOT/'initial_states'/f'weak_{nx}x{ny}_{rep:02d}.txt'
    if not initial.exists():
        rng=np.random.default_rng(812300+rep)
        n=np.rint(rng.normal(0,.02,(ny,nx))/p['eps']).astype(int)
        np.savetxt(initial,np.column_stack([n.ravel(),np.full(n.size,9)]),fmt='%d')
    p['init']='../initial_states/'+initial.name
    cfg=ROOT/'configs'/(name+'.cfg');cfg.write_text(''.join(f'{k} {v}\n' for k,v in p.items()))
    source_name='projected_control' if protocol=='projected_control' else 'autonomous_jump'
    metadata=p|dict(family=family,rep=rep,protocol=protocol,tag=tag,initial_sha256=hashlib.sha256(initial.read_bytes()).hexdigest(),source_sha256=hashlib.sha256((ROOT/'src'/(source_name+'.cpp')).read_bytes()).hexdigest())
    Path(str(out)+'_params.json').write_text(json.dumps(metadata,indent=2)+'\n')
    t=time.monotonic();run=subprocess.run([str(ROOT/'build'/source_name),str(cfg),str(out)],capture_output=True,text=True)
    if run.returncode:raise RuntimeError(name+': '+run.stderr)
    convert(out.with_suffix('.bin'),out.with_suffix('.npz'),p['eps'])
    result=json.loads(Path(str(out)+'_engine.json').read_text())|dict(id=name,protocol=protocol,family=family,rep=rep,wall_seconds=time.monotonic()-t)
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--families',nargs='+',default=list(FAMILIES),choices=FAMILIES)
    p.add_argument('--protocols',nargs='+',default=['autonomous'])
    p.add_argument('--replicates',type=int,default=4);p.add_argument('--duration',type=float,default=1200)
    p.add_argument('--nx',type=int,default=48);p.add_argument('--ny',type=int,default=56)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--tag',default='main')
    p.add_argument('--overrides',default='{}',help='JSON object of physical parameter overrides')
    args=p.parse_args();overrides=json.loads(args.overrides)
    if args.replicates<1 or args.duration<=0 or args.workers<1:
        raise ValueError('Require positive replicates, duration, and workers.')
    records=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs=[pool.submit(execute,f,i,pr,args.duration,args.nx,args.ny,args.tag,overrides) for pr in args.protocols for f in args.families for i in range(args.replicates)]
        for future in as_completed(jobs):
            x=future.result();records.append(x);print(json.dumps(x),flush=True)
            (ROOT/'validation'/f'{args.tag}_execution.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
