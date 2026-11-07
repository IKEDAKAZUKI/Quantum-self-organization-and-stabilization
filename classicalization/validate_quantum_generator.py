#!/usr/bin/env python3
"""Independent density-matrix check of the autonomous Lindblad construction.

The check uses one physical bosonic dimer and an environmental flag. Fixed
neighbor occupation/coherence values supply a nonzero local transport stencil.
It compares every basis-population derivative with the full quantum generator.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parent

def dissipator(A,rho):
    H=A.conj().T@A
    return A@rho@A.conj().T-(H@rho+rho@H)/2

def main():
    N=2;nf=N+2;df=nf*nf;dim=2*df
    def ket(a,b,h):
        z=np.zeros(dim,complex);z[h*df+a*nf+b]=1;return z
    def proj(z):return np.outer(z,z.conj())
    P=[];labels=[]
    for h in [0,1]:
        for n in range(N+1):
            L=ket(n,n+1,h);R=ket(n+1,n,h);Z=ket(n,n,h)
            if h==0:
                for s in [-1,1]:P.append(proj((L+s*R)/np.sqrt(2)));labels.append((h,n,s))
            else:P.append((proj(L)+proj(R))/2);labels.append((h,n,1))
            P.append(proj(Z));labels.append((h,n,0))
    ids={lab:i for i,lab in enumerate(labels)};B=np.stack(P)
    eps=.2;eta=1.;r=20.;kl=.02;kr=2.;eff=.91
    activation=.08;neighbor_gain=.6;Db=1e-3;smooth=.002
    C=np.zeros((dim,dim),complex);U=np.zeros_like(C);Hp=np.zeros_like(C)
    Tplus=np.zeros_like(C);source=np.zeros_like(C)
    leakL=np.zeros_like(C);leakR=np.zeros_like(C);repL=np.zeros_like(C);repR=np.zeros_like(C)
    jumps=[]
    for h in [0,1]:
        for a in range(nf):
            for b in range(nf):
                v=ket(a,b,h);U+=eps*((a+b-1)/2-1)*proj(v);Hp+=h*proj(v)
        for n in range(N+1):
            L=ket(n,n+1,h);R=ket(n+1,n,h);Z=ket(n,n,h)
            C+=np.outer(L,R.conj())+np.outer(R,L.conj())
            if h==0:
                plus=(L+R)/np.sqrt(2);minus=(L-R)/np.sqrt(2)
                source+=np.sqrt(r)*np.outer(plus,minus.conj())
            leakL+=np.sqrt(kl)*np.outer(Z,L.conj());leakR+=np.sqrt(kl)*np.outer(Z,R.conj())
            repL+=np.sqrt(kr/2)*np.outer(L,Z.conj());repR+=np.sqrt(kr/2)*np.outer(R,Z.conj())
            if n<N:
                for a,b in [(n,n+1),(n+1,n),(n,n)]:Tplus+=np.outer(ket(a+1,b+1,h),ket(a,b,h).conj())
    eye=np.eye(dim);H0=eye-Hp
    # A nonlinear occupation operator in the environment's activation channel.
    S=(U@U+.09*eye)/2
    aa=np.diag(S).real
    act=activation*aa**4/(aa**4+.04**4)+.015
    Aroot=np.diag(np.sqrt(act))
    for n in range(N+1):
        for a,b in [(n,n+1),(n+1,n),(n,n)]:
            J=np.outer(ket(a,b,1),ket(a,b,0).conj())@Aroot;jumps.append(J)
    na=np.diag([a for h in [0,1] for a in range(nf) for b in range(nf)])
    nb=np.diag([b for h in [0,1] for a in range(nf) for b in range(nf)])
    jumps.extend([source,np.sqrt(2*eta)*na,np.sqrt(2*eta)*nb,leakL,leakR,repL,repR])
    gain=H0@C+eff*Hp
    F=-.3*U-U@U@U+.7*gain@U+.1*eye
    values,V=np.linalg.eigh(F)
    for sign,T in [(1,Tplus),(-1,Tplus.conj().T)]:
        rates=(np.hypot(values,smooth)+sign*values)/(2*eps)+Db/eps**2
        root=(V*np.sqrt(rates))@V.conj().T;jumps.append(T@root)
    Q=np.zeros((len(B),len(B)))
    def add(col,label,rate):
        if rate:Q[ids[label],col]+=rate;Q[col,col]-=rate
    for j,(h,n,s) in enumerate(labels):
        u=eps*(n-1-(s==0)/2);g=eff if h else s
        f=-.3*u-u**3+.7*g*u+.1
        for direction in [-1,1]:
            if 0<=n+direction<=N:
                q=(np.hypot(f,smooth)+direction*f)/(2*eps)+Db/eps**2
                add(j,(h,n+direction,s),q)
        if s==0:
            if h:add(j,(h,n,1),kr)
            else:
                add(j,(h,n,-1),kr/2);add(j,(h,n,1),kr/2)
        else:
            add(j,(h,n,0),kl)
            if not h:add(j,(h,n,-s),eta+(r if s==-1 else 0))
        if not h:
            a=(u*u+.09)/2
            add(j,(1,n,0 if s==0 else 1),activation*a**4/(a**4+.04**4)+.015)
    errors=[]
    for j,rho in enumerate(B):
        exact=sum((dissipator(J,rho) for J in jumps),np.zeros_like(rho))
        pop=np.einsum('i,ijk->jk',Q[:,j],B)
        errors.append(float(np.max(np.abs(exact-pop))))
    rho0=np.zeros(len(B));rho0[ids[(0,1,-1)]]=rho0[ids[(0,1,1)]]=.5
    rows=[]
    for t in [0,.05,.2,1,10,50,100,200,400]:
        weights=expm(Q*t)@rho0
        joint=np.einsum('i,ijk->jk',weights,B)
        physical=joint[:df,:df]+joint[df:,df:]
        pt=physical.reshape(nf,nf,nf,nf).transpose(0,3,2,1).reshape(df,df)
        neg=float(np.maximum(-np.linalg.eigvalsh(pt),0).sum())
        inactive=sum(w for w,(h,n,s) in zip(weights,labels) if h==0 and s!=0)
        coherence=float(np.trace(C@joint).real)
        leakage=sum(w for w,(h,n,s) in zip(weights,labels) if s==0)
        rows.append(dict(t=t,negativity=neg,negativity_lower=max(0,coherence-2*leakage)/4,negativity_upper=inactive/2,
                         trace=float(np.trace(joint).real),min_eigenvalue=float(np.linalg.eigvalsh(joint).min())))
    env_adjoint=sum((J.conj().T@U@J-(J.conj().T@J@U+U@J.conj().T@J)/2 for J in jumps[:3*(N+1)]),np.zeros_like(U))
    assert np.max(abs(env_adjoint))<1e-12
    result=dict(activation_occupation_residual=float(np.max(abs(env_adjoint))),local_hilbert_dimension=dim,population_states=len(B),max_generator_residual=max(errors),max_column_sum=float(abs(Q.sum(0)).max()),max_negative_probability=max(0.,float(-weights.min())),evolution=rows)
    assert max(errors)<1e-11,result
    assert all(x['negativity']<=x['negativity_upper']+1e-10 and x['negativity']>=x['negativity_lower']-1e-10 for x in rows)
    (ROOT/'validation'/'density_matrix_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
