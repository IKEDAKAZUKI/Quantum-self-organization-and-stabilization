"""Independent weighted partial-transpose and density-matrix checks."""
from pathlib import Path
import json, itertools
import numpy as np
from scipy.sparse import csr_matrix, eye, kron, diags
from scipy.sparse.linalg import expm_multiply
from scipy.linalg import eigvalsh
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]

def dissipator(j):
    j = csr_matrix(j, dtype=complex)
    d = j.shape[0]
    I = eye(d, format='csr')
    K = j.conj().T @ j
    return kron(j.conjugate(), j) - 0.5 * kron(I, K) - 0.5 * kron(K.T, I)

def liouv(jumps, H=None):
    out = sum((dissipator(j) for j in jumps))
    if H is not None:
        h = csr_matrix(H)
        I = eye(h.shape[0], format='csr')
        out += -1j * (kron(I, h) - kron(h.T, I))
    return out.tocsr()

def witness_check():
    d = 9
    eps = 0.08
    n0 = 4
    D = d * d
    u = np.array([eps * ((a + b - 1) / 2 - n0) for a in range(d) for b in range(d)])
    P = np.diag([int(a == b) for a in range(d) for b in range(d)])
    C = np.zeros((D, D))
    for n in range(d - 1):
        C[n * d + n + 1, (n + 1) * d + n] = 1
        C[(n + 1) * d + n, n * d + n + 1] = 1
    rng = np.random.default_rng(574)
    rows = []
    for j in range(30):
        neighbors = rng.uniform(-0.6, 0.6, 4)
        f = ((u[:, None] - neighbors) ** 2).sum(axis=1)
        W = P @ np.diag(2 * f + 4 * eps ** 2 / 2) - C @ np.diag(f)
        pt = W.reshape(d, d, d, d).transpose(0, 3, 2, 1).reshape(D, D)
        rows.append({'case': j, 'min_partial_transpose_eigenvalue': float(eigvalsh(pt)[0]), 'hermiticity': float(abs(pt - pt.T).max())})
    pd.DataFrame(rows).to_csv(ROOT / 'validation/weighted_leak_witness.csv', index=False)
    return min((r['min_partial_transpose_eigenvalue'] for r in rows))

def full_check():
    nlevels = 3
    eps = 0.1
    delta = 0.3
    D4 = 0.04
    D2 = 0.12
    nu = 0.15
    b = 1
    c = 0.1
    r = 2
    eta = 0.3
    loss = 0.02
    repair = 0.8
    Db = 0.0003
    sm = 0.01
    states = list(itertools.product(range(nlevels), [1, -1, 0]))
    ds = len(states)
    ind = {x: i for i, x in enumerate(states)}
    u = np.array([eps * (n - 1 - (0.5 if s == 0 else 0)) for n, s in states])
    s = np.array([s for n, s in states])
    R = np.zeros((ds, ds))
    T = np.zeros_like(R)
    na = np.zeros_like(R)
    nb = np.zeros_like(R)
    VL = np.zeros_like(R)
    VR = np.zeros_like(R)
    WL = np.zeros_like(R)
    WR = np.zeros_like(R)
    for n in range(nlevels):
        plus, minus, eq = (ind[n, x] for x in [1, -1, 0])
        R[plus, minus] = 1
        na[plus, plus] = na[minus, minus] = nb[plus, plus] = nb[minus, minus] = n + 0.5
        na[eq, eq] = nb[eq, eq] = n
        na[plus, minus] = na[minus, plus] = -0.5
        nb[plus, minus] = nb[minus, plus] = 0.5
        VL[eq, plus] = VL[eq, minus] = 1 / np.sqrt(2)
        VR[eq, plus] = 1 / np.sqrt(2)
        VR[eq, minus] = -1 / np.sqrt(2)
        WL[plus, eq] = WL[minus, eq] = 1 / np.sqrt(2)
        WR[plus, eq] = 1 / np.sqrt(2)
        WR[minus, eq] = -1 / np.sqrt(2)
        if n < nlevels - 1:
            for sg in [1, -1, 0]:
                T[ind[n + 1, sg], ind[n, sg]] = 1
    I = eye(ds, format='csr')
    D = ds * ds
    jumps = []

    def loc(j, x):
        return kron(csr_matrix(j), I) if x == 0 else kron(I, csr_matrix(j))
    for x in [0, 1]:
        for j in [np.sqrt(r) * R, np.sqrt(2 * eta) * na, np.sqrt(2 * eta) * nb, np.sqrt(loss) * VL, np.sqrt(loss) * VR, np.sqrt(repair / 2) * WL, np.sqrt(repair / 2) * WR]:
            jumps.append(loc(j, x))
    config = list(itertools.product(range(ds), repeat=2))
    G = np.zeros((D, D))
    for x in [0, 1]:
        ff = []
        for i, j in config:
            vals = np.array([u[i], u[j]])
            sign = [s[i], s[j]]
            L = np.array([[1.0, -1.0], [-1.0, 1.0]])
            f = -delta * vals - D4 * (L @ L @ vals) + D2 * 0.5 * (sign[0] + sign[1]) * (L @ vals) + nu * vals ** 2 - b * vals ** 3 - c * vals ** 5
            ff.append(f[x])
        ff = np.asarray(ff)
        rup = (np.hypot(ff, sm) + ff) / (2 * eps)
        rdn = (np.hypot(ff, sm) - ff) / (2 * eps)
        jumps.append(loc(T, x) @ diags(np.sqrt(rup)))
        jumps.append(loc(T.T, x) @ diags(np.sqrt(rdn)))
        jumps.append(np.sqrt(Db) / eps * loc(T, x))
        jumps.append(np.sqrt(Db) / eps * loc(T.T, x))
        for j, (ia, ib) in enumerate(config):
            cur = [ia, ib]
            n, sg = states[cur[x]]
            events = []
            if n < nlevels - 1:
                events.append(((n + 1, sg), rup[j] + Db / eps ** 2))
            if n > 0:
                events.append(((n - 1, sg), rdn[j] + Db / eps ** 2))
            if sg == 0:
                events.extend([((n, 1), repair / 2), ((n, -1), repair / 2)])
            else:
                events.extend([((n, -sg), eta + (r if sg == -1 else 0)), ((n, 0), loss)])
            for target, rate in events:
                c2 = cur.copy()
                c2[x] = ind[target]
                k = c2[0] * ds + c2[1]
                G[k, j] += rate
                G[j, j] -= rate
    L = liouv(jumps)
    p0 = np.zeros(D)
    for sa, sb in itertools.product([1, -1], repeat=2):
        p0[ind[1, sa] * ds + ind[2, sb]] = 0.25
    rho0 = np.diag(p0)
    ts = np.linspace(0, 2, 21)
    sol = expm_multiply(L, rho0.ravel(order='F'), start=0, stop=2, num=21)
    ps = expm_multiply(csr_matrix(G), p0, start=0, stop=2, num=21)
    rows = []
    for t, y, p in zip(ts, sol, ps):
        rho = y.reshape(D, D, order='F')
        rows.append(dict(t=t, max_density_difference=float(abs(rho - np.diag(p)).max()), trace_error=float(abs(np.trace(rho) - 1)), min_eigenvalue=float(eigvalsh((rho + rho.conj().T) / 2)[0])))
    pd.DataFrame(rows).to_csv(ROOT / 'validation/full_density_reduction.csv', index=False)
    return {'hilbert_dimension': D, 'max_density_difference': max((x['max_density_difference'] for x in rows)), 'min_eigenvalue': min((x['min_eigenvalue'] for x in rows)), 'generator_trace_residual': float(abs(G.sum(axis=0)).max())}
