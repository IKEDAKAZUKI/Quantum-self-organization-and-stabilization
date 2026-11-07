"""Independent physical-number pinching check, 81-dimensional Hilbert space.
Build microscopic jump matrices in the +-0 basis, transform to LR0, apply full
number pinching, and compare with independently enumerated phase-averaged rates.
Also check the five-cell local field sum with leakage labels.
"""
from pathlib import Path
import sys, itertools, json, hashlib
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix, eye, kron, diags
from scipy.sparse.linalg import expm_multiply
ROOT = Path(__file__).resolve().parents[1]
from .density import liouv

def full_projected():
    levels = 3
    eps = 0.1
    delta = 0.3
    D4 = 0.04
    D2 = 0.12
    nu = 0.15
    b = 1.0
    c = 0.1
    r = 2.0
    eta = 0.3
    loss = 0.02
    repair = 0.8
    Db = 0.0003
    sm = 0.01
    st = list(itertools.product(range(levels), [1, -1, 0]))
    ds = len(st)
    ix = {v: i for i, v in enumerate(st)}
    I = eye(ds, format='csr')
    R = np.zeros((ds, ds))
    T = R.copy()
    na = R.copy()
    nb = R.copy()
    VL = R.copy()
    VR = R.copy()
    WL = R.copy()
    WR = R.copy()
    B = R.copy()
    for n in range(levels):
        p, m, z = (ix[n, a] for a in [1, -1, 0])
        R[p, m] = 1
        na[p, p] = na[m, m] = nb[p, p] = nb[m, m] = n + 0.5
        na[z, z] = nb[z, z] = n
        na[p, m] = na[m, p] = -0.5
        nb[p, m] = nb[m, p] = 0.5
        VL[z, p] = VL[z, m] = 1 / np.sqrt(2)
        VR[z, p] = 1 / np.sqrt(2)
        VR[z, m] = -1 / np.sqrt(2)
        WL[p, z] = WL[m, z] = 1 / np.sqrt(2)
        WR[p, z] = 1 / np.sqrt(2)
        WR[m, z] = -1 / np.sqrt(2)
        B[p, p] = B[m, p] = B[p, m] = 1 / np.sqrt(2)
        B[m, m] = -1 / np.sqrt(2)
        B[z, z] = 1
        if n < levels - 1:
            for s in [1, -1, 0]:
                T[ix[n + 1, s], ix[n, s]] = 1

    def loc(a, x):
        return kron(csr_matrix(a), I) if x == 0 else kron(I, csr_matrix(a))
    configs = list(itertools.product(range(ds), repeat=2))
    D = ds * ds
    js = []
    vals = np.array([[eps * (st[i][0] - 1 - (0.5 if st[i][1] == 0 else 0)) for i in cfg] for cfg in configs])
    sg = np.array([[st[i][1] for i in cfg] for cfg in configs])
    lap = np.array([[1.0, -1.0], [-1.0, 1.0]])
    fs = -delta * vals - D4 * (vals @ lap @ lap) + D2 * (sg[:, 0] + sg[:, 1])[:, None] / 2 * (vals @ lap) + nu * vals ** 2 - b * vals ** 3 - c * vals ** 5
    for x in [0, 1]:
        for a in [np.sqrt(r) * R, np.sqrt(2 * eta) * na, np.sqrt(2 * eta) * nb, np.sqrt(loss) * VL, np.sqrt(loss) * VR, np.sqrt(repair / 2) * WL, np.sqrt(repair / 2) * WR]:
            js.append(loc(a, x))
        f = fs[:, x]
        up = (np.hypot(f, sm) + f) / (2 * eps)
        dn = (np.hypot(f, sm) - f) / (2 * eps)
        js.extend([loc(T, x) @ diags(np.sqrt(up)), loc(T.T, x) @ diags(np.sqrt(dn)), np.sqrt(Db) / eps * loc(T, x), np.sqrt(Db) / eps * loc(T.T, x)])
    BB = kron(csr_matrix(B), csr_matrix(B))
    jn = [BB.T @ j @ BB for j in js]
    L = liouv(jn)
    diagids = np.arange(D) * (D + 1)
    pinch = diags(np.isin(np.arange(D * D), diagids).astype(float))
    LP = (pinch @ L @ pinch).tocsr()
    K = L[diagids][:, diagids].real.toarray()
    cs = list(itertools.product(range(levels), [1, 0]))
    ci = {x: i for i, x in enumerate(cs)}
    dc = len(cs)
    C = dc * dc
    A = np.zeros((C, D))
    expander = np.zeros((D, C))
    for j, cfg in enumerate(configs):
        pair = [(st[q][0], int(st[q][1] != 0)) for q in cfg]
        row = ci[pair[0]] * dc + ci[pair[1]]
        A[row, j] = 1
        expander[j, row] = 2.0 ** (-sum((q[1] for q in pair)))
    Q = np.zeros((C, C))
    for j, (ia, ib) in enumerate(itertools.product(range(dc), repeat=2)):
        pair = [cs[ia], cs[ib]]
        u = np.array([eps * (n - 1 - (0.5 if active == 0 else 0)) for n, active in pair])
        f0 = -delta * u - D4 * (lap @ lap @ u) + nu * u * u - b * u ** 3 - c * u ** 5
        signs = list(itertools.product(*[[1, -1] if active else [0] for n, active in pair]))
        for x in [0, 1]:
            f = np.array([f0[x] + D2 * (s[0] + s[1]) / 2 * (lap @ u)[x] for s in signs])
            up = np.mean((np.hypot(f, sm) + f) / (2 * eps)) + Db / eps ** 2
            dn = np.mean((np.hypot(f, sm) - f) / (2 * eps)) + Db / eps ** 2
            n, active = pair[x]
            events = []
            if n < levels - 1:
                events.append(((n + 1, active), up))
            if n > 0:
                events.append(((n - 1, active), dn))
            events.append(((n, 1 - active), loss if active else repair))
            for target, rate in events:
                new = pair.copy()
                new[x] = target
                k = ci[new[0]] * dc + ci[new[1]]
                Q[k, j] += rate
                Q[j, j] -= rate
    rng = np.random.default_rng(66321)
    p0 = rng.random(D)
    p0 /= p0.sum()
    z0 = A @ p0
    ts = np.linspace(0, 2, 21)
    full = expm_multiply(LP, np.diag(p0).ravel(order='F'), start=0, stop=2, num=21)
    red = expm_multiply(csr_matrix(Q), z0, start=0, stop=2, num=21)
    rows = []
    for t, f, z in zip(ts, full, red):
        rho = f.reshape((D, D), order='F')
        p = np.real(np.diag(rho))
        rows.append(dict(t=t, marginal_error=float(np.max(abs(A @ p - z))), offdiagonal_residual=float(np.max(abs(rho - np.diag(p)))), trace_error=float(abs(p.sum() - 1))))
    pd.DataFrame(rows).to_csv(ROOT / 'validation/projected_density_timeseries.csv', index=False)
    return dict(hilbert_dimension=D, coarse_dimension=C, lumpability_residual=float(abs(A @ K - Q @ A).max()), column_sum_residual=float(abs(Q.sum(axis=0)).max()), max_marginal_error=max((x['marginal_error'] for x in rows)), max_offdiagonal=max((x['offdiagonal_residual'] for x in rows)))

def local_five():
    rng = np.random.default_rng(210621)
    rows = []
    for q in range(200):
        u = rng.uniform(-0.8, 0.8, 5)
        active = rng.random(5) > 0.2
        weights = rng.uniform(1.3, 1.6, 4)
        f0 = rng.uniform(-2, 2)
        eps = 0.01
        sm = 0.002
        fs = []
        for mask in range(32):
            sg = np.array([1 if mask >> i & 1 else -1 for i in range(5)]) * active
            fs.append(f0 + np.dot(weights, 0.5 * (sg[0] + sg[1:]) * (u[0] - u[1:])))
        fs = np.array(fs)
        up = np.mean((np.hypot(fs, sm) + fs) / (2 * eps))
        dn = np.mean((np.hypot(fs, sm) - fs) / (2 * eps))
        unique = []
        for sg in itertools.product(*[[1, -1] if a else [0] for a in active]):
            sg = np.array(sg)
            unique.append(f0 + np.dot(weights, 0.5 * (sg[0] + sg[1:]) * (u[0] - u[1:])))
        ff = np.array(unique)
        eu = np.mean((np.hypot(ff, sm) + ff) / (2 * eps))
        ed = np.mean((np.hypot(ff, sm) - ff) / (2 * eps))
        rows.append(dict(case=q, active_sites=int(active.sum()), rate_difference=max(abs(up - eu), abs(dn - ed)), drift_error=abs(eps * (up - dn) - f0), projected_variance=eps * eps * (up + dn), substituted_variance=eps * np.hypot(f0, sm)))
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / 'validation/projected_local_rate_tests.csv', index=False)
    return dict(cases=len(d), max_rate_error=float(d.rate_difference.max()), max_drift_error=float(d.drift_error.max()), max_variance_ratio=float((d.projected_variance / d.substituted_variance).max()))
