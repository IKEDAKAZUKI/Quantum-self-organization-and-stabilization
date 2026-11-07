"""Periodic morphology classification from occupation-contrast configurations.

Spectral nematic order, compact positive components and signed-domain geometry
implement the fixed criteria documented in docs/ANALYSIS_DETAILS.md.
"""
import numpy as np
from scipy.ndimage import label, maximum_filter
from scipy.spatial import cKDTree
PARAMS = dict(power_cut=0.005, shell_fraction_cut=0.4, nematic_cut=0.85, level_rms=0.2, compact_roundness=0.38, compact_aspect=2.3, spot_area_share=0.65, lab_area_share=0.4)

def periodic_labels(mask):
    a, n = label(mask, np.ones((3, 3), int))
    par = np.arange(n + 1)

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    def join(x, y):
        if x and y:
            x, y = (find(x), find(y))
            if x != y:
                par[max(x, y)] = min(x, y)
    for k in [-1, 0, 1]:
        for x, y in zip(a[0], np.roll(a[-1], k)):
            join(x, y)
        for x, y in zip(a[:, 0], np.roll(a[:, -1], k)):
            join(x, y)
    lut = np.array([find(x) for x in range(n + 1)])
    a = lut[a]
    vals = np.unique(a)
    vals = vals[vals > 0]
    out = np.zeros_like(a)
    for i, v in enumerate(vals, 1):
        out[a == v] = i
    return (out, len(vals))

def components(v, level):
    a, n = periodic_labels(v > level)
    ny, nx = v.shape
    total = (a > 0).sum()
    items = []
    for i in range(1, n + 1):
        m = a == i
        y, x = np.where(m)
        area = len(x)
        per = sum((int(np.sum(m & ~np.roll(m, d, axis=ax))) for ax in (0, 1) for d in (-1, 1)))
        cx = np.angle(np.exp(2j * np.pi * x / nx).mean()) * nx / (2 * np.pi)
        cy = np.angle(np.exp(2j * np.pi * y / ny).mean()) * ny / (2 * np.pi)
        dx = (x - cx + nx / 2) % nx - nx / 2
        dy = (y - cy + ny / 2) % ny - ny / 2
        if area >= 3:
            cov = np.cov(np.c_[dx, dy].T)
            ev = np.linalg.eigvalsh(cov)
            aspect = np.sqrt((ev[-1] + 0.25) / (ev[0] + 0.25))
        else:
            aspect = 1.0
        circ = 4 * np.pi * area / max(per * per, 1)
        compact = area >= 5 and circ >= PARAMS['compact_roundness'] and (aspect <= PARAMS['compact_aspect'])
        items.append((area, per, circ, aspect, compact))
    share = sum((i[0] for i in items if i[-1])) / max(total, 1)
    return dict(ncomponents=n, ncompact=sum((i[-1] for i in items)), compact_area_share=share, largest_frac=max([i[0] for i in items], default=0) / v.size, area_median=float(np.median([i[0] for i in items])) if items else 0.0)

def features(u, p, power=None, nematic=None):
    v = u - u.mean()
    ny, nx = u.shape
    M = u.size
    P = abs(np.fft.fft2(v) / M) ** 2
    var = np.mean(v * v)
    wx = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(nx))
    wy = 2 - 2 * np.cos(2 * np.pi * np.fft.fftfreq(ny))
    w = wx[None, :] + wy[:, None]
    mask = -p['delta'] + (p['D2'] + p['Dcl']) * w - p['D4'] * w * w - p['zeta'] * wy[:, None] > 0
    mask[0, 0] = False
    band = P[mask].sum() if power is None else power
    ang = np.arctan2(np.fft.fftfreq(ny)[:, None], np.fft.fftfreq(nx)[None, :])
    nem = abs((P[mask] * np.exp(2j * ang[mask])).sum()) / max(band, 1e-30) if nematic is None else nematic
    rms = np.sqrt(var)
    co = components(v, PARAMS['level_rms'] * rms)
    neg = components(-v, PARAMS['level_rms'] * rms)
    peak = (v == maximum_filter(v, size=5, mode='wrap')) & (v > 0.5 * rms)
    y, x = np.where(peak)
    xy = np.c_[x, y]
    psi = np.nan
    local = np.nan
    cv = np.nan
    if len(xy) >= 7:
        dist, ne = cKDTree(xy, boxsize=[nx, ny]).query(xy, k=7)
        dv = xy[ne[:, 1:]] - xy[:, None, :]
        dv = dv - np.rint(dv / [nx, ny]) * [nx, ny]
        z = np.exp(6j * np.arctan2(dv[:, :, 1], dv[:, :, 0])).mean(axis=1)
        psi = abs(z.mean())
        local = abs(z).mean()
        cv = dist[:, 1:].std() / dist[:, 1:].mean()
    shell = band / max(var, 1e-30)
    resolved = band >= PARAMS['power_cut'] and shell >= PARAMS['shell_fraction_cut']
    if not resolved:
        cl = 'unresolved'
    elif nem >= PARAMS['nematic_cut']:
        cl = 'stripe'
    elif co['compact_area_share'] >= PARAMS['spot_area_share'] and co['ncompact'] >= 4 and (co['largest_frac'] < 0.12):
        cl = 'spot'
    elif nem < 0.85 and co['compact_area_share'] < PARAMS['lab_area_share']:
        cl = 'labyrinth'
    else:
        cl = 'mixed'
    return dict(morphology=cl, resolved=bool(resolved), rms=rms, shell_fraction=shell, nematic=nem, psi6=psi, psi6_local=local, spacing_cv=cv, npeaks=len(xy), **co, negative_largest_frac=neg['largest_frac'])
