"""Fixed morphology and periodic spatial-extent observables.

All trials contribute to the ensemble; a failed target has pattern intensity zero.
The thresholds are the ones used to generate the reported confirmation statistics.
"""
import numpy as np
from numba import njit
from .morphology import features, PARAMS
from .observables import spectrum
GLOBAL_RULES = dict(stripe='at least one signed component has a noncontractible torus cycle', spot='compact positive component centers occupy at least six of nine equal tiles', labyrinth='one signed component has >=0.15 array area and a noncontractible torus cycle', labyrinth_min_area=0.15, spot_min_tiles=6, tiles_per_axis=3, connectivity=8, original_morphology=PARAMS)

@njit(cache=True)
def lifted_components(mask):
    ny, nx = mask.shape
    N = ny * nx
    labels = np.zeros((ny, nx), np.int32)
    X = np.zeros((ny, nx), np.int32)
    Y = np.zeros((ny, nx), np.int32)
    q = np.empty(N, np.int32)
    areas = np.zeros(N + 1, np.int32)
    wx = np.zeros(N + 1, np.bool_)
    wy = np.zeros(N + 1, np.bool_)
    lab = 0
    for y0 in range(ny):
        for x0 in range(nx):
            if not mask[y0, x0] or labels[y0, x0] != 0:
                continue
            lab += 1
            h = 0
            t = 1
            q[0] = y0 * nx + x0
            labels[y0, x0] = lab
            X[y0, x0] = x0
            Y[y0, x0] = y0
            while h < t:
                v = q[h]
                h += 1
                y = v // nx
                x = v % nx
                areas[lab] += 1
                for dy in range(-1, 2):
                    for dx in range(-1, 2):
                        if dy == 0 and dx == 0:
                            continue
                        yy = (y + dy) % ny
                        xx = (x + dx) % nx
                        if not mask[yy, xx]:
                            continue
                        tx = X[y, x] + dx
                        ty = Y[y, x] + dy
                        if labels[yy, xx] == 0:
                            labels[yy, xx] = lab
                            X[yy, xx] = tx
                            Y[yy, xx] = ty
                            q[t] = yy * nx + xx
                            t += 1
                        else:
                            if tx != X[yy, xx]:
                                wx[lab] = True
                            if ty != Y[yy, xx]:
                                wy[lab] = True
    return (labels, areas[:lab + 1], wx[:lab + 1], wy[:lab + 1])

def extent(u):
    v = np.asarray(u, float) - np.mean(u)
    rms = np.sqrt(np.mean(v * v))
    ny, nx = v.shape
    pos = v > PARAMS['level_rms'] * rms
    neg = -v > PARAMS['level_rms'] * rms
    a, areas, wx, wy = lifted_components(pos)
    b, barea, bwx, bwy = lifted_components(neg)
    wrapping = bool(np.any(wx[1:] | wy[1:]) or np.any(bwx[1:] | bwy[1:]))
    wa = np.r_[areas[1:][wx[1:] | wy[1:]], barea[1:][bwx[1:] | bwy[1:]]]
    largest_wrap = float(np.max(wa) / u.size) if wa.size else 0.0
    compact_tiles = set()
    compact_count = 0
    for i in range(1, len(areas)):
        m = a == i
        y, x = np.where(m)
        area = len(x)
        per = sum((int(np.sum(m & ~np.roll(m, d, axis=ax))) for ax in (0, 1) for d in (-1, 1)))
        cx = np.angle(np.exp(2j * np.pi * x / nx).mean()) % (2 * np.pi) * nx / (2 * np.pi)
        cy = np.angle(np.exp(2j * np.pi * y / ny).mean()) % (2 * np.pi) * ny / (2 * np.pi)
        dx = (x - cx + nx / 2) % nx - nx / 2
        dy = (y - cy + ny / 2) % ny - ny / 2
        if area >= 3:
            ev = np.linalg.eigvalsh(np.cov(np.c_[dx, dy].T))
            aspect = np.sqrt((ev[-1] + 0.25) / (ev[0] + 0.25))
        else:
            aspect = 1.0
        circ = 4 * np.pi * area / max(per * per, 1)
        if area >= 5 and circ >= PARAMS['compact_roundness'] and (aspect <= PARAMS['compact_aspect']):
            compact_count += 1
            compact_tiles.add((min(2, int(3 * cx / nx)), min(2, int(3 * cy / ny))))
    return dict(any_wrap=wrapping, largest_wrapping_fraction=largest_wrap, compact_tiles=len(compact_tiles), compact_count=compact_count)

def measure(u, p, family=None):
    f = features(u, p)
    ex = extent(u)
    f.update(ex)
    f['power'] = float(spectrum(u, p)['power'])
    for g in ('stripe', 'spot', 'labyrinth'):
        good = ex['any_wrap'] if g == 'stripe' else ex['compact_tiles'] >= 6 if g == 'spot' else ex['largest_wrapping_fraction'] >= 0.15
        hit = f['morphology'] == g and good
        f['chi_' + g] = int(hit)
        f['G_' + g] = float(hit) * f['power']
    if family is not None:
        f['chi'] = f['chi_' + family]
        f['G'] = f['G_' + family]
    return f
