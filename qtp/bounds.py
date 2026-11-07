"""Separable formation and maintenance bounds in Eqs. (28), (30), (38)--(42).

Times are dimensionless. Parameters use the saved simulation names: ``Db`` is
D_bg, ``smooth`` is s_0, and ``leak``/``repair`` are kappa_l/kappa_r.
"""
from __future__ import annotations
from typing import Mapping
import numpy as np

Parameters = Mapping[str, float]


def _times(t: object) -> np.ndarray:
    out = np.asarray(t, dtype=float)
    if np.any(~np.isfinite(out)) or np.any(out < 0):
        raise ValueError("Times must be finite and nonnegative.")
    return out


def decay_convolution(t: object, alpha: float, k: float) -> np.ndarray:
    """Return (exp(-k*t)-exp(-alpha*t))/(alpha-k), including equal rates.

    This evaluates the convolution of two decays without subtracting nearly
    equal exponentials. The min-rate form also avoids overflow for long times.
    At alpha == k the exact limit is t*exp(-alpha*t).
    """
    t = _times(t)
    if not np.isfinite(alpha) or not np.isfinite(k) or min(alpha, k) < 0:
        raise ValueError("Decay rates must be finite and nonnegative.")
    gap = abs(float(alpha) - float(k))
    if gap == 0.0:
        return t * np.exp(-alpha * t)
    with np.errstate(over="ignore", under="ignore"):
        return np.exp(-min(alpha, k) * t) * (-np.expm1(-gap * t)) / gap


def comparison_constants(p: Parameters, theta: float = 0.15) -> dict[str, float]:
    """Coefficients of dE/dt <= -alpha*E + c0 + c1*ell, Appendix B.2.

    Raise ValueError outside the coefficient assumptions of the stated bound.
    """
    required = ("eps", "U", "delta", "D2", "D4", "Dcl", "zeta", "nu", "b", "c",
                "disorder", "Db", "smooth", "leak", "repair")
    v = {name: float(p[name]) for name in required}
    if not all(np.isfinite(x) for x in v.values()):
        raise ValueError("All bound parameters must be finite.")
    if not 0 < theta < 1 or not 0 <= v["disorder"] < 1:
        raise ValueError("Require 0 < theta < 1 and 0 <= disorder < 1.")
    if min(v[x] for x in ("eps", "U", "delta", "D4", "b")) <= 0:
        raise ValueError("eps, U, delta, D4 and b must be positive.")
    if min(v[x] for x in ("D2", "Dcl", "zeta", "c", "Db", "smooth", "leak", "repair")) < 0:
        raise ValueError("Dissipative rates, nonlinear c and positive transport weights cannot be negative.")
    k = v["leak"] + v["repair"]
    if k <= 0:
        raise ValueError("The reported bound assumes leak + repair > 0.")
    eps, D4, b, c = v["eps"], v["D4"], v["b"], v["c"]
    dm = v["delta"] * (1 - v["disorder"])
    dp = v["delta"] * (1 + v["disorder"])
    J = v["D2"] * (1 + v["disorder"])
    nu = abs(v["nu"])
    d = 2 * dm - 16 * v["Dcl"] - nu * nu / (2 * b * (1 - theta)) - eps * nu
    if d <= 0:
        raise ValueError("Noncontracting comparison: the bound requires d > 0.")
    Cnl = 27*b*eps**4/(2048*theta**3) + c*5**5*eps**6/(6**6*2**5)
    Cdiff = eps**2 * (dp**2/d + (8*D4+v["Dcl"])**2/(4*D4)
                     + v["zeta"]/2 + 2*J**2/np.sqrt(d*D4)) + Cnl
    Ue = v["U"] + eps/2
    ell_inf = v["leak"]/k
    g = eps*Ue + eps*eps/4
    c0 = Cdiff + eps*v["smooth"] + 2*v["Db"] + g*v["leak"]
    c1 = 32*J*Ue*Ue + 2*J*eps*eps + g*v["repair"]
    alpha = d/2
    return dict(d=d, alpha=alpha, k=k, Cnl=Cnl, Cdiff=Cdiff, c0=c0, c1=c1,
                ell_inf=ell_inf, B_inf=(c0+c1*ell_inf)/alpha, J=J, Ue=Ue)


def formation_bound(t: object, p: Parameters, E0: float | None = None) -> np.ndarray:
    """Eq. (28): formation bound for initial mean squared contrast E0 and ell(0)=0."""
    t = _times(t)
    q = comparison_constants(p)
    E0 = float(p.get("E0", 0.0) if E0 is None else E0)
    if not np.isfinite(E0) or E0 < 0:
        raise ValueError("E0 must be finite and nonnegative.")
    return (E0*np.exp(-q["alpha"]*t) + q["B_inf"]*(-np.expm1(-q["alpha"]*t))
            - q["c1"]*q["ell_inf"]*decay_convolution(t, q["alpha"], q["k"]))


def maintenance_bound(tau: object, p: Parameters, ell_f: float) -> np.ndarray:
    """Eq. (30): maintenance bound conditional on the starting leakage fraction.

    ``tau`` is time after continuation starts. ``ell_f`` is the measured fraction
    of leaked cells in that starting configuration, not its stationary average.
    The conservative initial mean squared contrast is Ue**2, as in the article.
    """
    tau = _times(tau)
    ell_f = float(ell_f)
    if not np.isfinite(ell_f) or not 0 <= ell_f <= 1:
        raise ValueError("ell_f must be a finite leakage fraction in [0, 1].")
    q = comparison_constants(p)
    return (q["B_inf"] + (q["Ue"]**2-q["B_inf"])*np.exp(-q["alpha"]*tau)
            + q["c1"]*(ell_f-q["ell_inf"])*decay_convolution(tau, q["alpha"], q["k"]))


def source_moments(t: object, p: Parameters, e0: float = 0.0,
                   ell0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Mean coherence and leakage, Eqs. (12), (35)--(36), for constant rates."""
    t = _times(t)
    r, eta, leak, repair = (float(p[x]) for x in ("r", "eta", "leak", "repair"))
    if not all(np.isfinite(x) and x >= 0 for x in (r, eta, leak, repair)):
        raise ValueError("Source rates must be finite and nonnegative.")
    if not 0 <= ell0 <= 1 or abs(e0) > 1-ell0 + 1e-12:
        raise ValueError("Initial moments require 0 <= ell0 <= 1 and |e0| <= 1-ell0.")
    a, k = r + 2*eta + leak, leak + repair
    if k == 0:
        ell = np.full_like(t, ell0)
        e = (np.full_like(t, e0) if a == 0 else
             e0*np.exp(-a*t) + r*(1-ell0)*(-np.expm1(-a*t))/a)
    else:
        ell_inf = leak/k
        ell = ell0*np.exp(-k*t) + ell_inf*(-np.expm1(-k*t))
        if a == 0:
            e = np.full_like(t, e0)
        else:
            e = (e0*np.exp(-a*t) + r*(1-ell_inf)*(-np.expm1(-a*t))/a
                 - r*(ell0-ell_inf)*decay_convolution(t, a, k))
    return e, ell
