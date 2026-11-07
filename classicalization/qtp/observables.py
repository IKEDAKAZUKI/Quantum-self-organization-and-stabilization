"""Physical occupation fields and Fourier observables used in the article."""
from __future__ import annotations
from typing import Mapping
import numpy as np


def occupation_contrast(n: object, s: object, eps: float) -> np.ndarray:
    """Decode the centered integer n and leakage label s, u=eps*(n-I[s=0]/2)."""
    n, s = np.asarray(n), np.asarray(s)
    if n.shape != s.shape or n.ndim < 2:
        raise ValueError("n and s must have identical shape (..., ny, nx).")
    if not np.all(np.isin(s, [-1, 0, 1])):
        raise ValueError("Saved state labels must be -1, 0, or 1; 9 is only an initialization instruction.")
    if not np.isfinite(eps) or eps <= 0:
        raise ValueError("eps must be finite and positive.")
    return float(eps)*(n.astype(float)-0.5*(s == 0))


def mean_coherence(s: object, mode: int) -> np.ndarray:
    """Quantum mode 0: spatial mean of s. Projected mode 1: identically zero.

    In mode 1, stored s=1 marks a nonleaked cell, not relative-phase coherence.
    """
    s = np.asarray(s)
    if s.ndim < 2 or not np.all(np.isin(s, [-1, 0, 1])):
        raise ValueError("s must be a saved label array (..., ny, nx).")
    if mode == 0:
        return s.mean(axis=(-2, -1))
    if mode == 1:
        return np.zeros(s.shape[:-2], dtype=float)
    raise ValueError("mode must be 0 (quantum) or 1 (number projected).")


def fixed_band(p: Mapping) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Wavevectors, Laplacian eigenvalues and the fixed unstable band at e=1."""
    nx, ny = int(p['nx']), int(p['ny'])
    kx = 2*np.pi*np.fft.fftfreq(nx)[None, :]
    ky = 2*np.pi*np.fft.fftfreq(ny)[:, None]
    wx, wy = 2-2*np.cos(kx), 2-2*np.cos(ky)
    omega = wx+wy
    sigma = -p['delta']+p['D2']*omega-p['D4']*omega*omega-p['zeta']*wy+p['Dcl']*omega
    band = sigma > 0
    band[0, 0] = False
    return kx, ky, omega, band


def spectrum(u: object, p: Mapping) -> dict[str, np.ndarray]:
    """Band power, mean squared contrast and fixed-band orientation statistics.

    The band is held fixed when dephasing or source rates are changed. ``energy``
    denotes mean squared occupation contrast, not Hamiltonian energy.
    """
    u = np.asarray(u, dtype=float)
    if u.ndim < 2 or u.shape[-2:] != (int(p['ny']), int(p['nx'])):
        raise ValueError("Field shape must end in (ny, nx) from the parameters.")
    if not np.all(np.isfinite(u)):
        raise ValueError("The field contains a nonfinite value.")
    kx, ky, omega, band = fixed_band(p)
    if not np.any(band):
        raise ValueError("The fixed e=1 unstable band is empty for these parameters.")
    power_grid = abs(np.fft.fft2(u, axes=(-2, -1))/(p['nx']*p['ny']))**2
    power = power_grid[..., band].sum(axis=-1)
    energy = np.mean(u*u, axis=(-2, -1))
    angles = np.arctan2(ky, kx)
    nematic = np.abs((power_grid[..., band]*np.exp(2j*angles[band])).sum(axis=-1))/np.maximum(power, 1e-30)
    selected = np.where(band, power_grid, -1)
    shape = u.shape[:-2]
    flat = selected.reshape(-1, p['nx']*p['ny'])
    index = flat.argmax(axis=-1)
    return dict(power=power, energy=energy, nematic=nematic,
                k=np.sqrt((kx*kx+ky*ky).ravel()[index]).reshape(shape),
                omega=omega.ravel()[index].reshape(shape), peak=flat.max(axis=-1).reshape(shape))
