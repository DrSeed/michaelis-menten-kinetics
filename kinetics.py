# Michaelis-Menten parameter estimation: two linearisations and a proper nonlinear fit.
import numpy as np


def mm_rate(s, vmax, km):
    return vmax * s / (km + s)


def fit_lineweaver_burk(s, v):
    # 1/v = (Km/Vmax)(1/S) + 1/Vmax
    slope, intercept = np.polyfit(1 / np.asarray(s, float), 1 / np.asarray(v, float), 1)
    vmax = 1 / intercept
    return vmax, slope * vmax


def fit_eadie_hofstee(s, v):
    # v = Vmax - Km (v/S)
    s, v = np.asarray(s, float), np.asarray(v, float)
    slope, intercept = np.polyfit(v / s, v, 1)
    return intercept, -slope


def fit_nonlinear(s, v, p0=None, iters=200):
    """Levenberg-Marquardt least squares; returns (vmax, km, se_vmax, se_km)."""
    s, v = np.asarray(s, float), np.asarray(v, float)
    vmax, km = p0 if p0 is not None else (v.max(), float(np.median(s)))
    lam = 1e-3

    def sse(a, b):
        return float(np.sum((v - mm_rate(s, a, b)) ** 2))

    def jac(a, b):
        return np.column_stack([s / (b + s), -a * s / (b + s) ** 2])

    cur = sse(vmax, km)
    for _ in range(iters):
        J = jac(vmax, km)
        r = v - mm_rate(s, vmax, km)
        A = J.T @ J
        step = np.linalg.solve(A + lam * np.diag(np.diag(A)), J.T @ r)
        new_v, new_k = vmax + step[0], km + step[1]
        if new_k > 0 and sse(new_v, new_k) < cur:
            vmax, km, cur = new_v, new_k, sse(new_v, new_k)
            lam *= 0.3
        else:
            lam *= 10
        if np.max(np.abs(step)) < 1e-10 or lam > 1e10:
            break
    J = jac(vmax, km)
    cov = cur / max(len(s) - 2, 1) * np.linalg.inv(J.T @ J)
    return vmax, km, float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
