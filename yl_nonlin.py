#!/usr/bin/env python3
import math
import numpy as np


def _compute_gradients(u, dr, dtheta):
    u_r = np.zeros_like(u)
    u_theta = np.zeros_like(u)

    u_r[1:-1, :] = (u[2:, :] - u[:-2, :]) / (2.0 * dr)
    u_r[0, :] = (u[1, :] - u[0, :]) / dr
    u_r[-1, :] = (u[-1, :] - u[-2, :]) / dr

    u_theta[:, :] = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2.0 * dtheta)
    return u_r, u_theta


def _compute_c(u_r, u_theta, r):
    r_safe = r.copy()
    r_safe[0] = 1.0
    grad_sq = u_r * u_r + (u_theta / r_safe[:, None]) ** 2
    c = 1.0 / np.sqrt(1.0 + grad_sq)
    c[0, :] = 1.0 / np.sqrt(1.0 + u_r[0, :] * u_r[0, :])
    return c


def _surface_area(u, r, dr, dtheta):
    u_r, u_theta = _compute_gradients(u, dr, dtheta)
    r_safe = r.copy()
    r_safe[0] = 1.0
    integrand = np.sqrt(1.0 + u_r * u_r + (u_theta / r_safe[:, None]) ** 2)
    integrand[0, :] = np.sqrt(1.0 + u_r[0, :] * u_r[0, :])
    return float(np.sum(integrand * r[:, None]) * dr * dtheta)


def yl_nonlin(
    rho,
    sigma,
    g,
    omega_rpm,
    a_orbit,
    R,
    thetac_deg,
    hmax=None,
    n_theta=None,
    tol=1e-6,
    max_iter=2000,
    omega_relax=1.2,
):
    """
    Solve the nonlinear Young-Laplace equation in a circular domain using
    a polar-grid finite difference method.

    Returns:
        A: surface area of the interface
        hL: height difference between x=-R and x=+R along y=0
        pts: array of [x, y, z] points on the interface
    """
    lc = math.sqrt(sigma / (rho * g))
    omega = omega_rpm * 2.0 * math.pi / 60.0
    F = a_orbit * omega * omega / g
    thetac = thetac_deg * math.pi / 180.0
    g_bc = math.tan(math.pi / 2.0 - thetac)
    a_coeff = 1.0 / (lc * lc)

    if hmax is None:
        hmax = 0.02 * R
    dr = float(hmax)
    nr = int(math.ceil(R / dr)) + 1
    dr = R / (nr - 1)

    if n_theta is None:
        n_theta = int(math.ceil(2.0 * math.pi * R / dr))
        n_theta = max(n_theta, 32)
    if n_theta % 2 == 1:
        n_theta += 1
    dtheta = 2.0 * math.pi / n_theta

    r = np.linspace(0.0, R, nr)
    theta = np.linspace(0.0, 2.0 * math.pi, n_theta, endpoint=False)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    u = np.zeros((nr, n_theta), dtype=float)
    f = -F * a_coeff * (r[:, None] * cos_t[None, :])

    for _ in range(max_iter):
        u[0, :] = u[1, :]
        u_r, u_theta = _compute_gradients(u, dr, dtheta)
        c = _compute_c(u_r, u_theta, r)
        u[-1, :] = u[-2, :] + dr * g_bc / np.clip(c[-1, :], 1e-10, None)

        u_r, u_theta = _compute_gradients(u, dr, dtheta)
        c = _compute_c(u_r, u_theta, r)

        max_delta = 0.0
        for i in range(1, nr - 1):
            r_i = r[i]
            r_e = r_i + dr / 2.0
            r_w = r_i - dr / 2.0
            for j in range(n_theta):
                jp = (j + 1) % n_theta
                jm = (j - 1) % n_theta
                c_p = c[i, j]
                c_e = 0.5 * (c_p + c[i + 1, j])
                c_w = 0.5 * (c_p + c[i - 1, j])
                c_n = 0.5 * (c_p + c[i, jp])
                c_s = 0.5 * (c_p + c[i, jm])

                a_r = (r_e * c_e) / (r_i * dr * dr)
                a_l = (r_w * c_w) / (r_i * dr * dr)
                a_t = c_n / (r_i * r_i * dtheta * dtheta)
                a_b = c_s / (r_i * r_i * dtheta * dtheta)
                diag = a_r + a_l + a_t + a_b + a_coeff
                rhs = (
                    f[i, j]
                    + a_r * u[i + 1, j]
                    + a_l * u[i - 1, j]
                    + a_t * u[i, jp]
                    + a_b * u[i, jm]
                )
                u_new = rhs / diag
                updated = (1.0 - omega_relax) * u[i, j] + omega_relax * u_new
                delta = abs(updated - u[i, j])
                if delta > max_delta:
                    max_delta = delta
                u[i, j] = updated

        u[0, :] = u[1, :]
        u[-1, :] = u[-2, :] + dr * g_bc / np.clip(c[-1, :], 1e-10, None)
        if max_delta < tol:
            break

    area = _surface_area(u, r, dr, dtheta)

    r_sample = max(R - 0.5 * dr, r[0])
    i1 = int(np.searchsorted(r, r_sample))
    i0 = max(i1 - 1, 0)
    if i1 == i0:
        weight = 0.0
    else:
        weight = (r_sample - r[i0]) / (r[i1] - r[i0])

    j0 = 0
    jpi = n_theta // 2
    u_right = (1.0 - weight) * u[i0, j0] + weight * u[i1, j0]
    u_left = (1.0 - weight) * u[i0, jpi] + weight * u[i1, jpi]
    hL = float(u_left - u_right)

    x = r[:, None] * cos_t[None, :]
    y = r[:, None] * sin_t[None, :]
    pts = np.column_stack((x.ravel(), y.ravel(), u.ravel()))

    return area, hL, pts


__all__ = ["yl_nonlin"]
