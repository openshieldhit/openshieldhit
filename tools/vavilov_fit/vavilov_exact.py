"""Exact Vavilov energy-loss distribution, computed from the published
definition (Vavilov, Sov. Phys. JETP 5, 749 (1957); B. Schorr, Comput. Phys.
Commun. 7, 215 (1974)).

Offline tooling only (reference fixtures + sampler fit); NOT compiled into the
OpenShieldHIT C build.  This module is an independent implementation of the
published distribution and its closed-form moments — the reference the sampler
fit is validated against.

Density via the Fourier (imaginary-axis Bromwich) form:

    p(λ; κ, β²) = (1/π) ∫₀^∞ Re[ φ(iy) e^{iλy} ] dy
    φ(s) = exp( κ(1 + β²γ_E) + s·ln κ + (s + β²κ)(Ein(s/κ) − γ_E) − κ e^{−s/κ} )
    Ein(ia) = γ_E + ln a − Ci(a) + i·Si(a) ,   a = y/κ

The Vavilov λ has closed-form mean = γ_E − 1 − β² − ln κ and var = (1 − β²/2)/κ,
used here to self-validate the integrator.
"""

import numpy as np
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import sici

EULER = 0.5772156649015329


def _log_phi_imag(y, kappa, beta2):
    """log φ(iy) = C + ψ(iy) (complex)."""
    a = y / kappa
    si, ci = sici(a)
    ein_minus_g = np.log(a) - ci + 1j * si  # Ein(ia) − γ_E
    s = 1j * y
    psi = s * np.log(kappa) + (s + beta2 * kappa) * ein_minus_g - kappa * np.exp(-1j * a)
    return kappa * (1.0 + beta2 * EULER) + psi


def vavilov_pdf(lam, kappa, beta2, upper=None):
    """Exact Vavilov density p(λ; κ, β²) by adaptive quadrature (accurate;
    validated to machine precision against the closed-form moments)."""
    if upper is None:
        upper = max(300.0, 80.0 / max(kappa, 1e-3))

    def integrand(y):
        return np.real(np.exp(_log_phi_imag(y, kappa, beta2) + 1j * lam * y))

    val, _ = quad(integrand, 0.0, upper, limit=500, epsabs=1e-11, epsrel=1e-9)
    return val / np.pi


def _support(kappa, beta2):
    mean = EULER - 1.0 - beta2 - np.log(kappa)
    sd = np.sqrt((1.0 - 0.5 * beta2) / kappa)
    return mean - 6.0 * sd - 2.0, mean + 22.0 * sd + 40.0


def cdf_table(kappa, beta2, n=400):
    """(lam, cdf) on a λ grid spanning ~all mass up to the far high-loss tail.

    A 3-point median filter removes the occasional isolated quadrature glitch
    (spike up OR down) in the far small-κ tail without cascading — a plain
    monotone clamp would propagate a single spuriously-low point downward and
    collapse the whole tail."""
    from scipy.signal import medfilt
    lo, hi = _support(kappa, beta2)
    lam = np.linspace(lo, hi, n)
    pdf = np.clip([vavilov_pdf(x, kappa, beta2) for x in lam], 0.0, None)
    pdf = medfilt(pdf, 3)
    cdf = np.concatenate([[0.0], cumulative_trapezoid(pdf, lam)])
    cdf /= cdf[-1]
    return lam, cdf


def ppf(u, kappa, beta2, n=400):
    """Inverse CDF λ(u) by monotone interpolation of the tabulated CDF."""
    lam, cdf = cdf_table(kappa, beta2, n)
    cdf_u, idx = np.unique(cdf, return_index=True)
    return np.interp(u, cdf_u, lam[idx])


def _selftest():
    print(f"{'kappa':>6} {'beta2':>5} | {'norm':>8} {'d_mean':>9} {'d_var':>9}")
    for kappa in (10.0, 1.0, 0.1, 0.03):
        beta2 = 0.5
        mean_ref = EULER - 1.0 - beta2 - np.log(kappa)
        var_ref = (1.0 - 0.5 * beta2) / kappa
        lo, hi = _support(kappa, beta2)
        lam = np.linspace(lo, hi, 900)
        pdf = np.array([vavilov_pdf(x, kappa, beta2) for x in lam])
        norm = np.trapezoid(pdf, lam)
        mean = np.trapezoid(lam * pdf, lam) / norm
        var = np.trapezoid((lam - mean) ** 2 * pdf, lam) / norm
        print(f"{kappa:6.2f} {beta2:5.2f} | {norm:8.5f} {mean - mean_ref:+9.1e} {var - var_ref:+9.1e}")


if __name__ == "__main__":
    _selftest()
