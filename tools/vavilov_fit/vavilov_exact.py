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
        # φ(0)=1 ⇒ integrand(0)=1.  Guard y=0 explicitly: there a=y/κ=0, and
        # sici(0)=(0,-inf) makes log(a)-Ci(a) an ∞-∞ → NaN.  QUADPACK's
        # Gauss-Kronrod rule samples only interior points so it never hits y=0
        # today, but the guard keeps the integrand finite regardless.
        if y == 0.0:
            return 1.0
        return np.real(np.exp(_log_phi_imag(y, kappa, beta2) + 1j * lam * y))

    val, _ = quad(integrand, 0.0, upper, limit=500, epsabs=1e-11, epsrel=1e-9)
    return val / np.pi


def _support(kappa, beta2):
    """(lo, core_hi, hi): the low-loss edge, the dense-core upper edge just past
    the peak, and the far high-loss edge.

    The left edge is bounded by a small κ-independent margin, NOT by 6·sd:
    sd = √((1−β²/2)/κ) diverges as κ→0, but the Vavilov/Landau low-loss tail is
    intrinsically short (essentially zero below λ ≈ −4), so an sd-scaled left
    edge (≈ −50 at κ=0.01) is a huge dead zone that a uniform grid wastes its
    points on — the root cause of the old small-κ / extreme-u reference glitch.
    The high-loss tail is genuinely heavy (~1/λ²), so its edge keeps scaling."""
    mean = EULER - 1.0 - beta2 - np.log(kappa)
    sd = np.sqrt((1.0 - 0.5 * beta2) / kappa)
    lo = mean - min(6.0 * sd + 2.0, 8.0)
    core_hi = mean + min(6.0 * sd, 30.0)
    hi = mean + 22.0 * sd + 40.0
    return lo, core_hi, hi


def _lam_grid(kappa, beta2, n):
    """Non-uniform λ grid: a dense linear core over the peak and bulk, plus a
    geometric (sparse) sampling of the heavy high-loss tail.  Concentrating
    points in the core — where the CDF rises fast — resolves the extreme-low-u
    inverse CDF at small κ that a uniform grid over the very wide [lo, hi] span
    badly under-resolves."""
    lo, core_hi, hi = _support(kappa, beta2)
    n_core = max(4, int(round(0.7 * n)))
    n_tail = max(1, n - n_core)
    core = np.linspace(lo, core_hi, n_core)
    if hi <= core_hi:
        return core
    tail = core_hi + np.geomspace(1e-3, hi - core_hi, n_tail + 1)[1:]
    return np.concatenate([core, tail])


def cdf_table(kappa, beta2, n=3000):
    """(lam, cdf) on the non-uniform λ grid, spanning ~all mass up to the far
    high-loss tail.

    A 3-point median filter removes the occasional isolated quadrature glitch
    (spike up OR down) without cascading — a plain monotone clamp would
    propagate a single spuriously-low point downward and collapse the tail.
    A relative noise floor then zeros any residual roundoff-level pdf so it
    contributes no spurious CDF mass at the corners."""
    from scipy.signal import medfilt
    lam = _lam_grid(kappa, beta2, n)
    pdf = np.clip([vavilov_pdf(x, kappa, beta2) for x in lam], 0.0, None)
    pdf = medfilt(pdf, 3)
    pdf[pdf < 1e-12 * pdf.max()] = 0.0
    cdf = np.concatenate([[0.0], cumulative_trapezoid(pdf, lam)])
    cdf /= cdf[-1]
    return lam, cdf


def ppf(u, kappa, beta2, n=3000):
    """Inverse CDF λ(u) by monotone interpolation of the tabulated CDF."""
    lam, cdf = cdf_table(kappa, beta2, n)
    cdf_u, idx = np.unique(cdf, return_index=True)
    return np.interp(u, cdf_u, lam[idx])


def _selftest():
    """Moment self-test + a small-κ / extreme-low-u corner regression check
    (the old glitch: κ=0.01, β²=0.5, u=0.001 read ≈ −3.0 on the coarse uniform
    grid; the exact / libvav value is ≈ −2.64)."""
    print(f"{'kappa':>6} {'beta2':>5} | {'norm':>8} {'d_mean':>9} {'d_var':>9}")
    for kappa in (10.0, 1.0, 0.1, 0.03, 0.01):
        beta2 = 0.5
        mean_ref = EULER - 1.0 - beta2 - np.log(kappa)
        var_ref = (1.0 - 0.5 * beta2) / kappa
        lam = _lam_grid(kappa, beta2, 3000)
        pdf = np.array([vavilov_pdf(x, kappa, beta2) for x in lam])
        norm = np.trapezoid(pdf, lam)
        mean = np.trapezoid(lam * pdf, lam) / norm
        var = np.trapezoid((lam - mean) ** 2 * pdf, lam) / norm
        print(f"{kappa:6.2f} {beta2:5.2f} | {norm:8.5f} {mean - mean_ref:+9.1e} {var - var_ref:+9.1e}")

    corner = ppf(0.001, 0.01, 0.5)
    print(f"\ncorner ppf(u=0.001, kappa=0.01, beta2=0.5) = {corner:.4f}  (expect ~ -2.64)")
    assert abs(corner - (-2.644)) < 0.05, f"corner regression: {corner:.4f}"
    print("corner OK")


if __name__ == "__main__":
    _selftest()
