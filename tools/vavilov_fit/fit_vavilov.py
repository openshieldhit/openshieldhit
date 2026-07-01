"""Fit our own rational + Chebyshev inverse-CDF coefficients to the exact
Vavilov reference (cached by gen_reference.py).  Clean-room: coefficients are
determined solely by the exact distribution; the runtime C evaluator reuses the
vavinv *architecture* (region branches; per region λ = P(x')/Q(x') with each
coefficient a 2-D Chebyshev sum in (κ, β)) but none of Thomsen's numbers.

Per region: linearised rational least squares with IRLS reweighting (Loeb):
minimise Σ w·(λ·Q − P)², w = 1/Q_prev², Q normalised (constant term 1).  x' is
normalised to ~[-1,1] within each band to condition the least squares.
"""

import os
import numpy as np

HERE = os.path.dirname(__file__)
NKA = 4  # Chebyshev order in kappa (T0..T3), variable = normalised ln κ
NBE = 3  # Chebyshev order in beta  (T0..T2)
NT = NKA * NBE


def cheb(kn, bn):
    tk = [1.0, kn, 2 * kn * kn - 1, (4 * kn * kn - 3) * kn]
    tb = [1.0, bn, 2 * bn * bn - 1]
    return np.array([tk[i] * tb[j] for i in range(NKA) for j in range(NBE)])


def xform(u, kind):
    return {"nlog": -np.log(u), "u": u, "1mu": 1.0 - u,
            "inv1mu": 1.0 / (1.0 - u),
            "nlog1mu": -np.log1p(-u)}[kind]  # −log(1−u): linearises the Vavilov high-loss tail


def _norm_lnk(k, klo, khi):
    return (2 * np.log(k) - np.log(klo) - np.log(khi)) / (np.log(khi) - np.log(klo))


def fit_region(cache, u_lo, u_hi, kind, ki, blo, bhi, mdeg, ndeg, iters=8):
    K, B, U, LAM = cache["kappa"], cache["beta2"], cache["u"], cache["lam"]
    klo, khi = K[ki[0]], K[ki[-1]]
    ui = np.where((U >= u_lo - 1e-12) & (U <= u_hi + 1e-12))[0]
    xall = xform(U[ui], kind)
    xmid = 0.5 * (xall.min() + xall.max())
    xhalf = 0.5 * (xall.max() - xall.min()) or 1.0

    def xn(u):
        return (xform(u, kind) - xmid) / xhalf

    kn_l, bn_l, x_l, lam = [], [], [], []
    for i in ki:
        kn = _norm_lnk(K[i], klo, khi) if khi > klo else 0.0
        for j in range(B.size):
            bn = (2 * B[j] - blo - bhi) / (bhi - blo)
            for m in ui:
                kn_l.append(kn); bn_l.append(bn)
                x_l.append(xn(U[m])); lam.append(LAM[i, j, m])
    lam = np.array(lam); x_l = np.array(x_l)
    T = np.array([cheb(a, b) for a, b in zip(kn_l, bn_l)])

    def solve(qprev):
        cols = [(x_l ** k)[:, None] * T for k in range(mdeg + 1)]
        cols += [-(lam * x_l ** k)[:, None] * T for k in range(1, ndeg + 1)]
        A = np.concatenate(cols, axis=1)
        w = 1.0 / np.maximum(np.abs(qprev), 1e-3)
        c, *_ = np.linalg.lstsq(A * w[:, None], lam * w, rcond=None)
        na = (mdeg + 1) * NT
        ac = c[:na].reshape(mdeg + 1, NT); bc = c[na:].reshape(ndeg, NT)
        Q = np.ones_like(lam)
        for k in range(1, ndeg + 1):
            Q = Q + (x_l ** k) * (T @ bc[k - 1])
        return ac, bc, Q

    q = np.ones_like(lam); ac = bc = None
    for _ in range(iters):
        ac, bc, q = solve(q)

    md = 0.0
    for i in ki:
        kn = _norm_lnk(K[i], klo, khi) if khi > klo else 0.0
        for j in range(B.size):
            bn = (2 * B[j] - blo - bhi) / (bhi - blo)
            t = cheb(kn, bn); xx = xn(U[ui])
            P = sum((xx ** k) * (t @ ac[k]) for k in range(mdeg + 1))
            Q = 1.0 + sum((xx ** k) * (t @ bc[k - 1]) for k in range(1, ndeg + 1))
            md = max(md, float(np.max(np.abs(P / Q - LAM[i, j, ui]))))
    return ac, bc, md, xmid, xhalf


def kbands_by_index(nk, nbands, overlap=1):
    """Contiguous κ-index bands (>= NKA nodes each) with 1-node overlap."""
    edges = np.linspace(0, nk, nbands + 1).astype(int)
    bands = []
    for b in range(nbands):
        lo = max(0, edges[b] - (overlap if b else 0))
        hi = min(nk, edges[b + 1] + (overlap if b < nbands - 1 else 0))
        bands.append(np.arange(lo, hi))
    return bands


if __name__ == "__main__":
    cache = np.load(os.path.join(HERE, "ref_cache.npz"))
    K, U = cache["kappa"], cache["u"]
    ubands = [(U.min(), 0.05, "nlog"), (0.05, 0.6, "u"),
              (0.6, 0.95, "1mu"), (0.95, U.max(), "inv1mu")]
    degs = {"nlog": (4, 4), "u": (3, 3), "1mu": (4, 3), "inv1mu": (5, 5)}
    for ki in kbands_by_index(K.size, 3):
        for (ulo, uhi, kind) in ubands:
            m, n = degs[kind]
            ac, bc, md, _, _ = fit_region(cache, ulo, uhi, kind, ki, 0.05, 0.99, m, n)
            print(f"k[{K[ki[0]]:.3g},{K[ki[-1]]:.3g}] nk={ki.size:2d} "
                  f"u[{ulo:.3g},{uhi:.3g}]({kind:6s}) deg({m},{n}): maxdev={md:.3e}")
