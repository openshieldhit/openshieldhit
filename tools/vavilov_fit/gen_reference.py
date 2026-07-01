"""Precompute a dense exact-Vavilov inverse-CDF reference table and cache it to
ref_cache.npz, so the coefficient fit does not re-run the (slow) quadrature.
Offline; the cache is a build artifact (git-ignored), regenerable from the exact
distribution in vavilov_exact.py."""

import os
import numpy as np
from vavilov_exact import ppf

HERE = os.path.dirname(__file__)
KAPPA = np.geomspace(0.01, 10.0, 22)
BETA2 = np.linspace(0.05, 0.99, 7)
# u grid: dense where λ varies fast (tails); capped at the vavinv valid range
# u ≤ 0.995 (above that the Vavilov high-loss tail is out of the sampler's range
# and κ < 0.01 is handled by the Landau sampler anyway).
U = np.unique(np.concatenate([
    np.geomspace(1e-3, 0.05, 10),
    np.linspace(0.05, 0.95, 31),
    1.0 - np.geomspace(5e-3, 0.05, 8)[::-1],  # 0.95..0.995
]))


def main():
    lam = np.empty((KAPPA.size, BETA2.size, U.size))
    for i, ka in enumerate(KAPPA):
        for j, be in enumerate(BETA2):
            lam[i, j] = ppf(U, float(ka), float(be), n=400)
        print(f"  kappa {ka:.4f} done ({i + 1}/{KAPPA.size})", flush=True)
    out = os.path.join(HERE, "ref_cache.npz")
    np.savez_compressed(out, kappa=KAPPA, beta2=BETA2, u=U, lam=lam)
    print(f"wrote {out}: lam shape {lam.shape}")


if __name__ == "__main__":
    main()
