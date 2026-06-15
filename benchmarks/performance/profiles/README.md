# Committed callgrind flat profiles

Reference instruction-count profiles for representative scenarios, produced by

```bash
python3 benchmarks/performance/run_bench.py \
    --filter c1_p100_dose,c4_p150_nucre_filters,c6_c12_let,c8_z1000 \
    --nstat-scale 0.2 --repeats 1 --callgrind --prof-nstat-divisor 8
```

callgrind counts instructions deterministically, so unlike wall-clock numbers
these tables are reproducible and largely machine-independent — that is why
they are committed while wall-clock results are not. They were generated at
reduced nstat (callgrind is ~30x slower than native); instruction *shares*
converge after a few primaries, absolute counts scale with nstat.

Headline shares (Ir %, gcc 13 -O3, x86-64):

| Symbol / bucket | C1 (clean CSDA) | C4 (NUCRE+filters) | C6 (C-12 LET) | C8 ~1000 zones |
|---|---|---|---|---|
| `eval_distance` | 18.8 | 15.0 | 13.6 | **70.9** |
| `osh_ray_transform` | 12.4 | 9.9 | 9.0 | — |
| `osh_gemca_runtime_get_zone` | 2.1 | 1.6 | 1.5 | **17.6** |
| libm transcendentals (log/pow/exp/cbrt/…) | ~27 | ~27 | ~26 | — |
| `osh_transport_ion_step` | 10.8 | 8.8 | 8.0 | — |
| `osh_scoring_score_step` | 5.7 | 5.5 | 9.2 | — |
| `osh_raytrace_traverse` | 3.9 | 3.3 | 3.7 | — |
| RNG (`osh_rng_*`) | ~4.0 | ~3.5 | ~2.5 | — |

Regenerate after any transport/geometry/physics change that aims at these
buckets, and compare with `benchmarks/performance/compare.py` on the wall-clock side.
