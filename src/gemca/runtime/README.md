# gemca/runtime — Design Notes

## What this layer does

`osh_gemca_runtime` compiles the internal compatibility `gemca_workspace`
(pointer-linked structs produced by the GEMCA prepare layer) into a flat,
cache-friendly representation that the transport kernel can query in a tight
inner loop.

Three contiguous arrays replace the cold layout:

| Array | Content |
|---|---|
| `surfaces[]` | All body surfaces concatenated, fixed-size parameter slots |
| `bodies[]` | Transform matrix, coordinate system tag, offset into `surfaces[]` |
| `zones[]` | Per-zone flat RPN instruction array, material index |

Zone membership is evaluated by a small integer stack machine that walks the
RPN array left-to-right.  The stack depth is bounded at compile time
(`OSH_GEMCA_RT_MAX_STACK`), so no dynamic allocation occurs in the hot path.

The batch API (`get_zone_batch`, `get_distance_batch`) accepts structure-of-arrays
(SoA) position and direction arguments so that SIMD and GPU implementations can
operate directly on contiguous memory without transposing.

---

## SIMD acceleration (AVX2)

`osh_gemca_runtime_avx2.c` is compiled with `-mavx2 -mfma` when the build
system detects support.  It is linked into the same library and dispatched at
runtime via `__builtin_cpu_supports("avx2")` in `get_zone_batch`.

The implementation processes **4 particles simultaneously** using `__m256d`
(4×64-bit double) vectors:

- Zone-outer, 4-particle-inner loop.
- `_in_body_avx2` — BZALIGN affine transform via 9 FMA instructions + surface
  checks, early exit via `_mm256_testz_si256`.
- `_eval_membership_avx2` — `__m256i` boolean stack; CSG operators become
  single `or`/`and`/`andnot` instructions.
- 4-bit scalar `resolved` bitmask tracks which lanes have found their zone;
  terminates the zone loop early when all 4 are resolved.
- Scalar tail (0–3 particles) falls back to `osh_gemca_runtime_get_zone`.

---

## GPU migration path (future)

The design was written with a future GPU port in mind.  Most of the structure
is already GPU-compatible; one gap remains.

### What is already GPU-ready

- **Flat contiguous arrays** (`surfaces[]`, `bodies[]`, `zones[]`) — three
  plain buffers, trivially copied to device memory (`cudaMemcpy` / `hipMemcpy`
  / SYCL USM).
- **SoA batch API** — `x[], y[], z[], ux[], uy[], uz[]` is exactly the layout
  for coalesced GPU reads.
- **Fixed-depth RPN stack** — `OSH_GEMCA_RT_MAX_STACK 32` slots fit in GPU
  registers; no dynamic allocation in the hot path.
- **Pure-math evaluators** — no I/O, no allocations, no callbacks.  The
  surface formulas and CSG operators are already written as stateless,
  side-effect-free computations — essentially shaders.
- **No function pointers** — dispatch is a `switch` on an integer type tag,
  which is warp-divergence rather than indirect calls (the acceptable form on
  GPU).

### The one structural gap: `insns[]` is a per-zone heap pointer

```c
struct gemca_rt_zone {
    struct gemca_rt_insn *insns;  /* heap pointer — not followable on GPU */
    int ninsns;
};
```

A GPU kernel cannot chase a host heap pointer.  The fix is a second flat layout
added alongside the current one — additive, no breakage of the CPU path:

```c
struct gemca_runtime {
    /* ... existing fields ... */

    /* GPU-portable flat instruction store (populated by setup_zones): */
    struct gemca_rt_insn *insns_flat;  /* all zones' instructions concatenated */
    int                  *insn_begin;  /* insn_begin[j] = offset into insns_flat for zone j */
};
```

`setup_zones` fills both.  The GPU kernel uses `insns_flat + insn_begin[j]`;
the CPU path keeps using `zones[j].insns`.  This is the only structural change
needed before a GPU kernel can be written.

### Minor host-only details to clean up before GPU

| Item | Location | Action for GPU |
|---|---|---|
| `osh_error()` in hot path | `eval_membership_batch_active` | Replace with per-particle error flag or silent no-op in GPU kernel |
| `struct ray` in scalar tail | `get_zone_batch_avx2` | Drop — no scalar tail on GPU; all lanes run in parallel |
| `__builtin_cpu_supports` dispatch | `get_zone_batch` | Host-only; kernel is launched unconditionally from host code |
| `workspace` pointer in `gemca_runtime` | diagnostics | Not copied to device |

### Recommended GPU kernel structure

```
__global__ void get_zone_batch_gpu(
    gemca_rt_surface const *surfaces,   /* device pointer */
    gemca_rt_body    const *bodies,     /* device pointer */
    gemca_rt_zone    const *zones,      /* device pointer (ninsns only) */
    gemca_rt_insn    const *insns_flat, /* device pointer */
    int              const *insn_begin, /* device pointer */
    double const *x, *y, *z, *ux, *uy, *uz,
    size_t n,
    size_t *zone_out)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    /* one thread per particle — evaluate RPN with register stack */
}
```

One thread per particle, the RPN stack lives in registers, and the flat
`insns_flat` array is streamed with L2/texture cache.  The existing
`_eval_membership_avx2` logic maps almost directly to the kernel body, with
`__m256d` replaced by scalar `double` (the GPU does its own SIMD across
threads in a warp).
