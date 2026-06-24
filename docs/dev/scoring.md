# Scoring in OpenShieldHIT

This page documents the scoring layer (`src/scoring/`): how cold detector and
filter definitions are compiled into runtime **accumulators**, how transport
deposits into them, how results are normalised and post-processed, and how the
accumulator / deposit / merge seams make history-level parallelism a small,
local change rather than a rewrite.

It complements the short orientation README next to the code
(`src/scoring/README.md`) — start there for the file map.

!!! info "Where the code lives"
    Cold API: `src/scoring/osh_scoring.{c,h}`.
    Runtime: `src/scoring/runtime/` — `osh_scoring_compile.*` (build),
    `osh_scoring_accumulator.*` (storage + deposit + merge),
    `osh_scoring_step.*` (hot path), `osh_scoring_postprocess.*`.
    Output writers: `src/scoring/save/`.

---

## 1. Pages, accumulators, and the deposit seam

Compilation (`osh_scoring_compile()`) turns each scored quantity into a
**page**. A page is split into two parts on purpose:

- The **descriptor** (`struct osh_scoring_page_runtime`) — geometry indices,
  bin strides, differential-axis configuration, scorer kind. Read-only during
  transport.
- The **accumulator** (`struct osh_scoring_accumulator`, embedded as
  `page->acc`) — the mutable, per-history storage.

The accumulator arrays are:

| Array | Meaning |
|---|---|
| `data` | Primary accumulator (always allocated). |
| `data2` | Secondary weight accumulator for two-pass averages (D/T-LET, D/T-Qeff); `NULL` for simple scorers. |
| `data_var`, `data2_var` | Running sum-of-squares for variance; `NULL` until the variance feature wires them. |
| `len` | Element count of every allocated array. |

Bins are row-major (`idx = ix + nx*(iy + ny*iz)`), with differential axes
appended as outer strides.

Every tally goes through one inline seam:

```c
static inline void osh_score_deposit(double *arr, size_t idx, double value) {
    arr[idx] += value;   /* the single place a write policy can change */
}
```

Today it inlines to a plain array store (zero overhead). Concentrating all
deposits here means a future parallel backend can swap in an atomic fetch-add,
a per-worker private store, or a locked update without touching any hot-path
call site.

## 2. Why storage is separate from the descriptor (parallel scoring)

Giving the per-history-written storage its own identity makes the accumulator
the unit a parallel worker can own **privately**: each worker scores its
assigned histories into a private accumulator set, then folds them into the
master with a single reduce:

```c
enum osh_status
osh_scoring_accumulator_merge(struct osh_scoring_accumulator *dst,
                              struct osh_scoring_accumulator const *src);
/* dst += src, element-wise over data / data2 / data_var / data2_var */
```

Because Monte Carlo histories are independent and the deposits commute, merging
in any order yields the same totals (up to floating-point summation order). The
same primitive serves threads (join then merge), MPI (`MPI_Reduce(MPI_SUM)`),
and — eventually — GPU block reductions. In the current single-worker build
each page simply owns one accumulator inline and no merge happens.

See [issue #161](https://github.com/openshieldhit/openshieldhit/issues/161) for
the parallelism roadmap these seams unblock.

## 3. Output formats: normalisation and multi-run merging

`osh_scoring_save()` dispatches on each output's `fileformat`. The two
general-purpose formats take deliberately different positions on normalisation;
DICOM RTDOSE is a specialised round-trip writer.

### ASCII (`text`, `txt`, `ascii`, `dat`)

For quick inspection and single-run use. Values are normalised **per primary
particle** (divided by `nstat`) at write time, except averaged quantities
(DLET, TLET) which are already physical means after post-processing and are
written as-is. ASCII output is **not** suitable for merging partial results:
once divided by `nstat`, each run's absolute weight is lost.

### BDO 2019 (`bdo`, `bdo2019`, `binary`, `bin`; default)

For production use and multi-run accumulation. Values are written **exactly as
accumulated** (raw sums) and the primary count is embedded as the
`OSHBDO_RT_NSTAT` tag; the reader is responsible for normalisation. This
supports two merging paths:

1. **Native multi-worker** (future): each worker accumulates into its own
   accumulator set, the simulation layer merges them (§2), and one BDO file is
   written with the total `nstat`.
2. **Embarrassingly-parallel / multi-node**: independent runs each emit a
   `.bdo`, and a merge tool weights per scorer kind using each file's `nstat`:
     - `NORM` (DOSE, FLUENCE, ENERGY, …): `X = (sum_j x_j) / (sum_j nstat_j)`
     - `AVER` (DLET, TLET): `X = (sum_j x_j * nstat_j) / (sum_j nstat_j)`
     - `SUM` (COUNT, …): `X = sum_j x_j`
     - `APPEND` (MCPL): concatenation

The `OSHBDO_PAG_NORMALIZE` tag records each page's postproc mode so a merge
tool need not re-derive it from the scorer type.

### DICOM RTDOSE (`rtdose`)

A round-trip writer for clinical interchange. It reads the RTDOSE template
recorded in `rtdose_template_path` (set by the app when the output format is
RTDOSE, and carried from compile through the geometry runtime to the save
layer), replaces the pixel data with scored values normalised by `nstat` and
`dose_grid_scaling`, and writes a new `.dcm`. All other DICOM metadata is
preserved unchanged. Only single-page (single-quantity) outputs are supported;
multi-page output returns `OSH_ENOTSUP`.

## 4. Unit handling and post-processing

`osh_scoring_postprocess()` applies physics transforms once per bin before
saving:

- **DOSE**: converts MeV/g → Gy (`× OSH_MEVG2GY`).
- **DLET, TLET, DQEFF, TQEFF**: finalises the two-pass average as `data / data2`
  per bin; the `data2` denominator is not written to any output file.

After post-processing, page buffers are in one of these states:

| Kind | State after postprocess | Save-layer normalisation |
|---|---|---|
| NORM (DOSE, ENERGY, FLUENCE, …) | raw accumulated sum | divide by `nstat` |
| AVER (DLET, TLET, DQEFF, TQEFF) | physical mean (`data ÷ data2` done) | none — written as-is |
| SUM (COUNT, …) | raw count | none |

**Merging caveat**: AVER pages cannot be naively summed across BDO files
because `data2` is discarded at post-process; merging requires re-weighting by
each file's `nstat`.

## 5. Lifecycle and module boundaries

```
osh_scoring_workspace_create()   allocate cold workspace
osh_scoring_compile()            allocate accumulators, compile scoring geometry
osh_scoring_postprocess()        unit conversion / two-pass finalisation
osh_scoring_save()               write output
osh_scoring_workspace_free()     release cold workspace
```

`osh_scoring_compile`, `osh_scoring_postprocess` and `osh_scoring_save` are
invoked by `src/simulation/` as part of `osh_simulation_create`,
`osh_simulation_run` and `osh_simulation_save`; app code never calls them
directly.

Boundaries:

- This directory owns cold scoring definitions, runtime compilation,
  deposition, post-processing and output support.
- File-format loading and path resolution happen in the app layer
  (`src/apps/osh/`), not here.
- Parsed scoring definitions stay separate from runtime scoring buffers.
- Save / output code consumes the scoring workspace plus scoring runtime, not
  transport internals.
