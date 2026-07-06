# Scoring in OpenShieldHIT

This page documents the scoring layer (`src/scoring/`): the structures that hold
tallied results, how transport writes into them, how per-worker results combine
for parallelism, and how results are normalised, post-processed, and saved.

It complements the short orientation README next to the code
(`src/scoring/README.md`) — start there for the file map.

!!! info "Where the code lives"
    Cold API: `src/scoring/osh_scoring.{c,h}`.
    Runtime: `src/scoring/runtime/` — `osh_scoring_compile.*` (build),
    `osh_scoring_accumulator.*` (storage + deposit + merge),
    `osh_scoring_step.*` (hot path), `osh_scoring_postprocess.*`.
    Output writers: `src/scoring/save/`.

---

## 1. Overview: the data hierarchy, and "cold" vs "hot"

A single OpenShieldHIT run produces one or more **outputs**, each written to its
own file. Within an output are one or more **pages**:

```
run ─┬─ output (file) ─┬─ page  (quantity + particle filters, over a shared geometry)
     │                 ├─ page
     │                 └─ …
     ├─ output (file) ─── page
     └─ …
```

All pages of an output share the same **scorer geometry** — the binning grid,
e.g. a 100×1×1 mesh along depth, or a cylindrical R×Z grid. A page then picks
*which physical quantity* it tallies (dose, fluence, LET, …) and *which particles
count* (filters on charge, generation, …). Some pages also carry a
**differential axis** — for example dose *as a function of energy* — which
multiplies the bin count by the number of differential bins; such a page has a
higher dimensionality than its plain siblings (one or two differential axes are
supported).

The code runs in two phases, at two "temperatures":

- **Cold** code runs once at setup. It parses the user's detector and filter
  definitions and *compiles* them into runtime structures. Clarity matters here,
  speed does not.
- **Hot** code runs inside the transport loop, once per particle step — millions
  to billions of times per run. This is where scoring writes results. Every
  instruction counts.

Three terms used throughout:

- A **tally** is the running total a bin holds: the sum of every contribution
  written into it so far.
- A **deposit** is one such contribution (e.g. the energy a single step leaves in
  a voxel). Depositing is the hot-path write.
- An **accumulator** is the block of storage holding one page's tallies.

`osh_scoring_compile()` is the cold step that turns each scored quantity into a
page. The rest of this document follows a page from compilation, through
hot-path deposits, optional parallel merge, post-processing, and save.

## 2. Pages and accumulators

A page is split into two parts on purpose:

- The **descriptor** (`struct osh_scoring_page_runtime`) — geometry indices, bin
  strides, differential-axis configuration, scorer kind. Read-only during
  transport (cold data consumed on the hot path).
- The **accumulator** (`struct osh_scoring_accumulator`, embedded as `page->acc`)
  — the mutable, per-history storage written by deposits.

The accumulator arrays are:

| Array | Meaning |
|---|---|
| `data` | Primary accumulator (always allocated): running sum of per-history deposits Σx per bin. |
| `data2` | Secondary weight accumulator for two-pass averages (D/T-LET, D/T-Qeff); `NULL` for simple scorers. |
| `data_var`, `data2_var` | Welford **M2** (sum of squared deviations from the mean) for `data` / `data2`; `NULL` until the variance feature wires them. See §4. |
| `len` | Element count of every allocated array. |
| `weight`, `nbatch` | Per-accumulator scalars: statistical weight `W` (history count) and number of batches `B` folded in. See §4. |

Bins are row-major (`idx = ix + nx*(iy + ny*iz)`), with differential axes
appended as outer strides.

Every deposit funnels through one tiny inline function. We call this a **seam**:
a single point all writes pass through, so the write *policy* can be changed in
one place instead of at every call site.

```c
static inline void osh_score_deposit(double *arr, size_t idx, double value) {
    arr[idx] += value;   /* the single place a write policy can change */
}
```

Today it inlines to a plain array store (zero overhead). Concentrating all
deposits here means a future parallel backend can swap in an atomic fetch-add, a
per-worker private store, or a locked update without touching any hot-path call
site.

## 3. Private accumulators and merging (parallel scoring)

Giving the per-history-written storage its own identity (separate from the
read-only descriptor) makes the accumulator the unit a parallel worker can own
**privately**: each worker scores its assigned histories into a private
accumulator set, then folds them into the master with a single reduce:

```c
enum osh_status
osh_scoring_accumulator_merge(struct osh_scoring_accumulator *dst,
                              struct osh_scoring_accumulator const *src);
/* dst ⊕ src: sums add element-wise; variance state combines (§4). */
```

Because Monte Carlo histories are independent and the deposits commute, merging
in any order yields the same totals (up to floating-point summation order). The
same primitive serves threads (join then merge), MPI (`MPI_Reduce` on the raw
sums), and — eventually — GPU block reductions. In the current single-worker
build each page simply owns one accumulator inline and no merge happens.

See [issue #161](https://github.com/openshieldhit/openshieldhit/issues/161) for
the parallelism roadmap this descriptor/storage split unblocks.

## 4. Statistical uncertainty

A Monte Carlo result needs an error bar. The representation chosen for it is what
the merge in §3 must compute, so the two are designed together (decision recorded
in [issue #169](https://github.com/openshieldhit/openshieldhit/issues/169)).

### Background: online (running) variance

Computing a variance the textbook way needs every sample in memory: take the
mean, then sum the squared deviations. **Welford's online algorithm**
(Welford, 1962) avoids storing samples: it keeps only a running count, a running
mean, and **M2** — the running *sum of squared deviations from the current mean*
— and updates all three with each new sample in a single pass. The variance is
then `M2 / (n − 1)`. So whenever this document says **M2**, it means that one
running quantity from Welford's method.

Two extensions make it fit a parallel Monte Carlo:

- **West (1979)** generalised the update to *weighted* samples — each sample
  carries a weight, not just an implicit count of one.
- **Schubert & Gertz (2018)** gave the numerically stable rule for **combining
  two independently-accumulated partials** — two `(weight, mean, M2)` results —
  into one. That is exactly a parallel reduction.

OpenShieldHIT stores `M2` per bin (in `data_var`) and combines partials with the
Schubert–Gertz formula below.

### Why batch-means (not per-history)

The textbook MC uncertainty is the variance of the per-primary mean over
**independent histories**. Computing it exactly needs each history's complete
per-bin tally before squaring — an `nbins × live-histories` scratch buffer, which
is prohibitive for the wavefront/pool transport (many histories are in flight at
once, depositing per step). So instead a **batch** — one independent unit of work
— is treated as a single (weighted) observation, and the spread *between batches*
estimates the error. A batch is any of:

- a parallel worker's history range (threads / MPI ranks — still to come),
- a periodic partial-dump / checkpoint interval — the batch loop folds one
  observation per checkpoint (see the run-control work, #170),
- an internal sub-split of a single serial run — a plain `VARIANCE` run derives a
  count cadence that yields `N` batches out of the box (issue #209),
- a `--score-replicas N` range, or a separate independent run.

Because every bin is exposed to the **same** set of histories (most depositing
zero), the observation count is identical across bins — so the batch weight `W`
and batch count `B` are **per-accumulator scalars** (`weight`, `nbatch`), not
per-bin arrays. The only per-bin variance state is one extra array, the Welford
`M2` (`data_var`).

### The merge contract

Per bin the accumulator stores the running sum `data` (= Σx) and the Welford
`M2` (`data_var`); the per-primary mean is `data / weight`. Two batches A and B
combine with the Schubert–Gertz formula (the single-batch fold is West's weighted
update with one observation):

```
w  = wA + wB
δ  = meanB − meanA            (means = data / weight)
M2 = M2_A + M2_B + δ² · wA·wB / w      ← the cross-term a plain += cannot express
data, data2, weight, nbatch  → additive
```

`osh_scoring_accumulator_merge()` applies exactly this: the raw sums and the
batch bookkeeping are additive, but the `M2` arrays use the cross-term. A blanket
`+=` over `data_var` would silently drop `δ²·wA·wB/w` and corrupt the variance —
the reason the merge is representation-aware rather than a flat element-wise add.
The `wA·wB/w` weighting handles **unequal-size batches by construction**, which is
the normal case under heterogeneous CPU cores and arbitrary dump boundaries.

### Degrees of freedom and finalisation

At save time the standard error of the per-primary mean in a bin is

```
SE = sqrt( M2 / ((nbatch − 1) · weight) )
```

so a single batch (`nbatch == 1`, e.g. a plain serial run with no sub-splitting)
has **zero degrees of freedom and no error estimate** — at least two batches are
required.

`osh_scoring_finalize_errors()` (issue #209) applies this once, **before**
`osh_scoring_postprocess()` rescales `data` or collapses the two-pass ratios,
because it reads the raw sums. Rather than the absolute `SE` it stores the per-bin
**relative** error `SE / |mean| = sqrt(M2 · weight / (nbatch−1)) / |data|` back
into `data_var`; the save layer then emits the absolute error column as
`|value| · data_var`, so whatever per-primary / physical-mean / unit scaling the
writer applies to the value applies to its error automatically — the stored error
is normalisation-invariant. For the two-pass **AVER** quantities (DLET/TLET/Qeff)
the reported value is the ratio `data / data2`, so the numerator and denominator
relative errors are combined in quadrature (`rel² = rel_num² + rel_den²`); this
ignores their strong positive correlation and is therefore deliberately
**conservative** — an over-estimate, never an under-estimate. A bin with
`nbatch < 2`, a zero mean, or zero weight finalises to a `0` error.

The feature is off by default and enabled per run by the global `VARIANCE [N]`
card in `detect.dat` (see [`detect.dat.md`](../user/detect.dat.md)), which
allocates the companion `M2` (`data_var`) arrays via
`osh_scoring_accumulator_alloc_variance()`. With it off, the `M2` arrays are never
allocated and the accumulators and hot path are byte-for-byte unchanged.

> **MPI / GPU note.** The additive fields can ride `MPI_Reduce(MPI_SUM)`, but the
> `M2` arrays cannot — a rank reduction needs a custom `MPI_Op_create` that
> applies the combine above to the whole accumulator (sums + weight + M2
> together). On GPU, per-worker private accumulators reduced on the host with this
> same merge avoid needing float atomics entirely.

## 5. Output formats: normalisation and multi-run merging

`osh_scoring_save()` dispatches on each output's `fileformat`. The two
general-purpose formats take deliberately different positions on normalisation;
DICOM RTDOSE is a specialised round-trip writer.

### ASCII (`text`, `txt`, `ascii`, `dat`)

For quick inspection and single-run use. Values are normalised **per primary
particle** (divided by `nstat`) at write time, except averaged quantities
(DLET, TLET) which are already physical means after post-processing and are
written as-is. ASCII output is **not** suitable for merging partial results:
once divided by `nstat`, each run's absolute weight is lost.

When the `VARIANCE` card is active each quantity gains a paired `NAME_ERR` column
immediately after its value — the batch-means standard error from §4, in the same
units (the writer emits `|value| · data_var`). BDO 2019 tolerates variance pages
but writes values only for now; its standard-error field is a planned addition.

### BDO 2019 (`bdo`, `bdo2019`, `binary`, `bin`; default)

For production use and multi-run accumulation. Values are written **exactly as
accumulated** (raw sums) and the primary count is embedded as the
`OSHBDO_RT_NSTAT` tag; the reader is responsible for normalisation. This
supports two merging paths:

1. **Native multi-worker** (future): each worker accumulates into its own
   accumulator set, the simulation layer merges them (§3), and one BDO file is
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

## 6. Unit handling and post-processing

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

## 7. Lifecycle and module boundaries

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

## References

- Welford, B. P. (1962). Note on a method for calculating corrected sums of
  squares and products. *Technometrics* 4(3), 419–420.
  doi:[10.1080/00401706.1962.10490022](https://doi.org/10.1080/00401706.1962.10490022)
- West, D. H. D. (1979). Updating mean and variance estimates: an improved
  method. *Communications of the ACM* 22(9), 532–535.
  doi:[10.1145/359146.359153](https://doi.org/10.1145/359146.359153)
- Schubert, E. & Gertz, M. (2018). Numerically stable parallel computation of
  (co-)variance. *Proc. 30th Int. Conf. on Scientific and Statistical Database
  Management (SSDBM '18)*, Article 10, 1–12.
  doi:[10.1145/3221269.3223036](https://doi.org/10.1145/3221269.3223036)
