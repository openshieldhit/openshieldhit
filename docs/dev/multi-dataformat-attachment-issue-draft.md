<!--
Ready-to-file GitHub issue for the multi-format-per-estimator feature.
Title:  SCORING: attach multiple output formats to one estimator (dump BDO + TXT at once) — study & discussion
Kept in-repo as a draft so the write-up is reviewable in the diff; paste the
body below into a new issue when ready. Full analysis: docs/dev/multi-dataformat-attachment-study.md
-->

## Motivation

Today one scoring `Output` block (an *estimator*) writes exactly **one file in
one format** — a single `FileFormat` keyword picks ASCII **or** BDO 2019 **or**
DICOM RTDOSE. To get the same scored quantities out as both BDO **and** TXT you
have to duplicate the whole `Output` block (this is exactly what
`tests/cases/04_simple_loaddedx/detect.dat` does).

That workaround is worse than it looks. Pages are **not** deduplicated across
`Output` blocks: each block compiles its own accumulators (`compile.c:988`–
`1020`, one page per *(output, quantity)* pair) and the transport hot path
deposits into **every** page of a *(geometry, score_kind)* group
(`osh_scoring_geometry_runtime.h:35`). So a second copy of the block to get a
second format **doubles the accumulator memory and doubles the per-step deposit
work** for those quantities — and invites the two quantity lists to silently
drift apart.

Letting one estimator dump several formats at once would score **once** and write
**N** files — at **1× the accumulator memory regardless of format count**. It's a
nice-to-have, and it turns out to be cheap. Full write-up committed as
`docs/dev/multi-dataformat-attachment-study.md` on branch
`claude/multi-dataformat-attachment-study-4nhb9m`.

## The memory point (the priority)

The large allocation is the per-page accumulator — `data[]` (+ optional
`data2[]`/variance), sized `spatial_bins × diff1 × diff2 × sizeof(double)`. An
*output descriptor* is trivially cheap: a filename, a format keyword, and a
`size_t[]` list of page indices into the shared `rt->pages[]`.

- **Duplicate-block workaround: O(formats × accumulators)** — TEXT+BDO holds two
  identical `data[]` arrays and deposits into both every step.
- **Multi-format estimator: O(accumulators), independent of formats** — one
  page-set, one more cheap descriptor per format, all pointing at the same pages.

So the ask is a small **redesign of the output system**: separate the *scored
page-set* (built once) from the *(format, filename) write targets* (many, cheap).
`osh_scoring_estimate_memory()` must then count a block's pages **once**, never
× the format count.

## Why it's feasible (grounded in the code)

The output model already separates *what is scored* from *how it's written*, and
the save layer is a per-output dispatch on a format string:

```c
/* src/scoring/save/osh_scoring_save.c */
fileformat = ws->outputs[output_idx].fileformat;
if (fileformat_is_ascii(fileformat))   return osh_scoring_save_ascii_output(ws, rt, nstat, output_idx);
if (fileformat_is_bdo2019(fileformat)) return osh_scoring_save_bdo2019_output(ws, rt, nstat, output_idx);
if (fileformat_is_rtdose(fileformat))  return osh_scoring_save_rtdose_output(ws, rt, nstat, output_idx);
```

Two properties make fan-out easy:

1. **Writers are non-destructive.** `osh_scoring_postprocess()` runs once
   (`osh_simulation.c:736`); each writer then reads `page->acc.data[...]` and
   applies its scaling at write time (ASCII ÷ `nstat`; BDO writes raw sums + the
   `nstat` tag — `docs/dev/scoring.md` §5). A second writer can run straight
   after over the same pages.
2. **Writers barely use `ws`.** ASCII/BDO touch `ws` only for a bounds check;
   RTDOSE not at all. Everything substantive comes from `rt->outputs[i]`, which
   **already owns its own `filename` and `fileformat`** (`compile.c:991,996`) and
   references shared accumulators through `size_t` `page_indices`.

## Proposed shape (recommendation — open to discussion)

**Syntax: terse format list `FileFormat TEXT BDO`.**

```text
Output
    Filename NB_msh          # stem (extension optional)
    FileFormat TEXT BDO      # -> NB_msh.dat  +  NB_msh.bdo
    Geo MyMesh
    Quantity Energy
    Quantity Fluence
```

ASCII/BDO use the filename **verbatim** today (only RTDOSE appends `.dcm`), so
with several formats sharing one `Filename` the writer derives a path per format.
Backward-compatible rule:

- **Single format → filename verbatim** (zero behaviour change; every existing
  `detect.dat` and `tests/cases/` fixture is byte-for-byte unaffected).
- **Multiple formats → `Filename` is a stem**: strip a recognised extension
  (`.dat .txt .bdo .bdz .bin .dcm`) if present, append the canonical one per
  format — `.dat` (TEXT), `.bdo`/`.bdz` (BDO), `.dcm` (RTDOSE).

Optional escape hatch if wanted: `FileFormat TEXT NB_legacy.out` (a filename
after a single keyword) overrides the derived name for that target.

**Wiring (B1): fan-out at compile time.**
Build the block's pages/accumulators **once**, then emit one
`osh_scoring_output_runtime` per format sharing the same `page_indices`. The hot
path and all three writers are untouched; the save dispatcher reads the format
from `rt->outputs[i].fileformat` (already populated) and the
`ws->noutputs == rt->noutputs` invariant relaxes to `rt->noutputs`. (A more
explicit page-set/target split — B2 in the study — is possible if we later want
to dedup identical page-sets *across* blocks, but it isn't needed to hit the
memory goal.)

## Blast radius

Touched: input parser (`osh_scoring_parse_output.c` — accept a format list),
cold struct (`osh_scoring_output_def` in `include/openshieldhit/scoring.h` — hold
a format list + filename stem), compile (`osh_scoring_compile.c` — one runtime
output per format sharing page indices, derive per-format filenames), save
dispatch (`osh_scoring_save.c` + relax the `ws->noutputs` bounds checks in the
ASCII/BDO writers), memory estimate (`osh_scoring_estimate_memory` — count a
block's pages **once**, the guarantee), `docs/user/detect.dat.md`, and tests
(TEXT+BDO from one block; assert byte-identical to the duplicate-block output and
unchanged memory estimate).

**Not** touched: the transport hot path, `osh_scoring_step.*`, the
accumulator/merge machinery, the DICOM writer core. Snapshot/periodic dumps
(`osh_scoring_sink.c`) inherit the feature for free via the same dispatch.

## Open questions

1. **Extension policy** — confirm the canonical extension per format (`.dat`
   TEXT, `.bdo`/`.bdz` BDO, `.dcm` RTDOSE) and the stem-stripping set. Do we also
   want the optional per-target filename override (`FileFormat TEXT foo.out`)?
2. **RTDOSE in a mixed block** — the RTDOSE writer needs a voxel geometry and
   exactly one page (`osh_scoring_save_rtdose.c:45`). Reject a mixed block, or
   allow it and document that the RTDOSE target consumes only the dose page?
3. **Partial-write failure** — if the k-th target's `fopen` fails, is that fatal
   for the run, or best-effort with a warning?
4. **Wiring depth** — B1 (compile-time fan-out, recommended: least code, meets
   the memory goal) vs B2 (explicit page-set/target split, opening cross-block
   accumulator dedup later)? Is cross-block dedup a goal worth the extra
   generality now?
5. **Default** — strictly per-`Output`, or is there appetite for a run-wide
   "also emit BDO" switch?

## Related

- `docs/dev/scoring.md` §5 — per-format normalisation / multi-run merge (why the
  same accumulators can legitimately feed ASCII ÷nstat and BDO raw-sums at once).
- #238 — native graphics export: also adds writers behind the same `FileFormat`
  dispatch and raises the same "one file per Output vs combined report" UX split.
- #209 — MC standard error: shared accumulators emit variance consistently
  across a block's targets.
