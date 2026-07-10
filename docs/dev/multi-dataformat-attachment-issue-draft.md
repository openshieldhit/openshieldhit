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
**N** files. It's a nice-to-have, and it turns out to be cheap. Full write-up
committed as `docs/dev/multi-dataformat-attachment-study.md` on branch
`claude/multi-dataformat-attachment-study-4nhb9m`.

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

**Syntax (A1): repeatable `FileFormat <fmt> [filename]`.**

```text
Output
    Geo MyMesh
    Quantity Energy
    Quantity Fluence
    FileFormat TEXT NB_msh.dat
    FileFormat BDO  NB_msh.bdo
```

Each `FileFormat` line pushes one *(format, filename)* target. A lone
`FileFormat TEXT` + `Filename …` stays the single-target shorthand, so every
existing `detect.dat` is unaffected and BDO stays the default. The two-arg form
is needed because ASCII/BDO use the filename **verbatim** — no extension is
derived (only RTDOSE appends `.dcm`), so two formats need two names.

**Wiring (B1): fan-out at compile time.**
Build the block's pages/accumulators **once**, then emit one
`osh_scoring_output_runtime` per target sharing the same `page_indices`. The hot
path and all three writers are untouched; the save dispatcher reads the format
from `rt->outputs[i].fileformat` (already populated) and the
`ws->noutputs == rt->noutputs` invariant relaxes to `rt->noutputs`.

## Blast radius

Touched: input parser (`osh_scoring_parse_output.c`), cold struct
(`osh_scoring_output_def` in `include/openshieldhit/scoring.h`), compile
(`osh_scoring_compile.c`), save dispatch (`osh_scoring_save.c` + bounds checks in
the ASCII/BDO writers), memory estimate (`osh_scoring_estimate_memory` — count a
block's pages **once**), `docs/user/detect.dat.md`, and tests.

**Not** touched: the transport hot path, `osh_scoring_step.*`, the
accumulator/merge machinery, the DICOM writer core. Snapshot/periodic dumps
(`osh_scoring_sink.c`) inherit the feature for free via the same dispatch.

## Open questions

1. **Syntax** — A1 (`FileFormat TEXT foo.dat`, recommended), A2 (base name +
   auto-derived extension), or A3 (a new `Emit <file> <fmt>` verb)?
2. **RTDOSE in a mixed block** — the RTDOSE writer needs a voxel geometry and
   exactly one page (`osh_scoring_save_rtdose.c:45`). Reject a mixed block, or
   allow it and document that the RTDOSE target consumes only the dose page?
3. **Partial-write failure** — if the k-th target's `fopen` fails, is that fatal
   for the run, or best-effort with a warning?
4. **Filename policy** — if we ever add the terse A2 form, what's the canonical
   extension per format, and what happens when the base name already has one?
5. **Default** — strictly per-`Output`, or is there appetite for a run-wide
   "also emit BDO" switch?

## Related

- `docs/dev/scoring.md` §5 — per-format normalisation / multi-run merge (why the
  same accumulators can legitimately feed ASCII ÷nstat and BDO raw-sums at once).
- #238 — native graphics export: also adds writers behind the same `FileFormat`
  dispatch and raises the same "one file per Output vs combined report" UX split.
- #209 — MC standard error: shared accumulators emit variance consistently
  across a block's targets.
