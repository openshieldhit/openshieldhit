# Multi-format output per estimator — feasibility study

> Status: **study for discussion**, not an accepted proposal. Companion to the
> GitHub issue of the same name. Written against the tree on branch
> `claude/multi-dataformat-attachment-study-4nhb9m`.

## Summary

Today one scoring `Output` block (an *estimator*) writes exactly **one file in
one format**: a single `FileFormat` keyword selects ASCII, BDO 2019, or DICOM
RTDOSE. To get the *same* scored quantities out as both BDO **and** TXT, a user
must duplicate the whole `Output` block. That duplication is not just verbose —
it **scores everything twice**: each block compiles its own accumulators and the
transport hot path deposits into every one of them, so a second format costs a
second full copy of the accumulator memory and a second deposit per step.

Letting one estimator dump several formats at once is **feasible and low-risk**,
and — done right — it costs **1× the accumulator memory regardless of how many
formats are written**. The save layer is already a per-output dispatch keyed on
a format string, the writers are non-destructive readers of the post-processed
accumulators, and the runtime already stores each output's filename and format
independently of the cold workspace. The change is contained to the input
parser, the cold output struct, the compile step, and the save dispatch; the
transport hot path is untouched and — for the multi-format case — actually does
*less* work than the duplicate-block workaround it replaces.

The design goal, stated up front, is a small **redesign of the output system**
that separates the two things a `FileFormat` line conflates today:

- the **scored page-set** — geometry + quantities + filters compiled into
  accumulators (`rt->pages[]`); this is the only large allocation, and it must be
  built **once**;
- the **write targets** — the *(format, filename)* pairs that read that page-set
  and emit files; these are cheap (a string, a keyword, and a `size_t` index
  list) and there can be many per page-set.

The desired user syntax is the terse `FileFormat TEXT BDO` (a list of formats on
one line), not one line per file.

## Current architecture (how a format attaches today)

### Cold (parsed) form

`detect.dat` `Output` blocks parse into `struct osh_scoring_output_def`
(`include/openshieldhit/scoring.h:157`):

```c
struct osh_scoring_output_def {
    char *filename;      /* Output file name. */
    char *geometry_name; /* Referenced geometry name. */
    char *fileformat;    /* Optional format keyword; NULL = default BDO. */
    struct osh_scoring_page_def *pages;
    size_t npages;
};
```

`filename` and `fileformat` are **single scalars**. The parser
(`src/apps/osh/osh_scoring_parse_output.c`) frees and overwrites each on every
`Filename` / `FileFormat` line, so repeating either today just keeps the last
value:

```c
static enum osh_status output_fileformat(struct osh_scoring_output_def *out, ...) {
    ...
    free(out->fileformat);
    out->fileformat = strdup(words[1]);
    osh_lower_inplace(out->fileformat);
    ...
}
```

### Runtime form

Compile (`src/scoring/runtime/osh_scoring_compile.c:988`) copies each cold output
into one `struct osh_scoring_output_runtime`
(`src/scoring/runtime/osh_scoring_output_runtime.h:102`):

```c
struct osh_scoring_output_runtime {
    char *filename;       /* owned */
    char *fileformat;     /* owned; NULL = default BDO */
    size_t geometry_idx;
    size_t *page_indices; /* ordered indices into the flat rt->pages[] */
    size_t npages;
};
```

Crucially, **the runtime output already owns its own copies of `filename` and
`fileformat`** (`compile.c:991`, `compile.c:996`) — they are not aliases into the
cold workspace. The page *accumulators* live in the flat `rt->pages[]` array;
`page_indices` are just `size_t` indices into it.

### Save dispatch

`osh_scoring_save_outputs()` (`src/scoring/save/osh_scoring_save.c:23`) loops over
outputs and dispatches each to exactly one writer by string match:

```c
fileformat = ws->outputs[output_idx].fileformat;
if (fileformat_is_ascii(fileformat))   return osh_scoring_save_ascii_output(ws, rt, nstat, output_idx);
if (fileformat_is_bdo2019(fileformat)) return osh_scoring_save_bdo2019_output(ws, rt, nstat, output_idx);
if (fileformat_is_rtdose(fileformat))  return osh_scoring_save_rtdose_output(ws, rt, nstat, output_idx);
return OSH_ENOTSUP;
```

Two facts make multi-format easy:

1. **The writers are non-destructive** on the accumulators. `osh_scoring_postprocess()`
   runs once (`src/simulation/osh_simulation.c:736`); each writer then reads
   `page->acc.data[...]` and applies its own scaling *at write time* (ASCII
   divides by `nstat`; BDO writes raw sums plus the `nstat` tag —
   `docs/dev/scoring.md` §5). Nothing a writer does prevents a second writer from
   running immediately after over the same post-processed pages.
2. **The writers barely use `ws`.** ASCII and BDO touch `ws` only for a
   `output_idx < ws->noutputs` bounds check; RTDOSE ignores `ws` entirely. Every
   substantive field — filename, geometry, page list — is read from
   `rt->outputs[output_idx]`. The dispatcher's read of
   `ws->outputs[output_idx].fileformat` has an identical sibling at
   `rt->outputs[output_idx].fileformat`.

## The real cost of the workaround

The current way to emit both formats is to duplicate the block (this is exactly
what `tests/cases/04_simple_loaddedx/detect.dat` does):

```text
Output
   Filename NB_msh.dat
   Fileformat TEXT
   Geo MyMesh
   Quantity Energy
   Quantity Fluence

Output
   Filename NB_msh.bdo      # default BDO
   Geo MyMesh
   Quantity Energy
   Quantity Fluence
```

Pages are **not** deduplicated across outputs. Each `Output` block contributes
its own pages to the flat `rt->pages[]` array (`compile.c:988`–`1020`, one
prepared page per *(output, quantity)* pair), and each page owns its own
`struct osh_scoring_accumulator`. The hot path groups pages by
*(geometry, score_kind)* and deposits into **every** page of a group
(`osh_scoring_geometry_score_group`, `osh_scoring_geometry_runtime.h:35`;
loop in the estimator step handlers). So the duplicate-block pattern above:

- **doubles the scoring accumulator memory** for `Energy` + `Fluence` on
  `MyMesh` (two independent `data[]` arrays holding identical numbers), and
- **doubles the per-step deposit work** for those quantities during transport.

It is also an easy correctness trap: the two quantity lists can silently drift
apart when someone edits one block and forgets the other.

A first-class multi-format estimator scores **once** and writes **N** files,
eliminating both costs.

## Feasibility verdict

**Feasible, low blast radius, hot path untouched, 1× memory.** The output model
already separates *what is scored* (geometry + quantities + filters →
pages/accumulators in `rt->pages[]`) from *how it is written* (`filename` +
`fileformat` on the output descriptor). Multi-format is the natural consequence
of letting one scored page-set feed more than one writer — and because the
accumulators live in the shared flat `rt->pages[]` array and outputs reference
them only by `size_t` index, adding formats never allocates another accumulator.

## Memory: the guarantee

This is the headline requirement, so state it precisely. The large,
configuration-driven allocation is the per-page accumulator — `data[]` (plus
optional `data2[]`/variance) sized `spatial_bins × diff1_bins × diff2_bins ×
sizeof(double)` per page (`osh_scoring_estimate_memory`, `include/openshieldhit/
scoring.h:92`). An **output descriptor** is by comparison free: a filename
string, a lowercased format keyword, a geometry index, and a `size_t[]` list of
page indices.

- **Today's duplicate-block workaround: O(formats × accumulators).** Each block
  compiles its own pages, so TEXT+BDO on one mesh holds *two* identical `data[]`
  arrays and deposits into both every step.
- **Multi-format estimator: O(accumulators), independent of formats.** The block
  compiles one page-set; each format is one more cheap descriptor pointing at the
  *same* pages. Two formats, ten formats — same accumulator bytes, same deposits.

`osh_scoring_estimate_memory()` must reflect this: count a block's pages **once**,
never multiply by the number of formats.

## Design

Two orthogonal decisions: **(A) user-facing syntax** and **(B) internal wiring**.

### (A) Input syntax — the terse format list

Target syntax is a **list of formats on the `FileFormat` line**:

```text
Output
    Filename NB_msh          # stem (extension optional)
    FileFormat TEXT BDO      # -> NB_msh.dat  +  NB_msh.bdo
    Geo MyMesh
    Quantity Energy
    Quantity Fluence
```

The blocking constraint is filenames: ASCII and BDO write `out->filename`
**verbatim** — no extension is derived (only RTDOSE appends `.dcm`,
`osh_scoring_save_rtdose.c:88`). With several formats sharing one `Filename`, the
writer must derive a distinct path per format. Rule that keeps existing files
working:

- **Single format (the common case): filename used verbatim — zero behaviour
  change.** `FileFormat TEXT` + `Filename NB_msh.dat` still writes exactly
  `NB_msh.dat`.
- **Multiple formats: `Filename` is a stem.** Strip a recognised trailing
  extension (`.dat` `.txt` `.bdo` `.bdz` `.bin` `.dcm`) if present, then append
  the canonical extension per format:

  | Format keyword(s) | Canonical extension |
  |---|---|
  | `TEXT` / `ASCII` / `TXT` / `DAT` | `.dat` |
  | `BDO` / `BDO2019` / `BINARY` / `BIN` | `.bdo` (`.bdz` when compressed) |
  | `RTDOSE` | `.dcm` |

  So `Filename NB_msh` + `FileFormat TEXT BDO` → `NB_msh.dat` + `NB_msh.bdo`, and
  `Filename NB_msh.dat` + `FileFormat TEXT BDO` → the same (the `.dat` stem is
  stripped first).

Two optional escape hatches, if the maintainer wants them (not required for the
core feature):

- **Explicit per-format filename** — `FileFormat TEXT NB_msh.dat` (a filename
  after a single format keyword) overrides the derived name for that one target.
  Useful when a legacy name doesn't match the canonical extension.
- Keeping the single-format `FileFormat`/`Filename` pair as-is means every
  current `detect.dat` — and every fixture under `tests/cases/` — is byte-for-byte
  unaffected.

### (B) Internal wiring — decouple page-set from targets

Both viable strategies build the page/accumulator set **once**; they differ in
where the "one page-set → many files" fan-out lives.

**B1. Fan-out at compile time (recommended — smallest blast radius).**
Keep `osh_scoring_output_runtime` as-is (one file, one format, one page-index
list). In compile, expand an `Output` block carrying *K* formats into *K* runtime
outputs that **share the same page indices** — the pages (and their accumulators)
are built once and referenced by all *K* outputs (`page_indices` is a cheap
`size_t` copy; nothing in `rt->pages[]` is duplicated). Then:

- the transport hot path is byte-for-byte unchanged;
- the ASCII/BDO/RTDOSE writers are **unchanged** — they already read everything
  from `rt->outputs[i]`;
- the save dispatcher reads the format from `rt->outputs[i].fileformat`
  (already populated) instead of `ws->outputs[i].fileformat`, decoupling the
  runtime output list from the cold one;
- the `ws->noutputs == rt->noutputs` invariant (`osh_scoring_save.c:37`) and the
  writers' `output_idx < ws->noutputs` bounds checks relax to compare against
  `rt->noutputs`.

**B2. Fan-out at save time (an explicit page-set / target split).**
The more thorough "redesign" reading: give the runtime a flat list of **page-sets**
(the accumulators) and a separate flat list of **write targets**, each target
referencing a page-set by index. `osh_scoring_output_{def,runtime}` grows an array
of *(filename, format)* targets over one shared `page_indices`, and the save loop
becomes "for each output, for each target, dispatch a writer." This models the
decoupling most explicitly and leaves room to later dedup *identical* page-sets
across different `Output` blocks, but it changes the three writer signatures (they
must be handed the specific target's filename+format instead of reading
`out->filename`). More churn now for headroom later.

**Recommendation: B1.** It already delivers the 1× memory guarantee and the terse
syntax with the least code moved (parse → cold struct → compile → dispatch), and
leaves the writers and hot path untouched. B2's extra generality (cross-block
page-set dedup) is a separate, later improvement and need not block this feature.

## Blast radius (files touched, for the list-syntax + B1 plan)

| File | Change |
|---|---|
| `include/openshieldhit/scoring.h` | `osh_scoring_output_def`: hold a list of format keywords (`char **fileformats` + count) alongside the existing `filename` stem; keep the single `fileformat` scalar as the one-element shorthand. |
| `src/apps/osh/osh_scoring_parse_output.c` | `output_fileformat` accepts several keywords on one line (push each); optional trailing filename overrides the derived name for a single-format line. |
| `src/apps/osh/osh_scoring_parse.c` | Free the new format list on teardown. |
| `src/scoring/runtime/osh_scoring_compile.c` | Emit one `rt->outputs` entry per format, sharing the block's page indices; derive per-format filenames (stem + canonical extension) when >1 format; count pages once. |
| `src/scoring/save/osh_scoring_save.c` | Dispatch on `rt->outputs[i].fileformat`; relax the `ws/rt` count invariant. |
| `src/scoring/save/osh_scoring_save_ascii.c`, `…_bdo2019.c` | Relax the `ws->noutputs` bounds check to `rt->noutputs` (or drop the `ws` arg). |
| `src/scoring/osh_scoring.c` (`osh_scoring_estimate_memory`) | Count a block's pages once regardless of format count (no per-format multiplier) — the memory guarantee. |
| `docs/user/detect.dat.md` | Document the `FileFormat f1 f2` list syntax and the stem/extension rule. |
| `tests/unit/`, `tests/cases/` | A case that emits one geometry as TEXT+BDO from a single block; a unit test asserting both files appear, match the duplicate-block reference byte-for-byte, and that `osh_scoring_estimate_memory` is unchanged vs the single-format block (proving accumulators are shared, not doubled). |

The transport hot path, `osh_scoring_step.*`, the accumulator/merge machinery,
and the DICOM writer core are **not** touched.

## Edge cases and constraints

- **Filename collisions.** Two targets must resolve to different paths; the
  compile step should reject duplicates with a clear diagnostic (matching the
  existing `OSH_DIAG_ERRORF` style). A2's auto-extension scheme makes this
  subtler than A1's explicit names.
- **RTDOSE is special.** The RTDOSE writer requires a voxel/DICOM geometry and
  exactly one page (`osh_scoring_save_rtdose.c:45`, `out->npages != 1`). A block
  that mixes `RTDOSE` with a multi-page `TEXT`/`BDO` target must either be
  rejected or documented as "the RTDOSE target sees only the single dose page."
  Worth an explicit rule, not silent behaviour.
- **Snapshot / periodic dumps.** The mid-run dump path
  (`osh_scoring_sink.c` → `osh_scoring_save_outputs`) shares the same dispatch,
  so multi-format "just works" for `--dump-every*` / `NSTAT nsave` once the
  dispatch fans out. No separate work.
- **Memory estimate.** `osh_scoring_estimate_memory()` sums pages per block; with
  shared page-sets it must **not** multiply by the target count — the whole point
  is that memory is independent of how many formats are written.
- **Variance / error bars.** Per-page variance state (issue #209) is a property
  of the shared accumulators, so it is written consistently by whichever targets
  support it (ASCII today; BDO's σ field is still a planned addition per
  `docs/dev/scoring.md` §5). No new interaction, but worth a test.
- **Ordering / atomicity.** Writing N files is N separate `fopen`s; a failure on
  the 2nd should surface which target failed and should not leave the 1st file
  implying success of the whole set. Define whether a partial failure is fatal.

## Recommendation

1. Adopt the terse **`FileFormat TEXT BDO`** list syntax with the stem +
   canonical-extension rule (verbatim filename preserved for the single-format
   case, so nothing existing changes).
2. Wire it with **B1** (compile-time fan-out into shared-page runtime outputs):
   **1× accumulator memory independent of format count**, hot path and writers
   untouched.
3. Keep **BDO the default** and single-format the common case — the feature is a
   pure superset; existing `detect.dat` files are unaffected.
4. Land it with a `tests/cases` fixture that emits one geometry as **TEXT + BDO
   from one block**, asserts the outputs are byte-identical to today's
   duplicate-block outputs, **and** asserts `osh_scoring_estimate_memory()` is
   unchanged from the single-format block — the memory guarantee, enforced by a
   test.

## Open questions (for discussion with Niels)

1. **Extension policy for the list form** — confirm the canonical extension per
   format (`.dat` for TEXT, `.bdo`/`.bdz` for BDO, `.dcm` for RTDOSE) and the
   stem-stripping set. Do we also want the optional `FileFormat TEXT foo.dat`
   per-target filename override, or is the derived name always sufficient?
2. **RTDOSE in a mixed block** — the RTDOSE writer needs a voxel geometry and
   exactly one page (`osh_scoring_save_rtdose.c:45`). Reject a mixed block, or
   allow it and document that the RTDOSE target consumes only the dose page?
3. **Partial-write failure** — is a failure on the k-th target fatal for the
   whole run, or best-effort with a warning?
4. **Wiring depth** — B1 (compile-time fan-out, recommended: least code, meets
   the memory goal) vs B2 (explicit page-set/target split, more churn but opens
   cross-block accumulator dedup later). Is cross-block dedup a goal worth the
   extra generality now?
5. **Run-wide switch?** — strictly per-`Output`, or is there appetite for a
   global "always also emit BDO/TEXT" flag?

## Related work

- `docs/dev/scoring.md` §5 — per-format normalisation and multi-run merge rules
  (why ASCII ÷ nstat vs BDO raw-sums-plus-nstat matters when the *same*
  accumulators feed both).
- Issue #238 — native graphics export; also adds writers behind the same
  `FileFormat` dispatch and raises the parallel "one file per Output vs one
  combined report" UX question.
- Issue #209 — MC standard-error tracking; shared accumulators mean variance is
  emitted consistently across a block's targets.

## Source references

- `include/openshieldhit/scoring.h:157` — cold `osh_scoring_output_def`.
- `src/apps/osh/osh_scoring_parse_output.c` — `Filename`/`FileFormat`/`Quantity`
  parsing (single-scalar overwrite today).
- `src/scoring/runtime/osh_scoring_output_runtime.h:102` — runtime output.
- `src/scoring/runtime/osh_scoring_compile.c:988` — cold → runtime output copy.
- `src/scoring/save/osh_scoring_save.c` — per-output format dispatch.
- `src/scoring/save/osh_scoring_save_{ascii,bdo2019,rtdose}.c` — the writers
  (non-destructive; read `rt->outputs[i]`).
- `src/scoring/save/osh_scoring_sink.c` — snapshot/dump reuse of the dispatch.
- `tests/cases/04_simple_loaddedx/detect.dat` — the duplicate-block workaround in
  the wild.
