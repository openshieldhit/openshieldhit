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

Letting one estimator dump several formats at once is **feasible and low-risk**.
The save layer is already a per-output dispatch keyed on a format string, the
writers are non-destructive readers of the post-processed accumulators, and the
runtime already stores each output's filename and format independently of the
cold workspace. The change is contained to the input parser, the cold output
struct, the compile step, and the save dispatch; the transport hot path is
untouched and — for the multi-format case — actually does *less* work than the
duplicate-block workaround it replaces.

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

**Feasible, low blast radius, hot path untouched.** The output model already
separates *what is scored* (geometry + quantities + filters → pages/accumulators
in `rt->pages[]`) from *how it is written* (`filename` + `fileformat` on the
output descriptor). Multi-format is the natural consequence of letting one scored
page-set feed more than one writer.

## Design

Two orthogonal decisions: **(A) user-facing syntax** and **(B) internal wiring**.

### (A) Input syntax — how a user asks for several formats

The blocking constraint: the ASCII and BDO writers use `out->filename`
**verbatim** — no extension is derived (only RTDOSE appends `.dcm`,
`osh_scoring_save_rtdose.c:88`). Two formats therefore need two distinct
filenames; the syntax must supply them.

| Option | Example | Notes |
|---|---|---|
| **A1. Repeatable `FileFormat <fmt> <file>`** | `FileFormat TEXT NB_msh.dat`<br>`FileFormat BDO NB_msh.bdo` | Explicit, unambiguous, order-independent. Each line pushes one *(format, filename)* target. Fully back-compat: a lone `FileFormat TEXT` + `Filename …` remains the single-target shorthand. **Recommended.** |
| **A2. Format list + base name + auto extension** | `Filename NB_msh`<br>`FileFormat TEXT BDO` | Terse, but *changes* filename semantics (today the name is verbatim, incl. extension). Needs a canonical extension per format and a rule for when a base already carries one → collision/overwrite risk. |
| **A3. New repeatable target keyword** | `Emit NB_msh.dat TEXT`<br>`Emit NB_msh.bdo BDO` | Same power as A1 with a distinct verb; clearer diff but adds a keyword. |

A1 is recommended: it reuses the existing `FileFormat` keyword, needs no
extension policy, and degrades to today's behaviour for the single-format case.
The one-arg form (`FileFormat TEXT`) keeps pairing with the block's `Filename`;
the two-arg form (`FileFormat TEXT foo.dat`) is self-contained.

### (B) Internal wiring — two strategies

Both build the page/accumulator set **once**; they differ in where the
"one page-set → many files" fan-out lives.

**B1. Fan-out at compile time (recommended — smallest blast radius).**
Keep `osh_scoring_output_runtime` exactly as it is (one file, one format, one
page-index list). In compile, expand an `Output` block carrying *K* targets into
*K* runtime outputs that **share the same page indices** — the pages (and their
accumulators) are built once and referenced by all *K* outputs (`page_indices`
is a cheap `size_t` copy; the accumulators in `rt->pages[]` are not duplicated).
Then:

- the transport hot path is byte-for-byte unchanged;
- the ASCII/BDO/RTDOSE writers are **unchanged** — they already read everything
  from `rt->outputs[i]`;
- the save dispatcher reads the format from `rt->outputs[i].fileformat`
  (already populated) instead of `ws->outputs[i].fileformat`, decoupling the
  runtime output list from the cold one;
- the `ws->noutputs == rt->noutputs` invariant (`osh_scoring_save.c:37`) and the
  writers' `output_idx < ws->noutputs` bounds checks relax to compare against
  `rt->noutputs`.

**B2. Fan-out at save time (targets on the descriptor).**
Change `osh_scoring_output_{def,runtime}` to hold an array of *(filename,
format)* targets over one shared `page_indices`. The save loop becomes
"for each output, for each target, dispatch a writer." This keeps
`rt->outputs` 1:1 with the cold blocks but changes the three writer signatures
(they must receive the specific target's filename+format rather than reading
`out->filename`). More writer churn for no user-visible gain over B1.

**Recommendation: A1 + B1.** B1 leaves the writers and the hot path alone and
localises the change to parse → cold struct → compile → dispatch.

## Blast radius (files touched, for A1 + B1)

| File | Change |
|---|---|
| `include/openshieldhit/scoring.h` | `osh_scoring_output_def`: replace the single `fileformat`/`filename` with a small `targets[]` array (`{char *filename; char *fileformat;}`), or keep the scalars as the "target 0" shorthand plus an `extra_targets[]`. |
| `src/apps/osh/osh_scoring_parse_output.c` | `output_fileformat`/`output_filename` push a target instead of overwriting; accept the optional filename arg on `FileFormat`. |
| `src/apps/osh/osh_scoring_parse.c` | Free the new target array on teardown. |
| `src/scoring/runtime/osh_scoring_compile.c` | Emit one `rt->outputs` entry per target, sharing the block's page indices; count pages once. |
| `src/scoring/save/osh_scoring_save.c` | Dispatch on `rt->outputs[i].fileformat`; relax the `ws/rt` count invariant. |
| `src/scoring/save/osh_scoring_save_ascii.c`, `…_bdo2019.c` | Relax the `ws->noutputs` bounds check to `rt->noutputs` (or drop the `ws` arg). |
| `src/scoring/osh_scoring.c` (`osh_scoring_estimate_memory`) | Count a block's pages once even with several targets (no per-format multiplier). |
| `docs/user/detect.dat.md` | Document the multi-format syntax. |
| `tests/unit/`, `tests/cases/` | A case that emits one geometry as TEXT+BDO from a single block; a unit test asserting both files appear and match the duplicate-block reference byte-for-byte. |

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

1. Ship **A1 + B1**: repeatable `FileFormat <fmt> [filename]`, fanned out into
   shared-page runtime outputs at compile time.
2. Keep **BDO the default** and single-format the common case — the feature is a
   pure superset; existing `detect.dat` files are unaffected.
3. Land it with a `tests/cases` fixture that emits one geometry as **TEXT + BDO
   from one block** and asserts the outputs are byte-identical to today's
   duplicate-block outputs, proving the accumulators were shared, not doubled.

## Open questions (for discussion with Niels)

1. **Syntax** — A1 (`FileFormat TEXT foo.dat`, recommended), A2 (base name +
   auto extension), or A3 (a new `Emit` verb)?
2. **Filename derivation** — if we ever want the terse A2 form, what is the
   canonical extension per format (`.dat`/`.txt`? `.bdo`? `.dcm`), and what
   happens when the base name already has an extension?
3. **RTDOSE in a mixed block** — reject, or allow and document that it consumes
   only the single dose page?
4. **Partial-write failure** — is a failure on the k-th target fatal for the
   whole run, or best-effort with a warning?
5. **Wiring** — B1 (compile-time fan-out, recommended) vs B2 (targets on the
   descriptor)? B1 keeps the writers untouched.
6. **Should the default stay one-format?** i.e. is there any appetite for a
   run-wide "always also emit BDO" switch, or is this strictly per-`Output`?

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
