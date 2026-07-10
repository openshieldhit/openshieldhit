<!--
Ready-to-file GitHub issue for the multi-format-per-estimator feature.
Title:  SCORING: attach multiple output formats to one estimator (dump BDO + TXT at once) — study & discussion
Kept in-repo as a draft so the write-up is reviewable in the diff; paste the
body below into a new issue when ready. Full analysis: docs/dev/multi-dataformat-attachment-study.md
-->

## Motivation

Today one scoring `Output` block (an *estimator*) writes exactly **one file in
one format** — a single `FileFormat` keyword picks ASCII **or** BDO 2019 **or**
DICOM RTDOSE. To get the same scored quantities out as both BDO **and** TXT, you
have to duplicate the whole `Output` block. That's exactly what
`tests/cases/04_simple_loaddedx/detect.dat` does today:

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

This works, but it is not "harmlessly verbose" — it **doubles the memory
allocated for scoring**, and the point of this issue is that it shouldn't have
to. The goal: attach several output formats to one estimator so it scores
**once** and writes **N** files, with the memory footprint independent of *N*.

## A worked example: what actually costs memory

Take a realistic scoring mesh: `100 × 100 × 200` bins (2 mm voxels over a
`20 × 20 × 40` cm volume), scoring two quantities, `Dose` and `Fluence`.

```
voxels     = 100 × 100 × 200        = 2,000,000
one page   = 2,000,000 × 8 bytes    = 16,000,000 bytes  = 16 MB   (Dose OR Fluence)
page-set   = 2 pages × 16 MB        = 32 MB                       (Dose AND Fluence)
```

That 32 MB — call it the block's **page-set** — is the *only* large,
configuration-driven allocation involved in scoring this mesh. Everything else
about an `Output` block (filename, format keyword, which pages it covers) is
tiny by comparison.

With the current duplicate-block workaround, getting this mesh out as **both**
TEXT and BDO means **two** `Output` blocks. Pages are never shared between
blocks — each block compiles its **own** 32 MB page-set — so the real cost is
64 MB, and the transport hot path deposits into *both* copies for every history
that crosses the mesh: double the per-step scoring work for these two
quantities, for the entire run, not just double the disk space.

| Formats requested | **Today** (duplicate `Output` blocks) | **Proposed** (one shared page-set) |
|---|---|---|
| 1 (TEXT only) | 32 MB | 32 MB |
| 2 (TEXT + BDO) | 64 MB | 32 MB + a few dozen bytes |
| 3 (TEXT + BDO + RTDOSE) | 96 MB | 32 MB + a few dozen bytes |
| 5 | 160 MB | 32 MB + a few dozen bytes |

**Today** grows linearly with the number of formats: *N* formats → *N* × 32 MB.
**Proposed** is flat: memory depends only on the geometry and quantities
scored — never on how many ways the result is written to disk. That's the
requirement this issue is really about: *scoring memory should not grow just
because more output formats are attached.*

### Where the "a few dozen bytes" comes from

An **output descriptor** — "write these pages to this file, in this format" —
never holds scored values itself. It holds a filename string, a lowercase
format keyword (`"text"`, `"bdo"`, …), and a short list of page indices (plain
integers, 8 bytes each) that point into the *one* shared page-set. Attaching a
second format to a block means adding one more of these tiny descriptors — not
a second 32 MB array.

### The "same photo, two prints" mental model

Think of a page-set as a photograph, developed once at the end of a run. Today,
wanting a JPEG print and a PNG print means taking the photo **twice** — one
exposure developed into a JPEG, a wholly separate exposure developed into a PNG.
Same image, double the film.

What we actually want is to take the photo **once** and hand that *same* photo
to two printers. Each printer looks at the same pixels and produces a different
file — one a JPEG encoder, one a PNG encoder — but neither printer owns its own
copy of the photo, and neither is allowed to touch the original pixels.

### How that maps onto code: loop over the array, don't copy it

A writer's core job is a **read-only loop** over a page's already-computed
values. Two writers can run that loop back-to-back over the identical array,
without either of them copying it first:

```c
/* TEXT writer (simplified): read-only pass over the shared accumulator */
for (i = 0; i < page->len; ++i) {
    double value = page->acc.data[i] / nstat;   /* ASCII normalises at write time */
    fprintf(fp, "%.6e\n", value);
}

/* BDO writer (simplified): a DIFFERENT read-only pass over the SAME array */
for (i = 0; i < page->len; ++i) {
    double value = page->acc.data[i];           /* BDO keeps the raw sum, tags nstat separately */
    fwrite(&value, sizeof(value), 1, fp);
}
```

Both loops read `page->acc.data` — the identical pointer, the identical
2,000,000 doubles. Nothing is copied to feed the second loop. The two writers
disagree only about **how to encode a value on the way out** (text row vs
binary token; ÷nstat vs raw-sum-plus-tag) — never about **where the value
lives**. That's exactly the mechanism this issue proposes to generalise: an
`Output` block with several `FileFormat` targets runs that loop once per
target, and every target reads from the one page-set built once, when the run
started, and this is genuinely how it already works today for a *single*
writer — the gap is only that today one page-set can feed exactly one writer.

## Why this is feasible right now (grounded in the code)

The good news: the codebase already keeps *what is scored* and *how it's
written* almost entirely separate. Save is a per-output dispatch on a plain
format string:

```c
/* src/scoring/save/osh_scoring_save.c */
fileformat = ws->outputs[output_idx].fileformat;
if (fileformat_is_ascii(fileformat))   return osh_scoring_save_ascii_output(ws, rt, nstat, output_idx);
if (fileformat_is_bdo2019(fileformat)) return osh_scoring_save_bdo2019_output(ws, rt, nstat, output_idx);
if (fileformat_is_rtdose(fileformat))  return osh_scoring_save_rtdose_output(ws, rt, nstat, output_idx);
```

Two properties make the "one page-set, many writers" model easy to add:

1. **Writers are non-destructive.** `osh_scoring_postprocess()` runs once
   (`osh_simulation.c:736`); each writer then reads `page->acc.data[...]` and
   applies its own scaling *at write time* (ASCII ÷ `nstat`; BDO writes raw
   sums + the `nstat` tag — `docs/dev/scoring.md` §5). A second writer can run
   immediately after over the same pages without either one mutating them.
2. **Writers barely touch `ws` (the cold workspace).** ASCII/BDO use it only for
   a bounds check; RTDOSE not at all. Everything substantive — filename,
   geometry, page list — is already read from `rt->outputs[i]`, which **already
   owns its own `filename`/`fileformat`** (`compile.c:991,996`) and reaches the
   shared accumulators only through `size_t` page indices.

## Proposed shape (recommendation — open to discussion)

**Syntax: a terse format list on one line.**

```text
Output
    Filename NB_msh          # stem (extension optional)
    FileFormat TEXT BDO      # -> NB_msh.dat  +  NB_msh.bdo
    Geo MyMesh
    Quantity Energy
    Quantity Fluence
```

ASCII/BDO write the filename **verbatim** today (only RTDOSE appends `.dcm`), so
with several formats sharing one `Filename` the writer needs to derive a
distinct path per format. Backward-compatible rule:

- **Single format → filename used verbatim.** Zero behaviour change; every
  existing `detect.dat` and every fixture under `tests/cases/` stays
  byte-for-byte unaffected.
- **Multiple formats → `Filename` is treated as a stem.** Strip a recognised
  extension (`.dat .txt .bdo .bdz .bin .dcm`) if present, then append the
  canonical one per format: `.dat` (TEXT), `.bdo`/`.bdz` (BDO), `.dcm` (RTDOSE).

  So `Filename NB_msh` + `FileFormat TEXT BDO` → `NB_msh.dat` + `NB_msh.bdo`.

Optional escape hatch, if wanted: `FileFormat TEXT NB_legacy.out` (a filename
after a single format keyword) overrides the derived name for that one target.

**Wiring: fan-out at compile time, not at save time.**
Build the block's pages/accumulators **once**, then emit one lightweight
runtime output *per format*, all sharing the same page-index list. The hot path
and all three existing writers stay untouched; the save dispatcher just reads
the format from each runtime output (already stores its own `fileformat`
today). This is the change that delivers the flat row of the table above — no
writer, and no accumulator, is ever duplicated because a second format was
requested.

## Blast radius

Touched: input parser (accept a format list on one `FileFormat` line), the cold
output struct (hold a list of format keywords alongside the filename stem),
compile (emit one lightweight runtime output per format, all sharing page
indices; derive per-format filenames when there's more than one), save dispatch
(read the format per runtime output; loosen a count check that currently
assumes one runtime output per cold block), the memory estimator (must count a
block's pages **once**, never × the number of formats — this is the actual
guarantee, and it needs a test), `docs/user/detect.dat.md`, and tests.

**Not** touched: the transport hot path, the per-step deposit code, the
accumulator/merge machinery, or the DICOM writer core. Snapshot/periodic dumps
inherit the feature for free since they reuse the same save dispatch.

## Open questions

1. **Extension policy** — confirm the canonical extension per format (`.dat`
   TEXT, `.bdo`/`.bdz` BDO, `.dcm` RTDOSE) and the stem-stripping set. Do we also
   want the optional per-target filename override (`FileFormat TEXT foo.out`)?
2. **RTDOSE in a mixed block** — the RTDOSE writer needs a voxel geometry and
   exactly one page. Reject a mixed block, or allow it and document that the
   RTDOSE target consumes only the dose page?
3. **Partial-write failure** — if the k-th target's file can't be opened, is
   that fatal for the whole run, or best-effort with a warning?
4. **How much redesign** — the minimal version shares page indices across
   lightweight per-format output records (delivers the flat-memory row above
   with the least code moved); a deeper version could also let *different*
   `Output` blocks share identical page-sets. Is that extra generality worth
   pursuing now, or is it a separate follow-up?
5. **Default** — strictly per-`Output`, or is there appetite for a run-wide
   "also emit BDO" switch?

## Related

- `docs/dev/scoring.md` §5 — per-format normalisation / multi-run merge rules
  (why the same accumulators can legitimately feed ASCII ÷nstat and BDO
  raw-sums at once — this is the existing precedent the proposal leans on).
- #238 — native graphics export: also adds writers behind the same
  `FileFormat` dispatch and raises the same "one file per Output vs combined
  report" UX question.
- #209 — MC standard error: shared accumulators emit variance consistently
  across a block's targets.

Full engineering write-up (struct-level detail, file-by-file blast radius,
alternative wiring): `docs/dev/multi-dataformat-attachment-study.md`.
