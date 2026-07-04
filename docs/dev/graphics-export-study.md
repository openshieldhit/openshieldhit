# Study: exporting scoring results as viewable graphics (PNG / PDF / SVG)

Status: **research / design study, not yet implemented.**
Author: investigation on branch `claude/osh-graphics-export-study`.
Scope: evaluate adding a *native* image/plot output format to OpenShieldHIT so a
run can drop a ready-to-view picture (Bragg curve, dose map, spectrum) next to
its `.bdo` / `.dat`, **without** pulling in a heavy dependency or bloating the
binary.

Today the only way to get a picture is the external Python path
(`tools/plot_dicom.py` → matplotlib, the pymchelper style). That works but
requires a Python + matplotlib + numpy environment, i.e. it is *not*
self-contained and *not* available on a machine that only has the `openshieldhit`
binary.

---

## 1. TL;DR / recommendation

For a lean, MIT-licensed, Windows-portable, statically-linked C project that
**already hand-writes its own DICOM** and currently ships with **zero external
library dependencies in the core binary**, the on-brand answer is *not* to adopt
a plotting library. It is to write a few hundred lines of output code, and vendor
at most one tiny public-domain single-header file.

Recommended, in phases (each is independently shippable):

1. **`FileFormat SVG` — hand-rolled, zero dependencies (~250 LOC).**
   Covers 1-D line plots (depth–dose / Bragg curves, differential spectra with
   log axes) immediately. SVG is plain text (`fprintf`), viewable in any browser,
   and text/labels are *free* (native `<text>`). Highest value-to-effort. Ship
   this first.
2. **`FileFormat PDF` — multipage, vector.** Either hand-roll a minimal PDF
   writer (~400–500 LOC, 0 deps, uncompressed content streams) or vendor
   **PDFGen** (public domain, single `.c`+`.h`). This is what delivers the
   *multipage* deliverable the request specifically calls for (one detector or
   one slice per page).
3. **`FileFormat PNG` — 2-D heatmaps (dose maps).** Vendor
   **`stb_image_write.h`** (MIT / public domain, single header, built-in
   compressor, no zlib needed). Draw into an RGBA buffer with our own
   line/rect/colormap primitives, then one call writes the PNG. Optionally embed
   the same PNG as an image page inside the PDF so the PDF becomes the single
   unified deliverable.

Avoid: libharu, Cairo, PLplot, gnuplot-pipe, matplotlib-cpp, ImPlot — all bring
either heavy transitive dependencies or an external runtime/process, and
contradict the "one self-contained static binary" goal (§3).

### Measured footprint (this is the key number)

The release binary today is **660 KB on disk, ~594 KB of code** (`.text`).
I compiled the two lead single-file candidates and measured the *actual* linked
code they add (gcc `-O2 -ffunction-sections -Wl,--gc-sections`, PNG-only / vector
path):

| Option | Vendored source | Added code (`.text`) | Growth vs 594 KB |
|---|---:|---:|---:|
| Hand-rolled SVG (our code) | 0 (≈250 LOC ours) | ~3–5 KB | <1 % |
| Hand-rolled minimal PDF (our code) | 0 (≈450 LOC ours) | ~6–10 KB | ~1–1.5 % |
| **`stb_image_write.h`** (PNG path) | 1,724 LOC | **9.4 KB** (measured) | **+1.6 %** |
| **PDFGen** (multipage vector) | 6,353 LOC | **24 KB** (measured) | **+4 %** |
| libharu (static) | ~50–70k LOC | ~300–600 KB (est.) | +50–90 % |
| Cairo (+pixman +freetype +libpng +zlib) | huge | multiple MB + shared libs | ×several |

Takeaway: the tiny options are **noise** against a Monte-Carlo transport binary.
The heavyweight libraries would visibly bloat it. There is no reason to reach for
the heavyweight options for "simple plots to visualise the data".

---

## 2. What we would actually be plotting

Grounded in the real data model (`src/scoring/…`). A scored result is a set of
**pages**, each a flat `double` array indexed by a compiled geometry:

* `OSH_SCORING_GEO_MESH` — Cartesian, axes `X,Y,Z`, flat index
  `ix + nx*(iy + ny*iz)` (`osh_scoring_save_ascii.c:310`).
* `OSH_SCORING_GEO_CYL` — radial/axial, axes `R,Z`, flat index `ir + nr*iz`.
* Optional **differential** axes (`diff_nbins`, `diff2_nbins`) expand a page into
  a spectrum (`detect.dat.md`, "Differential scoring").

Each page carries a `quantity` name and a unit string (Dose MeV/g, DoseGy Gy,
Fluence 1/cm², DLET/TLET MeV/cm, …). Normalisation is the same rule the ASCII
writer already applies: NORM/SUM ÷ nstat, AVER (DLET/TLET) left as the physical
mean (`osh_scoring_save_ascii.c:186,333`).

That gives four concrete, useful plot types — all "simple", exactly as
requested:

| # | Plot | Data shape | Best format |
|---|---|---|---|
| 1 | **Depth–dose / Bragg curve**, any 1-D profile | MESH with 2 singleton axes, or CYL Z-column | **vector** (SVG/PDF): crisp line, tiny file |
| 2 | **Differential spectrum** (dΦ/dE, dose vs LET), often log–log | a page's `diff` axis | **vector** (SVG/PDF), log axes |
| 3 | **2-D dose map / heatmap** (XZ/XY/YZ plane, or CYL R/Z) | a non-singleton 2-D slice | **raster** (PNG): 1 texel/bin + colormap + colorbar |
| 4 | **Contact sheet**: many detectors/quantities in one file | several pages/outputs | **multipage PDF** |

Note the split: **1-D → vector**, **2-D → raster**. This is why the
recommendation is not a single format but a small menu. `tools/plot_dicom.py`
already does exactly type 3 (reshape page → plane → `imshow` with an `inferno`
colormap); we would be re-implementing that one specific behaviour in C, plus the
1-D line case which is arguably the more common need.

---

## 3. Constraints this repository imposes (why "just link matplotlib" is wrong)

These come straight from `CLAUDE.md`, `DEVELOPER.md`, and `llms.txt`:

* **Zero-dependency core is a design value, not an accident.** The core binary
  links no third-party library today (SDL2 is *examples-only*; the documented
  `.bdz`/zlib path is not actually wired into `src/`). DICOM CT + RTDOSE are
  hand-written specifically to avoid a DICOM library (`llms.txt:87`). A plotting
  solution must respect that or it will not be accepted.
* **MIT licensed.** Anything vendored must be license-compatible: **public
  domain / Unlicense / MIT / zlib / BSD are fine; GPL/LGPL are not** (LGPL static
  linking is a distribution headache; GPL is a non-starter). This rules out
  PLplot (LGPL) and gnuplot (its own restrictive-ish license).
* **Static, single binary, must build on Windows/MSVC.** No POSIX-only APIs
  (`DEVELOPER.md §11.1` bans `<unistd.h>`, `<sys/stat.h>`, `strcasecmp`,
  `mkstemp`, …). A candidate that needs `popen`, a shared `.so`, or a runtime
  interpreter fails this. (stb, PDFGen, and hand-rolled writers use only
  `stdio/stdlib/string/math` → clean.)
* **C11 + house style** for *our* wrapper code (block comments only, `struct foo`
  spelled out, `double const *p`, top-of-block decls). Vendored third-party
  headers are exempt as external code but must be kept unmodified, isolated in
  their own directory, and **excluded from clang-format / clang-tidy** (they will
  not pass §1–§15 and should not be expected to).
* **Hot-path allocation ban (§10) does not apply here.** Output/save runs once at
  the end (or at a snapshot), on the cold path — `malloc` for a pixel buffer or a
  PDF object list is completely fine. This is why plotting is a comfortable place
  to add code: none of the hard real-time rules bite.

---

## 4. Format trade-offs (PNG vs PDF vs SVG)

| | **SVG** | **PDF** | **PNG** |
|---|---|---|---|
| Kind | vector (XML text) | vector (+ can embed raster) | raster |
| Library needed | **none** (fprintf) | none (hand-roll) or 1 single-file | 1 single-file header |
| Text / axis labels | **native, free** | standard 14 fonts, ~free | must rasterize a font ourselves |
| 1-D line plots | excellent | excellent | ok |
| 2-D heatmaps | poor (1 `<rect>` per bin bloats) | ok (embed image) | **excellent** |
| Multipage | no | **yes** (the reason to pick PDF) | no |
| Direct viewing | any browser | any PDF viewer / browser | anything |
| Compression / size | text, gzips well | uncompressed streams OK | built-in deflate |
| External deps to *view* | none | none | none |

**Hidden cost — text rendering.** For *labeled* scientific plots (axis numbers,
titles, units) the deciding factor is text:

* **SVG**: `<text x=… y=…>12.5</text>` — free. ✅
* **PDF**: reference a built-in font (Helvetica) and emit `(...) Tj` — essentially
  free, no font embedding. ✅
* **PNG**: there is no text; we must rasterize glyphs ourselves. That means
  embedding a tiny bitmap font (a 5×7 or 8×8 ASCII table, ~1–2 KB) or vendoring
  `stb_easy_font.h`. This is the main extra work PNG carries.

This nuance is why the recommendation front-loads **SVG and PDF** for labeled
plots and treats **PNG** as the specialist for 2-D heatmaps (where the labels are
sparse: a title, axis ends, and a colorbar).

---

## 5. Library / approach survey (license + amount of code)

Measured line counts are from the actual current files where I fetched them;
others are well-known approximate sizes and are marked *(est.)*.

### Tier A — write it ourselves (0 vendored dependencies) — *recommended*

| Approach | License | Code | Formats | Notes |
|---|---|---|---|---|
| **Hand-rolled SVG writer** | ours (MIT) | ~150–250 LOC | SVG (1-D + simple 2-D) | pure `fprintf`; text free; browser-viewable; convert to PDF/PNG with external tools *if* wanted |
| **Hand-rolled minimal PDF** | ours (MIT) | ~300–500 LOC | PDF, **multipage**, vector | uncompressed streams (spec-legal), standard Helvetica font, image XObject for heatmaps; no deps, no zlib |
| **Hand-rolled PPM/PGM** | ours (MIT) | ~20 LOC | PPM/PGM | trivial but not broadly "viewable"; stepping-stone only |

### Tier B — vendor one tiny single-file public-domain/MIT header

| Library | License | Code | Formats | Deps | Measured add |
|---|---|---:|---|---|---:|
| **stb_image_write.h** | MIT / Public Domain (dual) | **1,724 LOC** | PNG, BMP, TGA, JPG, HDR | none (built-in deflate) | **+9.4 KB** (PNG path) |
| **PDFGen** (`pdfgen.c`+`.h`) | Public Domain (Unlicense) | **6,353 LOC** | PDF, multipage, lines/rects/text/JPEG-PNG-PPM embed | none | **+24 KB** |
| svpng | permissive (zlib-ish) | ~32 LOC | PNG (**uncompressed**) | none | tiny; bigger output files |
| Nayuki TinyPngOut | MIT | ~250 LOC | PNG (real deflate, minimal) | none | small |
| lodepng | zlib | ~6 kLOC *(est.)* single file | PNG enc/dec (own zlib) | none | heavier than stb for write-only |
| miniz | Public Domain (Unlicense) | ~10 kLOC *(est.)* single file | PNG write + real zlib/ZIP | none | overkill for plots *but* could also power the documented-yet-unbuilt `.bdz` deflate |
| fpng | MIT | ~2–3 kLOC *(est.)* C++ | PNG (fast) | none | **C++**; project is C |

### Tier C — real external libraries — *reject for this project*

| Library | License | Code / size | Why not here |
|---|---|---|---|
| libharu (libhpdf) | zlib/libpng | ~50–70 kLOC; static lib ~300–600 KB *(est.)* | full PDF, but a real dependency and ~doubles the binary; far more than "simple plots" needs |
| libpng + zlib | libpng + zlib | two libraries | the "standard" raster path, but two external libs vs stb's zero |
| Cairo (+pixman/freetype) | LGPL/MPL | multiple MB + shared libs | heavyweight; LGPL static-link friction; contradicts single-binary goal |
| PLplot | LGPL | large, many backends | LGPL; huge; built for a different scale |
| gnuplot (pipe via `popen`) | gnuplot license | external **program**, not a lib | needs gnuplot installed at runtime; `popen` is POSIX-ish; not self-contained |
| matplotlib-cpp | MIT header… | header is tiny **but embeds CPython + matplotlib** | reintroduces the exact Python/matplotlib runtime we are trying to escape |
| ImPlot / Dear ImGui | MIT | large, C++ | interactive GUI needing an GL/GPU context; wrong tool for batch file output |

---

## 6. How it plugs into osh (integration design)

The save subsystem is already a clean plugin point. Adding a plot format mirrors
the existing **RTDOSE** writer almost exactly.

**Dispatch** — `src/scoring/save/osh_scoring_save.c` picks a writer purely by the
`FileFormat` string (`fileformat_is_ascii/bdo2019/rtdose`). We add one predicate
and one call:

```c
/* in save_one_output() */
if (fileformat_is_plot(fileformat))          /* "png" | "pdf" | "svg" */
    return osh_scoring_save_plot_output(ws, rt, nstat, output_idx);
```

**No parser change is required.** `FileFormat` is stored as a free lowercased
string (`osh_scoring_parse_output.c:187`), so `FileFormat PNG` already parses; the
dispatcher just has to recognise it. (BDO2019 remains the default when the format
string is absent — see `fileformat_is_bdo2019(NULL) == 1`.)

**New writer** — one self-contained translation unit, added to
`src/scoring/CMakeLists.txt` next to the others:

```
src/scoring/save/osh_scoring_save_plot.c   /* dispatch by png/pdf/svg */
src/scoring/save/osh_scoring_save_plot.h
```

Everything it needs is already reachable from `rt` / `ws`, and the ASCII writer
is a ready template for pulling it out:

* axis bounds / bin counts via `mesh_axis_index()` + `geo->axes[i].{lo,hi,nbins}`
  (`osh_scoring_save_ascii.c:36,102,235`);
* page values `page->acc.data[idx]`, `page->quantity`, unit string, and the
  NORM/AVER normalisation rule (`:186,333`);
* classify plot type by counting non-singleton spatial axes (→ 1-D line vs 2-D
  heatmap) and checking `diff_nbins` (→ spectrum) — the same logic
  `plot_dicom.py:extract_planes` uses.

**Third-party isolation** — if we vendor stb/PDFGen, put it under
`src/common/third_party/` (or top-level `third_party/`), compile it into
`osh_scoring` (cold path, §10 irrelevant), and exclude it from
`tools/clang-format-all.sh` and clang-tidy. Keep the header pristine so it can be
bumped.

**Optional-feature gate** — to keep the default build lean and the dependency
opt-in, guard PNG behind a CMake option modeled on the existing
`OSH_ENABLE_COVERAGE` pattern (`CMakeLists.txt:95`):

```cmake
option(OSH_ENABLE_PLOT_PNG "Enable native PNG heatmap output (vendored stb)" OFF)
```

SVG (and a hand-rolled PDF) need no gate at all — they add no dependency.

**Small primitives we would write once** (shared by all backends, ~250–400 LOC
total): pick "nice" axis ticks (~40 LOC), map a value → color via an embedded
256×3 viridis/inferno LUT (~0.8 KB, matches `plot_dicom.py`), Bresenham line into
a raster buffer (PNG only), and axis/frame drawing.

---

## 7. Multipage PDF specifics

The request explicitly wants multipage PDF. PDF is the *only* one of the three
formats that supports it, and it maps naturally onto osh's data:

* **one page per Output/detector**, or **one page per Z-slice** of a 3-D mesh
  (a "flip-book" through the dose volume), or **one page per quantity**
  (Dose, then LET, then Fluence).
* Mechanically: PDF is a set of numbered objects; each page is a `/Page` object
  pointing at a content stream, all listed in the `/Pages` kids array. Adding a
  page is appending objects and one kid — hand-rolled or via
  `pdf_append_page()` in PDFGen.
* **No compression is required** — PDF content streams may be stored raw, so a
  hand-rolled multipage vector PDF needs **no zlib**. (Compression only matters
  if pages embed large raster images; a `FlateDecode` stream would then want a
  deflate, which is the one place miniz/stb's compressor could be reused.)
* A mixed document is possible and attractive: vector line-plot pages for 1-D
  profiles + image-XObject pages for 2-D heatmaps → a single self-contained
  report file per run.

---

## 8. Suggested phasing & effort

| Phase | Deliverable | New deps | Effort | Binary cost |
|---|---|---|---|---|
| 1 | `FileFormat SVG`: 1-D line plots (Bragg, spectra, log axes) | none | ~1 day, ~250 LOC | <1 % |
| 2 | `FileFormat PDF`: multipage vector (hand-rolled) **or** vendor PDFGen | none / PDFGen (PD) | ~1–2 days | ~1.5 % / ~4 % |
| 3 | `FileFormat PNG`: 2-D heatmaps + colorbar (vendor stb + bitmap font) | stb (MIT/PD) | ~2 days | ~1.6 % |
| 3b | embed PNG heatmap pages into the PDF report | (reuses 2+3) | small | — |

Each phase is independently useful and independently revertible. Phase 1 alone
already removes the "you need Python+matplotlib to see anything" problem for the
most common plot (the depth–dose curve).

---

## 9. Licenses appendix

* **stb_image_write.h** — dual **MIT** / **Public Domain (Unlicense)**, choose
  either. Compatible with osh's MIT. No attribution burden under the PD option.
* **PDFGen** — **Public Domain (Unlicense)**. Compatible.
* **svpng** — permissive (zlib-style). Compatible.
* **Nayuki TinyPngOut** — **MIT**. Compatible.
* **lodepng** — **zlib** license. Compatible.
* **miniz** — **Public Domain (Unlicense)**. Compatible.
* **libharu** — **zlib/libpng** license. Compatible, but rejected on size/deps.
* **PLplot** — **LGPL**. Avoid (static-link friction).
* **gnuplot** — own license, and it is a program not a library. Avoid.

All Tier-A/Tier-B options above are license-clean for an MIT project. The
decision is therefore driven by **code size and dependency philosophy, not
licensing** — and on those axes the hand-rolled + one-tiny-header path wins.

---

## 10. Sources

* stb — https://github.com/nothings/stb (`stb_image_write.h`, v1.16, 1,724 lines,
  dual MIT/Public-Domain, built-in deflate)
* PDFGen — https://github.com/AndreRenaud/PDFGen (`pdfgen.c` 5,515 + `pdfgen.h`
  838 lines, Unlicense/Public-Domain, multipage, lines/text/rects/image embed)
* svpng — https://github.com/miloyip/svpng (≈32-line uncompressed PNG writer)
* Nayuki TinyPngOut — https://www.nayuki.io/page/tiny-png-output (MIT)
* lodepng — https://github.com/lvandeve/lodepng (zlib license)
* miniz — https://github.com/richgel999/miniz (Unlicense; PNG write + zlib/ZIP)
* libharu — https://github.com/libharu/libharu (zlib/libpng license)
* Internal: `src/scoring/save/*` (save dispatch + writers),
  `src/apps/osh/osh_scoring_parse_output.c` (FileFormat parsing),
  `src/scoring/runtime/osh_scoring_geometry_runtime.h` (geometry/page model),
  `tools/plot_dicom.py` (existing matplotlib heatmap behaviour to mirror).
