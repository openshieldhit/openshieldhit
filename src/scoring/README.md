# Scoring Module (`src/scoring/`)

This module owns scoring-domain data and logic: compiling cold detector and
filter definitions into runtime **accumulators**, depositing per-step
contributions during transport, post-run normalisation, and writing results
out (ASCII, BDO 2019, DICOM RTDOSE). It never touches files itself — the app
layer parses input and hands the save backends fully-resolved paths.

## What is where

| Path | Contents |
|---|---|
| `osh_scoring.{c,h}` | Cold scoring workspace API (`struct osh_scoring_workspace`) and the public detector/filter definition helpers. |
| `runtime/osh_scoring_compile.*` | Compile cold definitions into runtime accumulators and scoring geometry. |
| `runtime/osh_scoring_accumulator.*` | Owning accumulator storage (`data`/`data2`/variance), the single `osh_score_deposit()` write seam, and the `osh_scoring_accumulator_merge()` reduce. |
| `runtime/osh_scoring_step.*` | Per-step deposition called from transport (the scoring hot path). |
| `runtime/osh_scoring_postprocess.*` | Post-run unit conversion and two-pass average finalisation. |
| `runtime/osh_scoring_runtime.h` | Runtime layout consumed by transport and simulation. |
| `save/` | Output writers dispatched by `osh_scoring_save()`: ASCII (`text`/`txt`/`ascii`/`dat`), BDO 2019 (`bdo`/`bdo2019`/`binary`/`bin`; default), and DICOM RTDOSE (`rtdose`). |

## Why it is built this way (in brief)

- **Cold → hot split.** Parsed scoring definitions stay separate from the
  runtime accumulator buffers; the save layer consumes the workspace plus
  scoring runtime, never transport internals.
- **Storage separated from descriptor.** Each page's mutable, per-history
  arrays live in a `struct osh_scoring_accumulator` (`page->acc`), distinct
  from the page *descriptor* (geometry indices, strides, differential-axis
  config). That storage is the unit a future parallel worker owns privately and
  folds back with `osh_scoring_accumulator_merge()`.
- **One deposit seam.** Every tally funnels through `osh_score_deposit()` —
  today a plain `+=`, tomorrow the single place an atomic / private / locked
  write policy plugs in.
- **Two output stances on normalisation.** ASCII normalises per primary at
  write time (single-run inspection); BDO writes raw sums plus the `nstat` tag
  so files merge correctly across runs.

## Full documentation

The detailed treatment — the accumulator / deposit / merge model and its role
in parallel scoring, the **output-normalisation and multi-run merging** rules
per format, unit handling and the post-process state table, and the module
lifecycle and boundaries — lives in the developer guide:

➡️ **[`docs/dev/scoring.md`](../../docs/dev/scoring.md)**

Parallelism groundwork and roadmap:
[issue #161](https://github.com/openshieldhit/openshieldhit/issues/161) and the
project [`TODO.md`](../../TODO.md).
