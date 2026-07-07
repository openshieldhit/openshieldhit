# Bug hunt reports

A **bug hunt** is a deep, point-in-time audit of the codebase, looking beyond
what CI and the normal review process catch: physics/numerical correctness,
silent result bias, parallelization readiness, and architecture drift. Each
report is a snapshot pinned to one commit, not a living document — if a
finding gets fixed, the fix lands as a normal commit/issue, and the report
stays as a historical record of what was found and when.

## How to read a report

Findings are grouped per subsystem and given a stable ID (e.g. `G-1`, `N-2`,
`S-1`) plus two tags:

| Tag | Meaning |
|-----|---------|
| **Severity** | critical / high / medium / low / info — expected impact if the finding is real |
| **Confidence** | high = verified by running code or by construction; medium = strong code reading; low = suspicion worth checking |

Each report opens with a metadata block (commit audited, source PR, scope)
and closes with an executive summary and a suggested issue breakdown.

!!! warning "Provenance and how much to trust these findings"
    Bug hunt reports so far have been produced by an AI coding agent working
    autonomously from the source tree, not by a maintainer's line-by-line
    review. Every finding cites exact `file:line` evidence and, where a
    numeric physics claim is made, is checked against a running build or a
    standalone harness — but none of it is a substitute for maintainer
    review. Treat findings as **triage input**: high-confidence ones are
    usually worth filing as issues directly; medium/low-confidence ones are
    worth a second look before anyone spends time fixing them.

Follow-up issues get filed from a report's suggested breakdown over time;
when one is filed, cross-link it from the finding's section (or the report's
metadata block) so the two stay connected — check there before opening a
duplicate.

## Reports

| Date | Commit audited | Scope | Headline findings | Report |
|------|-----------------|-------|--------------------|--------|
| 2026-07-05 | [`e3c619a`](https://github.com/openshieldhit/openshieldhit/commit/e3c619a1328e9351bcbc1dc599321ac2770ad622) | Test infra, transport/EM physics, scoring, nuclear/neutron transport, RNG & parallel readiness, GEMCA geometry, parsers/IO, architecture | 1 critical, 5 high | [Deep audit →](2026-07-05-deep-audit.md) |

## Adding a new report

1. Add the report as its own page here, named `YYYY-MM-DD-<slug>.md`.
2. Add a row to the table above.
3. Register the page in `mkdocs.yml` under **Developer Guide → Bug hunts**.
4. Link every referenced issue/PR number as a Markdown reference link
   (`[#NNN]` in prose, with a `[#NNN]: https://github.com/openshieldhit/openshieldhit/issues/NNN`
   — or `/pull/NNN` — definition collected at the bottom of the file). Bare
   `#NNN` only autolinks on GitHub itself, not on this site.
5. Give every individual finding heading a stable anchor
   (`` {: #g-1 } `` after a `### G-1 (severity, confidence) ...` heading) so
   other reports, issues, and PRs can link straight to it.
