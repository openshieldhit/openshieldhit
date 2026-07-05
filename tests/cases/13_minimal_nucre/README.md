# 13_minimal_nucre

A minimal nuclear-reactions-on (`NUCRE 1`) case: 150 MeV protons in water, the
same geometry/material as `00_minimal`.  Enabling nuclear reactions routes ion
secondaries, sub-threshold fragment **point deposits**, and neutron scoring
through the deposit-target seam, so this is the nuclear counterpart used by the
`--score-replicas` reproducibility tests (issue #230) alongside the EM-only
`00_minimal`.

The bare `cases::13_minimal_nucre` entry runs it in `--dry-run` (parse only).
The real-transport reproducibility checks are registered separately by
`run_score_replicas.cmake` (see `tests/cases/CMakeLists.txt`).
