# Vendored MCPL core (issue #328)

`mcpl.h`, `mcpl.c`, `mcpl_fileutils.h`, `mcpl_fileutils.c` are the MCPL project's own
"core" reader/writer, vendored verbatim (Apache-2.0, see `LICENSE`) from
[mctools/mcpl](https://github.com/mctools/mcpl) at v2.2.8 — `mcpl.h` is generated from the
upstream `mcpl_core/include/mcpl.h.in` with its one template substitution
(`@MCPL_HOOK_FOR_ADDING_DEFINES@`) resolved to empty, since this build carries no extra
compile-time defines.

Compiled as plain source into a small static library (`src/thirdparty/mcpl/CMakeLists.txt`)
rather than pulled in via MCPL's own upstream CMake build: that build produces a shared
library plus CLI tools and tests, none of which `src/scoring/save/osh_scoring_save_mcpl.c`
needs, and compiling the four files directly avoids a shared-library deployment step.

The only external dependency this vendored core needs is zlib (`mcpl.c` includes `zlib.h`
unconditionally — this MCPL version has no build-time "no zlib" option, even though
OpenShieldHIT's own writer never calls a `gz*` function). `CMakeLists.txt` uses a system
zlib when `find_package(ZLIB)` succeeds, and otherwise fetches and builds one from source
(`madler/zlib`) so the build has no manual zlib-installation prerequisite on any of the
three CI platforms.

Not vendored: MCPL's own CLI tools (`mcpltool` etc.) and its SSW/PHITS converters
(`mcpl_extra/`) — neither is needed here, only the C API used directly by
`osh_scoring_save_mcpl.c`.
