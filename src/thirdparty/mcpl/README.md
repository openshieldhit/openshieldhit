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

## Licence files

`LICENSE` (Apache-2.0) and `NOTICE.md` are upstream's own, vendored alongside the
source: Apache-2.0 §4(a) and §4(d) require both to travel with any redistribution,
including binary ones, so the root `CMakeLists.txt` also installs them into
`share/openshieldhit/licenses/` for the `.deb`/`.tar.gz`/`.zip` packages. The
root `THIRD_PARTY_NOTICES.md` summarises all of this for a recipient who never
sees this directory.

Of the four source files only `mcpl.h` is modified (it is generated from
upstream's `mcpl.h.in` template, as described above); §4(b) requires that to be
stated, so the file carries a notice block at the top saying so. Re-vendoring a
newer MCPL means re-applying that block — nothing else here is patched.

## Building without MCPL

This directory is only added to the build when `OSH_ENABLE_MCPL` is `ON`, which
is the default. `-DOSH_ENABLE_MCPL=OFF` leaves every Apache-2.0 file here
uncompiled and unlinked; see the option's comment in the root `CMakeLists.txt`
for why that switch exists.
