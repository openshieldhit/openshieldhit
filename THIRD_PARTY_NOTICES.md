# Third-party notices

OpenShieldHIT itself is distributed under the MIT License (see `LICENSE`).

A default build bundles third-party source code that carries its own, different
licence. Each entry below names the switch that leaves it out, and a build
configured that way contains none of that entry's code — this file still ships
so that a recipient can tell which case they have.
Those licences are **not** superseded by OpenShieldHIT's MIT licence: the
bundled files remain under the terms below, and this file exists so that every
recipient of an OpenShieldHIT distribution — source tarball, `.deb`, `.tar.gz`,
`.zip`, or wheel — receives the required licences and attribution notices.

---

## MCPL — Monte Carlo Particle Lists

- **Bundled at:** `src/thirdparty/mcpl/` (`mcpl.h`, `mcpl.c`,
  `mcpl_fileutils.h`, `mcpl_fileutils.c`)
- **Version:** 2.2.8
- **Upstream:** <https://github.com/mctools/mcpl> — <https://mctools.github.io/mcpl/>
- **Licence:** Apache License, Version 2.0
- **Full licence text:** `src/thirdparty/mcpl/LICENSE`, also installed as
  `share/openshieldhit/licenses/LICENSE.mcpl-apache-2.0`
- **Upstream NOTICE:** `src/thirdparty/mcpl/NOTICE.md`, also installed as
  `share/openshieldhit/licenses/NOTICE.mcpl`
- **Optional:** configure with `-DOSH_ENABLE_MCPL=OFF` to build OpenShieldHIT
  without this dependency; no Apache-2.0 code is then compiled or linked.

Copyright 2015-2026 MCPL developers.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
these files except in compliance with the License. You may obtain a copy of the
License at <http://www.apache.org/licenses/LICENSE-2.0>.

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

### Attribution (from upstream `NOTICE.md`)

> This software was mainly developed at the European Spallation Source ERIC
> (ESS) and the Technical University of Denmark (DTU). This work was supported
> in part by the European Union's Horizon 2020 research and innovation
> programme under grant agreement No 676548 (the BrightnESS project).
>
> A substantial effort went into developing MCPL. If you use it for your work,
> we would appreciate it if you would use the following reference in your work:
>
> T. Kittelmann, et al., Monte Carlo Particle Lists: MCPL, Computer Physics
> Communications 218, 17-42 (2017), <https://doi.org/10.1016/j.cpc.2017.04.012>

OpenShieldHIT passes that citation request on: if you use OpenShieldHIT's MCPL
input or output, please cite the paper above alongside OpenShieldHIT itself.
The reference is also recorded in `CITATION.cff`.

### Modifications (Apache-2.0 section 4(b))

`mcpl.c`, `mcpl_fileutils.c` and `mcpl_fileutils.h` are vendored **verbatim**.

`mcpl.h` is **not** upstream source verbatim: upstream ships it as the CMake
template `mcpl_core/include/mcpl.h.in`, and the bundled copy was generated from
that template with its single substitution `@MCPL_HOOK_FOR_ADDING_DEFINES@`
resolved to empty. OpenShieldHIT additionally prepended a comment block stating
this. The MCPL source below that block is unmodified. A matching notice is
carried in the file itself.

Not bundled at all: MCPL's CLI tools (`mcpltool`, `mcpl-config`), its Python
bindings, and the SSW/PHITS converters under `mcpl_extra/`.

### Licence-compatibility note for downstream users

Apache-2.0 is compatible with the MIT terms OpenShieldHIT's own code uses, and
with GPLv3, but it is **not** compatible with GPLv2-only. A project that must
link OpenShieldHIT into a GPLv2-only work should configure with
`-DOSH_ENABLE_MCPL=OFF`, which leaves the resulting binary free of Apache-2.0
code.

---

## zlib

MCPL's core includes `<zlib.h>` unconditionally, so a build with
`-DOSH_ENABLE_MCPL=ON` (the default) links zlib. OpenShieldHIT bundles **no**
zlib source: `src/thirdparty/mcpl/CMakeLists.txt` uses the system zlib when
`find_package(ZLIB)` finds one, and otherwise fetches `madler/zlib` v1.3.1 at
build time. In the latter case the resulting binary contains zlib code covered
by the zlib licence (<https://zlib.net/zlib_license.html>); see the fetched
source tree under your build directory for its full text.
