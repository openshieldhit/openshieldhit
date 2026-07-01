# GitHub Copilot instructions

The agent brief for this repository lives in **[`../CLAUDE.md`](../CLAUDE.md)**.
Follow it (and the docs it links — [`DEVELOPER.md`](../DEVELOPER.md),
[`llms.txt`](../llms.txt)) when suggesting or reviewing code here.

The rules Copilot most often needs when generating C for this project:

- C11, but declare variables at the **top of each block** (K&R); none inside `for`/`if`.
- Use `struct foo` explicitly — **never** `typedef struct`.
- Block comments `/* */` only; `//` is reserved for temporary notes.
- `double const *p`, not `const double *p`.
- No non-portable APIs (POSIX-only or missing on MSVC): avoid `strcasecmp`/`<strings.h>`, `mkdtemp`, `getpid`, `<unistd.h>`, `<threads.h>`, … — Windows is a target.
- No heap allocation on the simulation hot path.
- Predicates return `int` (1 = yes); operations return `enum osh_status` (`OSH_OK = 0`).

See [`../CLAUDE.md`](../CLAUDE.md) for the full list and the build/test/format commands.
