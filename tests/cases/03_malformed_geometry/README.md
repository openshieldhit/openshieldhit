Malformed geometry regression case.

This case is expected to fail during geometry parsing.
The geometry references a non-existent body name in the zone definition,
so the executable should terminate with a non-zero configuration error.

