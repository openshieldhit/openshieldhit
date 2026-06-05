# Input files overview

Every simulation case is a directory with four plain-text input files.

| File | Purpose |
|------|---------|
| `beam.dat` | Primary particle, energy, beam geometry, physics switches, statistics |
| `geo.dat` | Geometry: bodies, zones, or DICOM CT |
| `mat.dat` | Material compositions and cross-section data |
| `detect.dat` | Scoring geometries and quantities |

All four files must be present.  openshieldhit reads them in the order listed; none
takes precedence over another — they are independent input domains.

## File format conventions

- Plain ASCII text, one keyword per line
- Arguments follow the keyword on the same line, whitespace-separated
- Lines beginning with `#` or `*` are comments (inline comments with `#` also
  work on keyword lines)
- Keywords are **case-insensitive**
- Unknown keywords produce a parse error (no silent ignoring)

## SH12A compatibility

openshieldhit uses the same file format as SHIELD-HIT12A.  Existing SH12A input
directories can generally be used without modification, with two caveats:

1. **USECBEAM coordinate convention** — spot X/Y positions are always
   isocenter coordinates in both SH12A and openshieldhit.  See [beam.dat](beam.dat.md#usecbeam).
2. Some rarely-used SH12A keys are not yet implemented; openshieldhit will error on them.
   Known stubs are listed in the respective reference pages.
