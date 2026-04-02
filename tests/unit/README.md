# Unit Tests

This directory contains C test programs for small, focused checks of library
behavior.

What belongs here:

- `test_*.c` source files
- tests for individual functions, parsers, helpers, and API behavior
- tests that run quickly and do not require a full OpenShieldHIT case

What does not belong here:

- full end-to-end input directories
- reference output datasets
- large collections of parser sample files

If a test needs input data, it should usually read it from `../fixtures/`.

