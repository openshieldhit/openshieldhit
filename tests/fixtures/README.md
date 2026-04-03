# Test Fixtures

This directory contains small reusable input files for unit tests.

What belongs here:

- tiny parser-focused sample inputs
- valid and invalid data fragments
- shared test data that multiple unit tests may read

What does not belong here:

- standalone C test programs
- full regression cases intended to be run like normal OpenShieldHIT jobs
- large reference output datasets

Use fixtures when a unit test needs controlled sample input without depending
on a full case directory from `../cases/`.

