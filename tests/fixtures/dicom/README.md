Complete datasets with CT,RS and RN files can go here.
Keep this lean, < 500 MB. Only add what is absolutely needed for test and development.

These fixtures are stored in Git LFS. Local checkouts and CI jobs that run the
DICOM unit tests must fetch LFS payloads rather than plain pointer files.
