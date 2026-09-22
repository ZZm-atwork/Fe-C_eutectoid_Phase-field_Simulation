# Import notes — 2026-09-21

Source archive: `v4_updated(2).zip`. SHA-256:
`933aee4ef0745e2f51d938f27d9cc0d0cea523fa6e51db0d705ffebd00ec19d2`.

Both C++ files, CMake configuration, input files, examples, and the supplied
MIT license are unchanged. The original README is preserved as `USER_GUIDE.md`.
A concise repository README, ignore rules, and text attributes were added.
The macOS metadata file was omitted. The workstation-specific absolute prefix
in one historical packaging JSON was replaced by `<original-workspace>`;
numeric results and filenames were retained.

The original SHA256SUMS is preserved as `IMPORTED_SHA256SUMS.txt`; it describes
the imported package before the documented packaging edits, including the old
README name. `SHA256SUMS.txt` covers the prepared repository files.

No new serial/MPI execution was performed: this environment lacked the required
MPI and HDF5 development dependencies. The imported `VALIDATION.md` records the
known local MPI4 discrepancy and the separately passing short Sol checks.
External research evidence referenced by those records was not in the archive.
