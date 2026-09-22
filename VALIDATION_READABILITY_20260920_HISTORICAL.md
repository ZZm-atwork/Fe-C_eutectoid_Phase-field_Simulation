# Local verification of the documented single-file update

Completed 2026-09-20 Phoenix / 2026-09-21 UTC. These checks cover the local
build, geometry integration, legacy short-run regression and restart behavior.
They are not a new Sol validation or a long-run growth calibration.

## Code and input preservation

Both complete C++ files build with C++17 and the installed HDF5 library; the
parallel file uses the installed MPICH. Independent source review found zero
changes to numeric tokens, operators, string literals or statement order in
the identifier/comment-only refactor. Existing numerical kernels remain
unchanged by the separately integrated filling module.

All original nine example files retain their parsed key/value strings.
`Input.in` provides the annotated original G137 reference input. New geometry
commands are accepted only with `initialization = microsim`; their ordered
parameters enter a versioned checkpoint signature. Legacy signatures retain
their old values. The embedded inherited model digest is not a hash of the
updated file; actual source hashes are in `SHA256SUMS.txt`.

The first serial build exposed a naming-script error involving a period at the
end of a comment. The script now ignores comments/literals when distinguishing
member access. Failure logs were preserved; both corrected builds pass.

## Filling verification

Fourteen complete-grid fixtures compare all three phase fields against nine
byte-identical original MicroSim filling-function bodies. All fixtures match
exactly. Forty invalid/configuration cases reject correctly. The full test
performs 547,733 assertions, including repeated commands, explicit opt-in,
unsupported phases/3D input, required parameters and signature sensitivity.

G137 sharp labels contain alpha/theta/gamma counts of **2520/357/31373**.
The product seed occupies rows 0–20, with 120 alpha and 17 theta columns.
Native and explicit-cube initialization use the same 100 capillarity-only
pseudo-steps. Their saved `smooth_energy_converged=0` diagnostic is preserved;
no claim of fully converged smoothing is made.

## Saved-field comparisons: inherited thermodynamics

Each full short trajectory reaches **0.015 µs / 1,000 accepted steps**, with
4,121 total CG iterations and zero rejected attempts. The parallel run uses
four ranks, 2 × 2, with `FI_PROVIDER=tcp` for the local MPI transport.

| Comparison | Alpha max abs | Theta max abs | Gamma max abs | Carbon max abs | Mu max abs |
|---|---:|---:|---:|---:|---:|
| Updated native serial vs preserved old serial | 0 | 0 | 0 | 0 | 0 |
| Explicit cubes serial vs native serial | 0 | 0 | 0 | 0 | 0 |
| Explicit cubes MPI4 vs serial | 1.388e-16 | 0 | 1.665e-16 | 6.939e-18 | 1.819e-12 |

Maxima include t=0, endpoint/checkpoint and applicable sharp/smoothed phase
fields. Serial matches are bitwise exact. MPI uses the unchanged original
Fe-C smoke absolute tolerances: phase 1e-8, carbon 1e-10, mu 1e-5. All fields
are finite; maximum phase-sum error is 2.221e-16 and independently reconstructed
closure error is 5.552e-17. No tolerance was loosened.

## Actual restart checks

- A deliberately stopped 500-step serial cube run restarts on four ranks and
  reaches step 1,000. The loaded fields and persisted attributes match the
  split checkpoint bitwise. Compared with uninterrupted serial, endpoint phase
  errors are zero, carbon error is 3.469e-18 and mu error is 9.095e-13. Summed
  CG counts are 2,121 + 2,000 = 4,121, with zero rejections.
- The updated serial executable loads the frozen original v4 checkpoint with
  bitwise-identical initial fields and attributes, then advances 100 steps
  successfully. This demonstrates actual loading compatibility; a separate
  old-binary continuation endpoint was not rerun.

## Newly supplied thermodynamics

The separate input `examples/fe_c_new_thermo_reference_only.in` evaluates the
supplied ABC at **999.68459774 K**, not at 10 K undercooling. Both local serial
and MPI4 runs complete 0.015 µs, but need 5,475 accepted steps and 755 rejected
attempts. This adaptive activity is a **REVIEW item for long-run readiness**.
No claim of constant accepted dt, three-day stability or verified 10 K physics
follows. See `NEW_FREE_ENERGY.md` and the project calibration report.

Per-step accepted timesteps and numerical true-CG residuals are not persisted
by this executable. Saved CG counts, saved closure diagnostics and independent
closure reconstructions are different quantities, not residual histories.

## Evidence in the project workspace

- `research/pearl_v04_readability_20260920/`: before snapshot, source diffs,
  identifier/integration scripts, independent source review, build logs,
  geometry/Config tests, immutable local inputs and all local run provenance.
- `local_checks/independent_hdf_audit/` and
  `local_checks/independent_restart_audit/` beneath that directory: full
  field-by-field values, tolerances and input/source/executable/output hashes.
- `research/pearl_v04_new_thermo_20260920/`: supplied-data preservation,
  independent fits, temperature-law analysis, reference-only calibration and
  the pending 72-hour Sol resource plan.

The earlier mechanical-packaging report is retained as
`PACKAGING_VERIFICATION_20260914_HISTORICAL.json`. It describes the previous
source hashes and does not describe this new filling/readability update.
