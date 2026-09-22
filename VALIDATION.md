# Current reference-temperature and phase-diagram validation

The current default input uses the supplied ABC at **999.68459774 K** and
phase-boundary slopes to construct a **989.68459774 K** (10 K undercooling) run.
It does not treat the supplied coefficients as lower-temperature data or reuse
the old grand-potential temperature derivatives. Read `NEW_FREE_ENERGY.md` for
the fixed-curvature approximation and `INPUT_PARAMETERS.md` for each input.

## Current results

The new temperature-law checks pass **171 assertions**. All 12 preserved legacy
input signatures match the preceding code exactly; all 11 existing example
files are byte-unchanged. Both combined C++ files compile locally. An actual
legacy-thermodynamics serial trajectory remains bitwise identical to the frozen
old result at t=0 and 0.015 microseconds. The numerical evolution kernels and
original numerical controls were not changed by this temperature-law addition.

The new 10 K local serial and MPI runs reach 0.015 microseconds with **5,455
accepted steps, 752 rejected attempts and 20,267 total CG iterations**. The
average accepted dt is 2.74977e-12 seconds; it is not a constant timestep.
One-rank MPI matches serial bitwise. The four-rank comparison is **FAIL** at the
unchanged chemical-potential limit:

| Field | Serial versus MPI4 maximum absolute error | Original limit | Result |
|---|---:|---:|---|
| alpha | 5.14890514e-9 | 1e-8 | PASS |
| theta | 5.14890530e-9 | 1e-8 | PASS |
| gamma | 5.85726606e-12 | 1e-8 | PASS |
| carbon | 7.02877756e-11 | 1e-10 | PASS |
| mu, J/mol | 0.0660804949 | 1e-5 | **FAIL** |

Initial fields match bitwise. Finite fields, phase closure, independent carbon
mass accounting and saved energy checks pass. A serial checkpoint at accepted
step 333 loads on four ranks with exact saved state; its continuation matches
the uninterrupted four-rank endpoint bitwise. This restart success does not
resolve the serial-versus-multiple-rank failure. Only six cells fail the mu
limit; high inverse mixture susceptibility amplifies small phase differences.
The origin of the perturbation is not yet established. No limit was relaxed.

Sol job **63725053** is the separate short platform verification: one node,
16 tasks, one CPU/task, 4x4 layout, HCOLL disabled, 15-minute allocation limit.
Its source/input/harness release is `phase_diagram_20260921_661051afa2cf`.
Job 63725053 completed with exit 0:0 in 205 seconds. All original Sol
comparisons and the 0.15 µs benchmark **PASS**; serial/MPI16 mu error is
9.09495e-11 J/mol and serial/restart mu error is 1.16415e-10 J/mol. The new-law
Sol runs each use 5,464 accepted steps, 753 rejections and 20,373 CG iterations
at 0.015 µs. The intended Sol environment is validated for this short scope;
the local MPI4 result remains separately qualified.
Production job **63725141** is running on sc067: 72-hour allocation limit, 200 µs
physical target, the same validated MPI executable, and unchanged remaining
input parameters. The local MPI4 failure remains recorded; this is a qualified
Sol validation, not a claim that all local comparisons passed.
See `PACKAGING_VERIFICATION.json` and the current project report for final status.

## Scope and preserved evidence

Per-step accepted-dt history and numerical true-CG residual history are not
persisted. Saved CG counts, the source-level true-residual convergence test,
saved closure diagnostics and reconstructed closure are distinct evidence;
none is a saved numerical true-residual trajectory. No long-time stability,
cooperative-growth, physical calibration or timestep-convergence claim follows
from these short tests.

The earlier readability/filling checks, including 547,733 assertions and 14
whole-grid original-MicroSim fixtures, are retained in
`VALIDATION_READABILITY_20260920_HISTORICAL.md`. That file describes the previous
revision and reference-only calibration, not the new 10 K result. Its original
source and all frozen v3/v4 data are preserved.

Current evidence: `../research/pearl_v04_phase_diagram_20260921/REPORT.md`,
`LOCAL_FIELD_VALIDATION.json`, `MPI1_DIAGNOSTIC.json`, `tests/TEST_RESULTS.json`,
`DIAGNOSIS_STATIC.json`, `DOCUMENTATION_CHECK.json` and the Sol evidence there.
