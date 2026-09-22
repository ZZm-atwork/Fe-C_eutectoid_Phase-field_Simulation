# Pearl: Fe–C pearlite phase-field solver

A two-dimensional, three-phase C++17 research solver for ferrite, cementite,
and austenite evolution. The current v0.4.0 source provides serial and MPI
executables, plain-text inputs, HDF5/XMF output, and checkpoint/restart support.

## Features

- Three-phase grand-potential formulation with carbon transport.
- Direct-temperature and eutectoid-reference thermodynamic inputs.
- Fixed-curvature phase-diagram temperature law in the current default input.
- MPI Cartesian decomposition in two dimensions.
- Adaptive time stepping, implicit diffusion, and acceptance diagnostics.
- MicroSim-style deterministic geometry initialization and moving-window support.

## Build

Requirements: a C++17 compiler, CMake 3.16+, HDF5 C development libraries, and an
MPI development environment for the parallel executable.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPEARL_ENABLE_MPI=ON
cmake --build build -j 4
```

Omit `-DPEARL_ENABLE_MPI=ON` for a serial-only build.

## Quick start

Use a new output directory for each run:

```sh
./build/pearl_serial examples/fe_c_phase_diagram_10K_smoke.in smoke_serial
mpirun -np 4 ./build/pearl_parallel \
  examples/fe_c_phase_diagram_10K_smoke.in smoke_mpi --px 2 --py 2
```

The smoke input targets 0.015 microseconds. It is a numerical startup check,
not a production growth simulation. Open the generated XMF files in ParaView.

## Documentation

| File | Contents |
| --- | --- |
| [USER_GUIDE.md](USER_GUIDE.md) | Complete build, execution, restart, and output guide |
| [INPUT_PARAMETERS.md](INPUT_PARAMETERS.md) | Inputs, units, and defaults |
| [FILLING_METHODS.md](FILLING_METHODS.md) | Geometry recipes |
| [NEW_FREE_ENERGY.md](NEW_FREE_ENERGY.md) | Thermodynamic construction and limitations |
| [CODE_MAP.md](CODE_MAP.md) | Numbered source sections |
| [VALIDATION.md](VALIDATION.md) | Qualified current validation record |
| [IMPORT_NOTES.md](IMPORT_NOTES.md) | Snapshot provenance and preparation checks |

## Validation status

The imported validation record reports successful short serial/MPI16/restart
comparisons on ASU Sol. It also records an unresolved local MPI4 comparison:
maximum chemical-potential error 0.0660804949 J/mol versus a 1e-5 J/mol limit.
That failure remains unresolved in this repository. Read [VALIDATION.md](VALIDATION.md)
for the scope and platform-specific results; short checks do not establish
long-time stability or physical calibration.

The repository preparation environment lacked MPI and HDF5 development
libraries, so no new Pearl compile or execution result is claimed here.
The original validation documents reference external research evidence that
was not included in the uploaded package. Their job-status statements are
historical snapshot statements, not live cluster monitoring.

## Version and license

This import uses `v4_updated(2).zip`, uploaded on 2026-09-21, and preserves both
solver CPP files byte-for-byte. It is the current v4 snapshot, not a synthetic
v3-to-v4 commit history. The supplied [MIT license](LICENSE) is retained.
