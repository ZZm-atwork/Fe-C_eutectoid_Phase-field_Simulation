# Pearl: build and run guide

Pearl is a two-dimensional, three-phase Fe-C phase-field solver for pearlite
growth. The phase order is alpha ferrite (`0`), cementite (`1`), and austenite
(`2`). Both executables read the same plain-text input format and write HDF5/XMF
fields that can be opened in ParaView.

## 1. Package contents

| Path | Purpose |
|---|---|
| `pearl_serial.cpp` | Complete single-process solver |
| `pearl_parallel.cpp` | Complete MPI solver |
| `CMakeLists.txt` | CMake build configuration |
| `Input.in` | Annotated G137 input: 999.68459774 K ABC plus phase-diagram slopes, 10 K undercooling, 1.5 µs example horizon |
| `INPUT_PARAMETERS.md` | Parameter meaning, units, defaults and compatibility |
| `FILLING_METHODS.md` | Accepted MicroSim-style 2D filling recipes and unsupported historical methods |
| `CODE_MAP.md` | Successor navigation for the current combined CPP sections |
| `examples/fe_c_microsim_cube.in` | Opt-in explicit cubes reproducing the G137 sharp labels |
| `examples/fe_c_eutectoid.in` | Fe-C pearlite example using reference-temperature thermodynamics |
| `examples/fe_c_smoke.in` | Short Fe-C run for checking a build |
| `examples/fe_c_new_thermo_reference_only.in` | New supplied free energies at 999.68459774 K; historical reference-only check |
| `examples/fe_c_phase_diagram_10K.in` | Same new reference-plus-diagram input as `Input.in` |
| `examples/fe_c_phase_diagram_10K_smoke.in` | Same new thermodynamics with a 0.015 µs short-check horizon |
| `examples/symmetric_eutectic.in` | Small symmetric numerical example |
| `examples/legacy_direct/*.in` | Inputs with free-energy coefficients supplied directly at the run temperature |

The two C++ files contain numbered comment sections and physical variable names.
`CODE_MAP.md`, `INPUT_PARAMETERS.md`, and `FILLING_METHODS.md` describe the current
code and input format. `Solver_Guide.pdf` is preserved as the prior guide; its
old numbered blocks are not the authoritative map for this documented revision.

## 2. Requirements

- C++17 compiler
- CMake 3.16 or newer
- HDF5 C library
- MPI C++ implementation for `pearl_parallel`

## 3. Build

From the package directory, build the serial executable:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
```

This creates:

```text
build/pearl_serial
```

Build both serial and MPI executables:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPEARL_ENABLE_MPI=ON
cmake --build build -j 4
```

This creates:

```text
build/pearl_serial
build/pearl_parallel
```

If CMake needs explicit compiler or HDF5 locations:

```sh
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DPEARL_ENABLE_MPI=ON \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DHDF5_ROOT=/path/to/hdf5
cmake --build build -j 4
```

Direct compilation is also possible:

```sh
h5c++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic \
  pearl_serial.cpp -o pearl_serial

mpicxx -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic \
  $(pkg-config --cflags hdf5) pearl_parallel.cpp \
  -o pearl_parallel $(pkg-config --libs hdf5)
```

## 4. Input-file format

Inputs use one `key = value` entry per line. Text following `#` is treated as a
comment. Ordinary parameter keys cannot be duplicated. With
`initialization = microsim`, repeated `FILL...` commands are the sole exception:
they form an ordered geometry recipe. Unknown keys/methods stop the run; FILL
commands paired with another initialization are rejected instead of ignored.

For dimensional Fe-C inputs, the principal units are:

| Quantity | Unit |
|---|---|
| Physical length, `dx`, `dy`, `epsilon`, `seed_height`, `radius` | m |
| MicroSim-style FILL positions, widths, lengths and radii | Grid indices/grid cells; see below |
| FILL arm angles | Degrees counterclockwise from +X |
| Time, `dt`, `end_time`, `output_dt` | s |
| Temperature | K |
| Carbon concentration | atomic fraction |
| Phase free energy and chemical potential | J/mol |
| Diffusivity, `D_*` | m2/s |
| Interfacial energy, `sigma` | J/m2 |
| Molar volume, `Vm` | m3/mol |

Native scalar geometry keys such as `seed_height` and `radius` remain in meters.
Coordinates inside `FILL... = {...}` use global grid indices and lengths in
grid cells, matching the original MicroSim geometry convention. Most position
arguments are integers; arm/ellipse junctions also accept fractional grid
coordinates. A shape is
clipped at the domain edge; periodic images are not painted automatically.
The symmetric benchmark instead uses its consistent nondimensional units.

### 4.1 Thermodynamic input

Choose one thermodynamic mode.

#### Reference-temperature mode

Use this mode when the three phase free-energy parabolas are supplied at the
eutectoid temperature. There are two separately named temperature laws.

The current `Input.in` uses the supplied **999.68459774 K** free energies and
phase-diagram slopes. It does not require ABC exported at 989.68459774 K:

```text
thermo_mode = eutectoid_reference
T_eutectoid = 999.68459774
undercooling = 10
max_undercooling = 10
temperature_law = phase_diagram_fixed_curvature

Aeq_alpha = 4724862.6373624019
Beq_alpha = 22964.296428571855
Ceq_alpha = -42255.03701292308
# Supply reference ABC and D for theta and gamma as well.

ceq_alpha = 0.00088690859310575857
ceq_theta = 0.25
ceq_gamma = 0.034458832049305582
slope_alpha_on_alpha_gamma = -4.6282639461434345e-06
slope_gamma_on_alpha_gamma = -0.00031198492411304955
slope_theta_on_theta_gamma = 0
slope_gamma_on_theta_gamma = 0.00010753713136114516
```

These `slope_*` inputs are **dc_eq/dT in carbon mole fraction/K**. They are the
reciprocal of the older MicroSim `dT/dc` convention where finite; theta's fixed
stoichiometry is entered as zero. They are separate from `dq/dT` inputs.
The solver keeps the supplied curvatures fixed, holds gamma ABC as the common
affine energy-reference convention, and shifts the alpha/gamma and theta/gamma
contact compositions along their separate branches. Product B has a linear
temperature increment; product C retains its full quadratic increment. The
supplied reference tangent residual is preserved.

The run remains isothermal at `T = T_eutectoid - undercooling`, here
**989.68459774 K**. Diagram slopes fitted above the eutectoid are extrapolated
10 K below it; exact agreement with the third, alpha/theta coexistence branch
is not guaranteed. [NEW_FREE_ENERGY.md](NEW_FREE_ENERGY.md) documents the fit,
equations and limitations. This selected construction needs no lower-temperature
free-energy file.

The original `linear_gp` law remains available when independent temperature
derivatives of the grand-potential coefficients are supplied:

```text
thermo_mode = eutectoid_reference
T_eutectoid = 999.76
undercooling = 10.0
max_undercooling = 10.0
temperature_law = linear_gp

Aeq_alpha = ...
Beq_alpha = ...
Ceq_alpha = ...
dq2_dT_alpha = ...
dq1_dT_alpha = ...
dq0_dT_alpha = ...
D_alpha = ...
```

For `linear_gp`, provide `Aeq_*`, `Beq_*`, `Ceq_*`, `dq2_dT_*`, `dq1_dT_*`,
`dq0_dT_*`, and `D_*` for all phases. Do not mix these `dq*_dT` keys with the
diagram law's `ceq_*` and `slope_*` keys. The run
temperature is

```text
T = T_eutectoid - undercooling
```

The external `dq2_dT_*`, `dq1_dT_*`, `dq0_dT_*`, `Aeq_*`, and related keys retain
their original spellings. Inside the CPP, `Config::T`, `temperature_slope`,
`mu_squared_coefficient`, `mu_coefficient`, and `constant_energy` describe the
physical quantities. These descriptive C++ names are **not additional input
aliases**; `INPUT_PARAMETERS.md` gives the mapping.

#### Direct-temperature mode

Use this mode when each free-energy parabola has already been evaluated at the
desired run temperature:

```text
thermo_mode = direct_at_temperature
temperature = 989.76

A_alpha = ...
B_alpha = ...
C_alpha = ...
D_alpha = ...
```

Provide `A_*`, `B_*`, `C_*`, and `D_*` for all three phases.

### 4.2 Grid, interface, and kinetics

| Key | Meaning | Requirement/default |
|---|---|---|
| `nx`, `ny` | Number of cells in x and y | Required; `nx >= 3`, `ny >= 2` |
| `dx` | Cell spacing in x | Required |
| `dy` | Cell spacing in y | Default: `dx` |
| `Vm` | Common molar volume used in the chemical driving-force density | Required |
| `epsilon` | Diffuse-interface length parameter | Required |
| `sigma` | Common pairwise interfacial energy | Required |
| `triple` | Three-phase junction energy coefficient | Required |
| `tau_at` | Alpha-theta relaxation coefficient | Required |
| `tau_ag` | Alpha-gamma relaxation coefficient | Required |
| `tau_tg` | Theta-gamma relaxation coefficient | Required |

The spatial boundary conditions are periodic in x and zero normal flux in y.

### 4.3 Time integration and solver controls

| Key | Meaning | Requirement/default |
|---|---|---|
| `dt` | Initial time step | Required |
| `dt_max` | Maximum accepted time step | Default: `dt` |
| `dt_min` | Minimum trial time step | Default: `dt / 1048576` |
| `end_time` | Physical end time | Required |
| `output_dt` | Physical time between field outputs | Required |
| `max_steps` | Absolute accepted-step limit | Default: `100000000` |
| `max_phase_change` | Maximum phase-fraction change in one accepted step | Default: `0.04` |
| `cg_rtol` | Relative tolerance for the implicit diffusion solve | Default: `1e-11` |
| `cg_atol` | Absolute tolerance for the implicit diffusion solve | Default: `1e-14` |
| `cg_max` | Maximum conjugate-gradient iterations | Default: `2000` |
| `closure_tol` | Maximum allowed composition/chemical-potential closure error | Default: `1e-10` |
| `energy_atol` | Absolute tolerance in the free-energy acceptance test | Computed from the domain and `sigma` |
| `mode` | `coupled` or `phase_only` | Default: `coupled` |

`coupled` advances both phase fields and carbon diffusion. `phase_only` advances
the phase fields while maintaining thermodynamic composition closure without a
diffusion solve.

### 4.4 Initial phase geometry

Choose one initialization:

| `initialization` | Geometry | Related keys |
|---|---|---|
| `lamella` | Diffuse alpha/theta product seed below a horizontal gamma region | `seed_height`, `theta_columns` |
| `sharp_smooth` | Sharp alpha/theta/gamma labels followed by capillarity-controlled smoothing | `seed_height`, `theta_columns`, `smooth_*` |
| `flat` | Planar interface between two selected phases | `seed_height`, `phase_a`, `phase_b` |
| `circle` | Circular region of `phase_a` in `phase_b` | `radius`, `phase_a`, `phase_b` |
| `microsim` | Ordered sharp-label FILL commands followed by the existing capillary-only smoothing | One or more `FILL...` commands and `smooth_*` |

Phase selections use `0 = alpha`, `1 = theta`, and `2 = gamma`.

The accepted geometry keys are:

| Key | Default |
|---|---|
| `seed_height` | `ny * dy / 4` |
| `theta_columns` | Integer `nx / 2` |
| `radius` | One quarter of the shorter domain length |
| `phase_a` | `0` |
| `phase_b` | `2` |
| `smooth_steps` | `100` |
| `smooth_min_steps` | `smooth_steps` |
| `smooth_energy_rtol` | `1e-10` |
| `smooth_max_phase_change` | `0.02` |

For example, the G137 sharp seed can be written explicitly as:

```text
initialization = microsim
FILLCUBE = {1,0,0,0,16,20,0};
FILLCUBE = {0,17,0,0,136,20,0};
```

The commands use `phase,xlo,ylo,zlo,xhi,yhi,zhi`, with inclusive bounds and
`zlo=zhi=0`. They assign 17 theta columns and 120 alpha columns in rows 0–20.
Gamma occupies the unpainted background. The ready-to-read
`examples/fe_c_microsim_cube.in` retains the current G137 model/time controls and
uses this recipe; its other recipes are comments only. The original nine
examples keep their original parsed parameters and initialization. `Input.in`
now selects the new reference-plus-diagram thermodynamics while retaining the
G137 geometry and numerical controls.

Seven deterministic methods are available: `FILLCUBE`, `FILLCENTERBOX2D`,
`FILLCYLINDER`, `FILLYJUNCTION2D`, `FILLYJUNCTIONLAMELLAE2D`,
`FILLADHEREDPEARLITEARM2D`, and `FILLPEARLITEELLIPSE2D`. The last two preserve
alpha/theta overlay geometry on one gamma background. Three distinct gamma
grains, 3D, random grain variants and fluid velocity fields are unsupported.
See `FILLING_METHODS.md` for complete arguments, units, overwrite rules and
historical limitations. Do not activate multiple alternative recipes together.

### 4.5 Initial carbon state

Supply exactly one of the following:

| Key | Initialization rule |
|---|---|
| `initial_mu` | Applies the specified uniform carbon chemical potential |
| `gamma_carbon` | Computes the uniform chemical potential from the stated gamma composition |
| `mean_carbon` | Solves for a uniform chemical potential that gives the stated domain-average carbon fraction |

After the phase geometry is prepared, the local carbon concentration is
calculated from the phase mixture and the common initial chemical potential.

### 4.6 Moving window

The moving window follows an advancing transformation front in the y direction.

| Key | Meaning | Default |
|---|---|---|
| `moving_window` | `0` disables and `1` enables shifting | `0` |
| `shift_trigger` | Front position that initiates a shift | `0.65 * ny * dy` |
| `shift_target` | Front position after shifting | `0.35 * ny * dy` |
| `shift_gamma_carbon` | Carbon fraction inserted with new gamma rows | Derived from the initial gamma state |

Window shifts are evaluated at output times. Rows are removed from the bottom and
new gamma rows are inserted at the top.

## 5. Run the serial solver

The output directory must not already exist.

```sh
./build/pearl_serial Input.in run_serial
```

For a short build check:

```sh
./build/pearl_serial examples/fe_c_phase_diagram_10K_smoke.in smoke_serial
```

## 6. Run the MPI solver

The MPI rank count must equal `px * py`, and each MPI block must contain at least
one grid cell.

```sh
mpirun -np 16 ./build/pearl_parallel \
  Input.in run_parallel --px 4 --py 4
```

For Open MPI/HPC-X environments that require HCOLL to be disabled:

```sh
mpirun --bind-to none -mca coll_hcoll_enable 0 -np 16 \
  ./build/pearl_parallel Input.in run_parallel \
  --px 4 --py 4
```

This illustrative Slurm template uses the requested 72-hour allocation and
4 × 4 decomposition. It is not a submitted job. A production run needs its own
frozen input and physical-time horizon; the 1.5 µs package example is separate
from the 72-hour wall-time limit:

```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --partition=public
#SBATCH --job-name=pearl

# Use the site-compatible HPC-X/HDF5 environment before this command.
mpirun --bind-to none -mca coll_hcoll_enable 0 -np 16 \
  ./build/pearl_parallel Input.in run_parallel --px 4 --py 4
```

## 7. Command-line options

```text
pearl INPUT_FILE NEW_OUTPUT_DIR [options]
```

| Option | Use |
|---|---|
| `--px N --py N` | Set the MPI Cartesian decomposition |
| `--restart FILE.h5` | Continue from a checkpoint |
| `--max-steps N` | Stop after N additional accepted steps and write a checkpoint |
| `--thermo-only` | Evaluate and write the thermodynamic report without evolving the fields |

Example thermodynamic evaluation:

```sh
./build/pearl_serial Input.in thermo_report --thermo-only
```

## 8. Stop and restart

`Ctrl-C`, `SIGINT`, or `SIGTERM` requests a graceful stop. The solver completes
its current accepted step, writes the final frame and `checkpoint.h5`, and creates
`STOPPED.txt`.

To continue a run:

1. Copy the input file and set `end_time` to a value later than the checkpoint
   time.
2. Keep the physical model, grid, and moving-window settings compatible with the
   checkpoint. Keep the same filling mode and ordered FILL recipe for a `microsim`
   checkpoint; its versioned filling signature differs from legacy initialization
   even when the initial mask is identical.
3. Choose a new output directory.
4. Start with `--restart`.

Serial restart:

```sh
./build/pearl_serial continuation.in continuation_serial \
  --restart run_serial/checkpoint.h5
```

MPI restart:

```sh
mpirun -np 16 ./build/pearl_parallel continuation.in continuation_parallel \
  --restart run_parallel/checkpoint.h5 --px 4 --py 4
```

The MPI decomposition may be changed for a restart when the new rank count and
`px * py` agree.

## 9. Output files

| Output | Contents |
|---|---|
| `Input.in` | Input used for the run |
| `thermodynamics.txt` | Evaluated temperature-dependent thermodynamic coefficients |
| `initialization_report.txt` | Phase fractions and interface-smoothing summary |
| `init_sharp.h5/.xmf` | Sharp phase geometry for `sharp_smooth` or `microsim` initialization |
| `init_smoothed.h5/.xmf` | Smoothed phase geometry for `sharp_smooth` or `microsim` initialization |
| `frame_*.h5/.xmf` | Time-indexed alpha, theta, gamma, carbon, and chemical-potential fields |
| `diagnostics.csv` | Time step, carbon balance, free energy, solver, and wall-time diagnostics |
| `checkpoint.h5/.xmf` | Final restart state |
| `COMPLETED.txt` | End time reached |
| `STOPPED.txt` | Graceful stop requested |

Open the `.xmf` files in ParaView. The XMF file provides the grid geometry and
links the five scalar datasets stored in the matching HDF5 file.

## 10. Recommended run sequence

1. Start from `examples/fe_c_phase_diagram_10K_smoke.in` and confirm that the executable produces a
   complete output directory.
2. Use `--thermo-only` after changing thermodynamic parameters and inspect
   `thermodynamics.txt`.
3. Run the intended input for a limited number of steps with `--max-steps`.
4. Inspect phase fractions, carbon, chemical potential, and `diagnostics.csv`.
5. Submit the full serial or MPI calculation.

## Current local verification and new free energies

See [VALIDATION.md](VALIDATION.md) for the completed local geometry, field and
restart checks and the scope of the new temperature-law tests.
[NEW_FREE_ENERGY.md](NEW_FREE_ENERGY.md) describes the selected reference-plus-
phase-diagram construction. Runtime reports and frozen per-run inputs provide
actual job status and duration; the latest submission is documented below.

## Current supplied-energy run (2026-09-20 Phoenix)

The default `Input.in` now reads the supplied reference ABC at 999.68459774 K
with phase-diagram dc/dT slopes and 10 K undercooling. Its 1.5 µs horizon is a
usage example. `examples/fe_c_phase_diagram_200us.in` is the frozen production
input used for Sol job **63725141**, requested for 16 ranks/4×4, one CPU/rank,
HCOLL disabled and a 72-hour wall-time limit after short Sol job63725053 passed.
The physical target is 200 µs; wall-time, step-budget, numerical or top-boundary
guards may stop it earlier with a saved checkpoint. Raw data stay on Sol.

Read `VALIDATION.md`: local MPI4 retains a mu comparison failure, while the
intended Sol serial/MPI16/restart comparisons pass their original limits.
These results do not establish long-time convergence or matched phase velocities.
The frozen original v4 examples remain available under their original names.

