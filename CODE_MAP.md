# Reading the combined Pearl CPP files

`pearl_serial.cpp` and `pearl_parallel.cpp` contain the same physical model and
input format. Each file is self-contained. The preamble selects the existing
serial or MPI implementation; the parallel path uses collective reductions and
halo exchange. Do not edit one physical implementation without applying the
same change to the other.

Search for the numbered section titles or function names below. Function names
are more stable navigation targets than line numbers. `Solver_Guide.pdf` is the
preserved earlier guide; this file and the current CPP comments describe the
updated section layout.

## Physical names

| CPP name | Physical or numerical meaning |
|---|---|
| `Config::T` | Constant run temperature in K |
| `Config::temperature_slope` | Per-phase derivatives of grand-potential coefficients for the existing `linear_gp` law |
| `PhaseDiagramInput::reference_composition` | Alpha, theta and gamma contact compositions at the supplied reference temperature |
| `PhaseDiagramInput::composition_slope` | Four dc_eq/dT slopes: alpha-on-AG, gamma-on-AG, theta-on-TG, gamma-on-TG |
| `free_energy_quadratic`, `free_energy_linear`, `free_energy_constant` | A, B and C in `f(carbon)=A*carbon²+B*carbon+C` |
| `mu_squared_coefficient`, `mu_coefficient`, `constant_energy` | Coefficients of `mu²`, `mu`, and 1 in the grand-potential polynomial |
| `State::phi`, `State::carbon`, `State::mu` | Phase fractions, carbon atomic fraction, and chemical potential |
| `mu_reference` | Reference chemical potential subtracted in the implicit diffusion solve |
| `chi` | Susceptibility `dc/dmu` |
| `K` | Diffusive mobility `sum(phi*D*chi)`, not diffusivity alone |
| `residual_vector`, `preconditioned_residual`, `search_direction` | Conjugate-gradient linear-solver vectors |
| `cg_step`, `direction_weight` | CG update coefficients; they are unrelated to alpha/theta phase IDs |

Existing input keys retain their original meanings. `dq2_dT_*`, `dq1_dT_*`,
and `dq0_dT_*` belong to `linear_gp` and are not phase-boundary slopes. The new
`phase_diagram_fixed_curvature` law separately accepts `ceq_*` and four
`slope_*_on_*` keys. These mean dc_eq/dT in mole fraction/K, the reciprocal of
the older MicroSim dT/dc where finite; theta's zero dc/dT is represented directly.
Neither family is a cooling rate or spatial temperature gradient. Descriptive
C++ names are not parser aliases.

## Section navigation

| Section in both CPPs | Main types/functions | What to inspect here |
|---|---|---|
| **1. THERMODYNAMICS** | `Parabola`, `GrandPotentialPolynomial`, `PhaseDiagramInput`, `mixture`, `composition`, `to_grand_potential`, `from_grand_potential`, `at_undercooling`, `at_phase_diagram_undercooling` | Both reference-to-run-temperature laws, shared chemical potential, phase carbon response, grand-potential driving force |
| **2. REFERENCE DIFFUSION** | `Grid`, `Field`, `solve_cg`, `diffuse` | Retained serial reference diffusion implementation; production distributed evolution is in section 6 |
| **3. CAPILLARITY** | `InterfaceParameters`, `capillary_energy_gradient`, `project_simplex` | Interface energy/gradient and nonnegative phase fractions summing to one |
| **4. PHASE KINETICS** | `PhaseParameters`, `phase_step` | Pair relaxation coefficients and explicit constrained phase updates |
| **5. SPATIAL GRID** | `Array`, `Domain`, `halo`, `gather`, `load_global` | Global/local indexing, ghost layers, MPI decomposition/reductions, serial equivalents |
| **6. COUPLED EVOLUTION** | `State`, `Controls`, `smooth_interface_only`, `propose_phi`, `implicit_diffusion`, `valid_state`, `advance` | Initialization preprocessing, physical coupled steps, diffusion solve, closure and accept/reject logic |
| **FILLING GEOMETRY** | `FillingOperation`, `parse_filling_operation`, `validate_filling_operations`, `filling_phase_at` | Ordered optional 2D MicroSim-style sharp-label masks, strict argument checks and overwrite semantics |
| **7. INPUT PARAMETERS AND FILLING** | `Config`, `physics_signature`, `read_config`, `initialize` | Read/validate parameters, choose geometry, smooth when requested, then initialize common-mu carbon |
| **8. OUTPUT AND RESTART** | `write_frame`, `read_checkpoint`, `write_initial_phases`, `shift_rows`, `write_thermodynamic_report` | HDF5/XMF output, provenance, restart validation and moving-window mass accounting |
| **9. MAIN PROGRAM** | `main` | Command-line options, fresh output directory, initialize/restart decision, output cadence, final checkpoint and completion marker |

## Data flow

1. `read_config` reads ordinary parameters and ordered FILL commands separately.
   It evaluates the fixed run-temperature thermodynamics and validates inputs.
2. On a new run, `initialize` creates phase geometry. Existing native modes retain
   their old behavior. `initialization=microsim` evaluates sharp labels at global
   integer indices; every MPI rank uses the same geometry definition.
3. `sharp_smooth` and `microsim` call the same capillarity-only preprocessing.
   This consumes no physical time and performs no carbon diffusion.
4. The initialized phase mixture and composition choice set a uniform chemical
   potential and its corresponding carbon field. A restart instead loads the
   saved fields and bookkeeping after checking the full physics signature.
5. `advance` proposes the phase change, solves conservative diffusion, verifies
   the coupled state and accepts or rejects the complete step. A rejection halves
   trial `dt`; output and end times can also shorten an accepted step.
6. Output and optional moving-window operations preserve the carbon ledger and
   global displacement needed for later analysis.

## Compatibility boundaries

The readability revision preserved existing equations and ordinary input
spellings. The subsequent phase-diagram extension adds a separately selected
temperature law; it does not reinterpret the historical `linear_gp` law or
direct-at-temperature mode. The original nine examples retain their parsed
values. `Input.in` now explicitly selects the supplied 999.68459774 K ABC plus
phase-diagram slopes with 10 K undercooling, while retaining the inherited G137
geometry and numerical controls.

`at_phase_diagram_undercooling` fixes reference curvatures and gamma ABC, moves
the alpha/gamma and theta/gamma contacts by their separate composition slopes,
and adds the resulting B/C increments to the supplied reference coefficients.
Its C increment is quadratic in temperature; it preserves reference ABC exactly
at zero undercooling and retains the supplied tangent residual. Gamma B/C is a
common affine reference choice; fixed curvature is a substantive model assumption.
Stable gamma-containing branches above Te are extrapolated below Te, and the
third alpha/theta boundary is a diagnostic rather than an independently fitted
constraint. See `NEW_FREE_ENERGY.md` for the formula and evidence.

The optional filling feature adds sharp geometry recipes. It does not add gamma
grain fields, three-dimensional evolution, or hydrodynamics. The default/native
initialization path remains available.

Legacy input signatures retain their existing layout and inherited model digest.
`initialization=microsim` appends `filling-v1` and the ordered method/coordinate
sequence. Keep that recipe when restarting. The model digest is not the current
CPP file hash; current file hashes and revision evidence are recorded separately.
The phase-diagram law records `phase-diagram-fixed-curvature-v1`, the three
reference contacts and all four slopes in its physics signature. Restarting
with another temperature law or changed contacts/slopes is not silently accepted.

The true CG residual is recomputed internally for convergence. This code does
not persist a numerical true-residual history or every accepted timestep.
Output-time diagnostics and aggregate counts must not be presented as those
missing histories.
