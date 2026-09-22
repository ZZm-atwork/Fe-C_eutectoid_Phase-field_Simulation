# Pearl v0.4 input parameters and filling

`Input.in` is the annotated G137 input using the supplied **999.68459774 K**
free-energy coefficients and four phase-diagram slopes, with 10 K undercooling.
It matches `examples/fe_c_phase_diagram_10K.in`: 137 × 250 cells, a 17-cell theta
seed, and an absolute stopping time of 1.5 µs. The short
`examples/fe_c_phase_diagram_10K_smoke.in` instead stops at 0.015 µs. These example
physical times are separate from a Slurm wall-time limit or production horizon.
The nine original example files retain every original parsed key and value.
`examples/fe_c_microsim_cube.in` is a separate
opt-in example: it changes only the filling selection to `microsim` and adds
the two G137 cube commands; the physical and numerical parameters are unchanged. Comments and grouping make the inputs easier to
read. The new default thermodynamics is an explicit scientific input change;
comments by themselves do not authorize running an example.

The input has two clearly marked parts: **INPUT PARAMETERS** for the model,
grid, time integration, and output; **FILLING** for the initial phase geometry,
carbon state, and optional interface preprocessing. Filling methods inherited
from MicroSim are documented separately in [FILLING_METHODS.md](FILLING_METHODS.md).
The tables below describe the original v0.4 parameters, the new separately named
phase-diagram law, native filling, and the opt-in `microsim` command recipe.
Existing input names and their meanings remain stable.

## Format, phases, and units

Each active line is `key = value`. A `#` starts a comment, including after a
value. Do not add section headings without `#`, unit suffixes to numeric values,
or duplicate ordinary keys. With `initialization=microsim`, repeated `FILL...`
commands are an ordered list and are allowed; no other key may repeat. Each
command must fit on one line. A final semicolon is optional for FILL lists only.
Names are case sensitive. Unknown keys/methods are errors. Ordinary
and scientific-notation numbers are accepted; nonfinite numbers are rejected.
Ordinary integer fields require exact integer values and magnitude no greater
than 100,000,000. FILL-list indices/counts use method-specific integer and domain
checks described in `FILLING_METHODS.md`; nonfinite or unresolved derived
geometry is rejected. Unless noted otherwise, an omitted required field is an error.

Phase IDs are **0 = alpha/ferrite**, **1 = theta/cementite**, and
**2 = gamma/austenite**. This is a three-phase model; the older MicroSim gamma
clone IDs are not additional phases in Pearl.

For the dimensional Fe-C examples, lengths are in m, times in s, carbon `c` is
atomic fraction rather than wt%, and chemical potential `mu = df/dc` and the
molar-code free energy `f` are in J/mol. `Vm` converts the chemical contribution
to J/m³. The 2D integrated energy is per unit out-of-plane thickness, in J/m.
These are the inherited code conventions, not a fresh certification of the
material-property normalization. `symmetric_eutectic.in` instead uses consistent
nondimensional benchmark units; its numbers are not Fe-C material data.

**FILL-list geometry is the units exception:** positions are global grid indices,
lengths/widths/radii are in grid cells, angles are in degrees counterclockwise
from +X, and fractions/aspect ratios are dimensionless. Native `seed_height` and
`radius` remain physical lengths in meters. Cube high bounds are inclusive;
center-box and Y bounding boxes have exclusive high bounds. See the method
catalog before converting an old filling file.

## Input parameters: temperature and thermodynamics

The solver evaluates one run temperature before initialization, then runs
isothermally. There is no cooling-rate input or spatial temperature-gradient
input. Changing temperature does not automatically refit diffusivities or phase
relaxation coefficients.

| Key | Meaning and units | Requirement/default |
|---|---|---|
| `thermo_mode` | Selects the interpretation of the supplied free-energy coefficients | Required: `eutectoid_reference` or `direct_at_temperature` |
| `T_eutectoid` | Reference eutectoid temperature, K | Required in reference mode; positive |
| `undercooling` | Positive reduction from the reference temperature, K | Required in reference mode; between 0 and `max_undercooling` |
| `max_undercooling` | Declared applicability limit of the supplied temperature law, K | Required in reference mode; nonnegative |
| `temperature_law` | Reference-to-run-temperature construction | Required in reference mode: `linear_gp` or `phase_diagram_fixed_curvature` |
| `temperature` | Run temperature, K; no additional correction is applied to supplied ABC | Required in direct mode; positive |

In reference mode, `T = T_eutectoid - undercooling` must be positive. For each
phase name `p` in `alpha`, `theta`, `gamma`, supply the following:

| Keys for phase `p` | Meaning | Units and requirement |
|---|---|---|
| `Aeq_p`, `Beq_p`, `Ceq_p` | Coefficients of `f(c,T_eutectoid) = Aeq*c² + Beq*c + Ceq` | J/mol; all required in reference mode; `Aeq_p > 0` |
| `dq2_dT_p` | Temperature slope of the coefficient multiplying `mu²` in the grand potential | mol/(J·K); required only for `linear_gp` |
| `dq1_dT_p` | Temperature slope of the coefficient multiplying `mu` in the grand potential | K⁻¹; required only for `linear_gp` |
| `dq0_dT_p` | Temperature slope of the constant grand-potential term | J/(mol·K); required only for `linear_gp` |
| `A_p`, `B_p`, `C_p` | Coefficients of `f(c,T) = A*c² + B*c + C` already at the run temperature | J/mol; all required in direct mode; `A_p > 0` |
| `D_p` | Phase carbon diffusivity; transport uses `K = sum(phi_p * D_p * dc_p/dmu)` | m²/s; required in both modes; nonnegative |

For example, `Aeq_p` means `Aeq_alpha`, `Aeq_theta`, or `Aeq_gamma`; the literal
suffix `_p` is not an accepted key. Only the coefficient family belonging to the
selected mode and law is accepted. Do not supply `A_alpha` and `Aeq_alpha`
together, and do not combine `dq*_dT` with the phase-diagram keys below.

### Reference ABC with phase-diagram slopes

`temperature_law = phase_diagram_fixed_curvature` is the current default input's
selected construction. Supply reference `Aeq_*`, `Beq_*`, `Ceq_*`, all `D_*`,
and these seven additional keys:

| Key | Meaning | Units and requirement |
|---|---|---|
| `ceq_alpha` | Alpha contact of the supplied reference alpha/gamma free-energy curves | Carbon mole fraction; required, between 0 and 1 |
| `ceq_theta` | Theta reference contact, 0.25 for the supplied surrogate | Carbon mole fraction; required, between 0 and 1 |
| `ceq_gamma` | Shared gamma reference contact of the supplied reference curves | Carbon mole fraction; required, between 0 and 1 |
| `slope_alpha_on_alpha_gamma` | Alpha-side `dc_eq/dT` on the alpha/gamma boundary | Mole fraction/K; required and finite |
| `slope_gamma_on_alpha_gamma` | Gamma-side `dc_eq/dT` on the alpha/gamma boundary | Mole fraction/K; required and finite |
| `slope_theta_on_theta_gamma` | Theta-side `dc_eq/dT` on the theta/gamma boundary | Mole fraction/K; required and finite; zero is valid |
| `slope_gamma_on_theta_gamma` | Gamma-side `dc_eq/dT` on the theta/gamma boundary | Mole fraction/K; required and finite |

The two gamma-side slopes belong to different coexistence branches. They are
not interchangeable. These inputs use **dc/dT**, while the older MicroSim
`slopes` table used **dT/dc**; take the reciprocal only when finite and after
converting mole percent to mole fraction. A stoichiometric theta branch uses
`slope_theta_on_theta_gamma = 0`, avoiding an infinite dT/dc value. The values
are composition slopes, not derivatives of grand-potential coefficients,
cooling rates, or spatial temperature gradients.

For product phase `p` (alpha or theta), the helper
`at_phase_diagram_undercooling` evaluates:

```text
delta_T = -undercooling
dp = slope_p_on_p_gamma * delta_T
dg = slope_gamma_on_p_gamma * delta_T
A_p(T) = Aeq_p
B_p(T) = Beq_p + 2*(Aeq_gamma*dg - Aeq_p*dp)
C_p(T) = Ceq_p + Aeq_p*dp*(2*ceq_p+dp)
                    - Aeq_gamma*dg*(2*ceq_gamma+dg)
```

Gamma ABC remains equal to its reference values; fixed gamma B/C chooses the
common affine energy-reference convention. Fixing every A is a separate model
assumption. Product C includes the **full quadratic** increment in delta_T.
This is the reference-anchored form of the MicroSim pair-boundary construction,
not a linear-C truncation or an implicit reinterpretation of `linear_gp`.
At zero undercooling the supplied ABC values are returned unchanged.

All three `ceq_*` and four extrapolated pair-contact compositions must be finite
and in [0,1]. The source preserves any supplied reference tangent discrepancy
instead of modifying ABC to remove it. The new input uses the actual fitted
alpha/gamma contacts and theta=0.25; their small differences from the phase
diagram's invariant compositions are documented in `NEW_FREE_ENERGY.md`.

The supplied slopes use an anchored fit to unique diagram points between Te
and Te+10 K. Their use below Te is a metastable extrapolation. The independently
predicted alpha/theta branch generally differs from the exported third boundary,
especially because the finite-curvature theta surrogate is not an exactly
stoichiometric thermodynamic phase. No automatic fit or slope estimation occurs
inside the solver.

### Existing linear grand-potential temperature law

The input names `dq2_dT_*`, `dq1_dT_*`, and `dq0_dT_*` are retained for existing
input compatibility. Here **q is a grand-potential coefficient; dq/dT is its
temperature slope**. It is not a phase-diagram tie-line slope:

```text
Psi(mu,T) = q2(T)*mu² + q1(T)*mu + q0(T)
q2 = -1/(4*A)
q1 = B/(2*A)
q0 = C - B²/(4*A)
qk(T) = qk(T_eutectoid) - undercooling * dqk_dT
```

The minus sign follows from `T - T_eutectoid = -undercooling`: the slopes are
derivatives with respect to **temperature**, not with respect to undercooling.
The evaluated `q2` must remain negative. At zero undercooling the original
reference ABC values are retained exactly. A single set of reference ABC values
does not determine temperature slopes; those are separate material inputs.

### External names and internal physical names

The input parser retains the original names and adds explicit keys for the
separate phase-diagram law. Descriptive C++ identifiers are not parser aliases.
For phase `p`, the mapping is:

| External input name | CPP quantity | Meaning |
|---|---|---|
| `temperature`, or `T_eutectoid-undercooling` | `Config::T` | Evaluated constant run temperature |
| `dq2_dT_p` | `temperature_slope[p].mu_squared_coefficient` | Temperature derivative of the `mu²` coefficient |
| `dq1_dT_p` | `temperature_slope[p].mu_coefficient` | Temperature derivative of the `mu` coefficient |
| `dq0_dT_p` | `temperature_slope[p].constant_energy` | Temperature derivative of the constant grand-potential term |
| `ceq_alpha`, `ceq_theta`, `ceq_gamma` | `phase_diagram.reference_composition` | Reference contact compositions in phase order |
| Four `slope_*_on_*` keys above | `phase_diagram.composition_slope` | Alpha-on-AG, gamma-on-AG, theta-on-TG, gamma-on-TG composition slopes |
| `A_p` / `Aeq_p` | `free_energy_quadratic` | Coefficient multiplying carbon² |
| `B_p` / `Beq_p` | `free_energy_linear` | Coefficient multiplying carbon |
| `C_p` / `Ceq_p` | `free_energy_constant` | Free-energy constant term |
| `D_p` | `diffusivity` | Phase carbon diffusivity |

`State::carbon` is the carbon field. `mu_reference` names the reference chemical
potential used by the shifted diffusion solve. CG names such as
`residual_vector`, `search_direction`, and `cg_step` describe linear-solver
quantities, not ferrite/cementite phases. `CODE_MAP.md` locates these sections.

## Input parameters: grid, interfacial energy, and kinetics

| Key | Meaning and units | Requirement/default |
|---|---|---|
| `nx` | Number of physical cells in x | Required; at least 3 |
| `ny` | Number of physical cells in y | Required; at least 2 |
| `dx` | Cell spacing in x, m | Required; positive |
| `dy` | Cell spacing in y, m | Defaults to `dx`; positive |
| `Vm` | Common molar volume converting molar chemical energy to energy density, m³/mol | Required; positive |
| `epsilon` | Diffuse-interface model length, m; not the measured 10–90% width | Required; positive |
| `sigma` | Equal pairwise interfacial energy, J/m² | Required; positive |
| `triple` | Three-phase overlap energy coefficient, J/m² | Required; nonnegative |
| `tau_at` | Alpha-theta phase relaxation coefficient, J·s/m⁴ | Required; positive |
| `tau_ag` | Alpha-gamma phase relaxation coefficient, J·s/m⁴ | Required; positive |
| `tau_tg` | Theta-gamma phase relaxation coefficient, J·s/m⁴ | Required; positive |

The grid is uniform, cell centered, periodic in x, and zero-flux at the two y
boundaries. Its physical width is `nx*dx` and height is `ny*dy`. These boundary
conditions are fixed in this version, not selectable input fields. The pairwise
relaxation coefficients enter `tau_effective * epsilon * dphi/dt = driving force`;
they are not automatically replaced by a calibration formula.

## Input parameters: time, solver controls, and output

| Key | Meaning and units | Requirement/default |
|---|---|---|
| `mode` | `coupled` advances phase and conservative carbon diffusion; `phase_only` advances phases at fixed chemical potential and resets carbon to thermodynamic closure | Default `coupled` |
| `dt` | Initial trial time step, s | Required; positive |
| `dt_max` | Cap on subsequent time-step proposals, s | Default `dt`; must be at least `dt` |
| `dt_min` | Smallest permitted trial time step, s | Default `dt/1048576`; positive and no greater than `dt` |
| `end_time` | Absolute stopping time in physical seconds, including time already present in a checkpoint | Required; positive |
| `output_dt` | Spacing of requested output times, s | Required; positive |
| `max_steps` | Absolute accepted-step limit, including steps already in a checkpoint | Default 100000000; at least 1 |
| `max_phase_change` | Allowed maximum local change of a phase fraction per accepted physical step | Default 0.04; greater than 0 and at most 0.2 |
| `cg_rtol` | Relative residual tolerance for the conjugate-gradient diffusion solve | Default 1e-11; positive |
| `cg_atol` | Absolute residual tolerance for the diffusion linear system in its code normalization | Default 1e-14; positive |
| `cg_max` | Maximum conjugate-gradient iterations per solve | Default 2000 |
| `closure_tol` | Maximum absolute error in the reconstructed carbon/chemical-potential relation; carbon-fraction units | Default 1e-10; positive |
| `energy_atol` | Absolute allowance in the reduced-energy acceptance test, J/m | Default `1e-13*sigma*max(nx*dx,ny*dy)`; nonnegative |

The integrator shortens steps to hit output times and the final time. Therefore
`dt = dt_max` and zero rejected steps do not by themselves guarantee every
accepted step has that duration. `output_dt` can affect the actual step sequence.
On rejection the trial time step is halved; an accepted step proposes up to 1.1
times its accepted duration, capped by `dt_max`.

The source also contains fixed internal controls, such as the energy relative
tolerance and retry limit. They are **not input keys**. Adding lines such as
`energy_rtol`, `max_retries`, or `check_energy` to an input file is an error.

The output directory, restart checkpoint, and MPI process-grid dimensions are
command-line arguments, not input parameters. The serial and parallel programs
accept the same physical input file. See `README.md` for build and invocation
details. Restarting loads stored fields and physical time; it does not rerun the
filling stage.

## Filling: native initial phase geometry

| Key | Meaning and units | Default and constraints |
|---|---|---|
| `initialization` | Native choice `lamella`, `sharp_smooth`, `flat`, `circle`, or opt-in command recipe `microsim` | Default `lamella` |
| `seed_height` | Height of the seed top or flat interface above y=0, m | Default `ny*dy/4`; strictly between 0 and `ny*dy` |
| `theta_columns` | Theta stripe width in x cells for lamellar filling | Default integer `nx/2`; at least 1 and less than `nx` |
| `radius` | Radius of the centered inclusion for `circle`, m | Default `min(nx*dx,ny*dy)/4` |
| `phase_a` | Phase below the flat interface or inside the circle | Default 0; valid IDs 0, 1, 2 |
| `phase_b` | Phase above the flat interface or outside the circle | Default 2; valid IDs 0, 1, 2; different from `phase_a` |

`lamella` creates an analytic diffuse alpha/theta seed below a gamma region.
Theta occupies the first `theta_columns` in x; alpha fills the remainder of the
product period. `sharp_smooth` starts from exact phase labels using the same
lamellar layout, then performs the interface-only preprocessing described below.
`flat` creates a two-phase analytic diffuse horizontal interface at `seed_height`.
`circle` creates a two-phase analytic diffuse circle centered at
`(nx*dx/2, ny*dy/2)`. The analytic diffuse profile uses `epsilon`.

In sharp filling, the center of zero-based cell `(i,j)` is
`((i+0.5)*dx,(j+0.5)*dy)`. A cell is in the product seed only when its y center is
strictly below `seed_height`. At the supplied G137 values, this gives rows 0–20
and theta columns 0–16; the remaining 120 seed columns are alpha. The discrete
theta share within the seed is `17/137`, not exactly 12%. Smoothed phase-volume
fractions can differ from these sharp labels.

The original parser reads the native geometry fields even when a selected
method does not use all of them. It still validates `seed_height`,
`theta_columns`, and the phase pair for every method, including `microsim`.
For `microsim`, the FILL commands determine the sharp geometry; the retained
native seed fields do not alter those labels. Retain valid values
for such fields; setting an unused seed dimension to zero is not a way to
disable it.

## Filling: opt-in MicroSim-style recipes

Choose `initialization = microsim` and add one or more FILL commands. The map
starts as pure gamma, operations apply in file order, and then the same existing
capillary-only smoothing runs before carbon initialization. Commands with any
other initialization are rejected. No FILL list is ignored silently.

```text
initialization = microsim
FILLCUBE = {1,0,0,0,16,20,0};
FILLCUBE = {0,17,0,0,136,20,0};
```

These G137 cubes reproduce the 17-theta/120-alpha sharp seed at rows 0–20.
`FILLCUBE`, `FILLCENTERBOX2D`, `FILLCYLINDER`, `FILLYJUNCTION2D`,
`FILLYJUNCTIONLAMELLAE2D`, `FILLADHEREDPEARLITEARM2D`, and
`FILLPEARLITEELLIPSE2D` are implemented. Full argument lists and mutually
alternative recipes are in `FILLING_METHODS.md` and the commented menu in
`examples/fe_c_microsim_cube.in`. No extra thermodynamic phase, gamma grain field,
3D domain or fluid velocity is introduced. Unsupported historical commands fail
with an explanation.

The new mode records a `filling-v1` suffix and the ordered operation names/values
in the physics signature. Keep that recipe for restart. Legacy inputs retain
their original signature; an equivalent-looking `microsim` mask deliberately
does not share a legacy initialization signature.

## Filling: initial carbon and chemical potential

Specify **exactly one** of the following keys. The initial chemical potential
is spatially uniform, and the actual carbon field is computed from the phase
mixture and the evaluated run-temperature thermodynamics.

| Key | Meaning | Units/default |
|---|---|---|
| `initial_mu` | Directly sets the uniform initial chemical potential | J/mol; no default choice |
| `gamma_carbon` | Uses pure-gamma carbon to calculate `initial_mu = 2*A_gamma*gamma_carbon + B_gamma` | Carbon atomic fraction; no default choice |
| `mean_carbon` | Solves for the uniform chemical potential giving this domain-average carbon on the initialized phase field | Carbon atomic fraction; no default choice |

For `sharp_smooth` and `microsim`, this composition initialization happens **after** smoothing.
The preprocessing has no carbon diffusion and consumes no physical time.
`gamma_carbon` is not a command to fill every phase with the same carbon value.

## Filling: capillary-only preprocessing

These controls affect `sharp_smooth` and `microsim`; they remain valid parsed
options for the other native methods but do not cause smoothing there.

| Key | Meaning | Default and constraints |
|---|---|---|
| `smooth_steps` | Maximum number of interface-only relaxation iterations | 100; at least 1 |
| `smooth_min_steps` | Minimum iteration before convergence-based stopping | Defaults to `smooth_steps`; 0 through `smooth_steps` |
| `smooth_energy_rtol` | Relative interface-energy change threshold | 1e-10; positive |
| `smooth_max_phase_change` | Maximum trial phase-fraction change used to select preprocessing step size | 0.02; greater than 0 and at most 0.2 |

Relaxation uses the existing gradient, obstacle, and three-phase interface
energies, without chemical driving or carbon transport. After the minimum
iteration, convergence requires both relative energy change below
`smooth_energy_rtol` and phase change below the fixed internal threshold 1e-5.
A state with exactly zero force can finish immediately. This stage is an
initialization operation, not a physical-time trajectory.

## Input parameters: optional moving window

| Key | Meaning and units | Default and constraints |
|---|---|---|
| `moving_window` | Enables y-direction window shifting | Integer default 0; any nonzero value enables it |
| `shift_trigger` | Front height within the current window that triggers a shift, m | Default `0.65*ny*dy` |
| `shift_target` | Target height after shifting, m | Default `0.35*ny*dy` |
| `shift_gamma_carbon` | Carbon atomic fraction of new gamma rows inserted at the top | Defaults to `gamma_carbon` in gamma-initialization mode, otherwise pure-gamma composition at the parsed initial chemical potential |

When enabled, the constraints are
`0 < shift_target < shift_trigger < ny*dy` and
`0 <= shift_gamma_carbon <= 1`. The code tracks cumulative displacement and
carbon added/removed. An analysis must restore that displacement to compare
physical front positions. In `mean_carbon` mode, the default
`shift_gamma_carbon` is evaluated while parsing, before the later domain-average
chemical-potential solve; explicitly choose this optional feed composition if
that mode is used with a moving window.

The supplied G137 inputs omit `moving_window`, so it remains disabled. Documenting
its options does not change that behavior.
