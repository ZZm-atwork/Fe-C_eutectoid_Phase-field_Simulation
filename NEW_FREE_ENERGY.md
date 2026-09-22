# Supplied reference free energies with phase-diagram slopes

The selected input method uses the supplied free-energy coefficients at
**999.68459774 K**, accompanied by composition slopes from the supplied phase
diagram. `Input.in` and `examples/fe_c_phase_diagram_10K.in` now use this method.
A 10 K undercooling still sets the isothermal simulation temperature to
**989.68459774 K**; the solver constructs the corresponding model parabolas from
the reference inputs. No free-energy export at that lower temperature is required
by this selected construction.

The law is named `phase_diagram_fixed_curvature`. The historical `linear_gp`
law and direct-at-temperature inputs remain available with their original
meanings. The earlier `fe_c_new_thermo_reference_only.in` is retained as a
separate historical diagnostic at Te.

## Reference data and contact compositions

The supplied convention is `G(c) = A*c*c + B*c + C`, with carbon mole fraction
`c` and coefficients in J/mol:

| Phase | Aeq | Beq | Ceq |
|---|---:|---:|---:|
| Alpha | 4724862.6373624019 | 22964.296428571855 | -42255.03701292308 |
| Theta surrogate | 100000000 | -49968654.66102279 | 6207741.2462954065 |
| Gamma | 172591.65731409981 | 19450.725112218221 | -42053.816371518253 |

An independent decimal least-squares reconstruction of the local data windows
reproduced the supplied alpha/gamma coefficients to relative 1.3e-12 or better.
The original CSV compositions are mole percent and were divided by 100. The
supplied text files already contain mole fractions.

The new input uses the actual alpha/gamma common tangent of these supplied
parabolas, with theta at c=0.25:

| Input | Selected fit contact | Phase-diagram invariant contact |
|---|---:|---:|
| `ceq_alpha` | 0.00088690859310575857 | 0.000886365 |
| `ceq_theta` | 0.25 | 0.25 |
| `ceq_gamma` | 0.034458832049305582 | 0.0344584275 |

The small differences between fit and diagram contacts are retained explicitly.
The gamma fit contact also supplies `gamma_carbon`, establishing one initial
chemical potential near 31345.33897721396 J/mol. This bulk composition is distinct
from the two gamma coexistence compositions shifted by undercooling.

The supplied theta curve lies about **8.234e-5 J/mol below** the reference
alpha/gamma tangent. Its minute derivative mismatch is also retained. No energy
offset or reference refit is applied to force exact equality. Theta curvature
A=1e8 is a chosen finite-susceptibility surrogate; one stoichiometric energy point
does not measure that curvature.

## Phase-diagram slopes

The selected slopes come from an anchored least-squares fit to all unique
exported points between Te and Te+10 K. The alpha/gamma export contains two
sampling passes; all distinct points in that window are included once. Exact
row numbers and residuals are recorded in the project evidence.

| Input key | dc_eq/dT (mole fraction/K) | Fit points including anchor | Maximum composition fit error |
|---|---:|---:|---:|
| `slope_alpha_on_alpha_gamma` | -4.6282639461434345e-6 | 8 | 5.48157e-8 |
| `slope_gamma_on_alpha_gamma` | -0.00031198492411304955 | 8 | 1.17793e-5 |
| `slope_theta_on_theta_gamma` | 0 | 4 | 0 |
| `slope_gamma_on_theta_gamma` | 0.00010753713136114516 | 4 | 5.30760e-7 |

These values are **dc/dT**, in carbon mole fraction/K. The older MicroSim
phase-diagram `slopes` convention was **dT/dc**: its reciprocal gives the new
quantity where finite. Theta's fixed stoichiometry is entered as dc/dT=0,
without requiring an infinite input. The two gamma-side slopes belong to
different phase pairs and must remain distinct.

These are not grand-potential coefficient derivatives. Do not copy them into
`dq2_dT_*`, `dq1_dT_*`, or `dq0_dT_*`. The parser requires the matching input
family and rejects mixed-law parameters.

## Reference-anchored MicroSim construction

The older MicroSim pair-boundary construction obtains product B and C from
shifted product/gamma coexistence compositions. This implementation uses the
change in that construction from Te, added to the supplied ABC, so the reference
fit and its residual are preserved.

For product phase p (alpha or theta), let `sp` and `sg` be the product- and
gamma-side slopes on that same branch:

```text
delta_T = T - Te = -undercooling
dp = sp * delta_T
dg = sg * delta_T
cp(T) = ceq_p + dp
cg_on_p_gamma(T) = ceq_gamma + dg

A_p(T) = Aeq_p
B_p(T) = Beq_p + 2*(Aeq_gamma*dg - Aeq_p*dp)
C_p(T) = Ceq_p + Aeq_p*dp*(2*ceq_p + dp)
                   - Aeq_gamma*dg*(2*ceq_gamma + dg)
```

All reference curvatures remain fixed. Gamma ABC remains equal to its reference
values: holding gamma's B/C increments at zero chooses a common affine energy
reference. Fixed curvature is an additional model assumption, separate from
that reference choice. Product C retains the **full quadratic temperature
increment**. Linearizing C or inserting its reference derivative into
`linear_gp` would be a different model.

The C++ helper is `at_phase_diagram_undercooling` in both combined files.
At zero undercooling it returns the supplied ABC unchanged. At nonzero
undercooling it checks the shifted contact compositions, preserves each D, and
does not alter Vm, tau, sigma, epsilon, the triple penalty, or numerical controls.
The new law and its contact/slope inputs are included in restart compatibility.

## Independent 10 K thermodynamic calculation

The reference-plus-slope construction gives:

| Phase | A(T), J/mol | B(T), J/mol | C(T), J/mol |
|---|---:|---:|---:|
| Alpha | 4724862.6373624019 | 23603.858102865859 | -42293.428311015959 |
| Theta | 100000000 | -49969025.861257277 | 6207753.8378329016 |
| Gamma | 172591.65731409981 | 19450.725112218221 | -42053.816371518253 |

At the initial common chemical potential, independent analysis gives
product-minus-gamma grand-potential differences of **-37.8457083 J/mol** for
alpha and **-80.2089479 J/mol** for theta. Both products therefore have a chemical
driving advantage in this initialized model. These are calculated thermodynamic
values, not saved growth-rate or solver-residual measurements.

The shifted pair contacts are approximately:

| Branch | Product carbon | Gamma carbon |
|---|---:|---:|
| Alpha/gamma | 0.00093319123256719 | 0.037578681290436 |
| Theta/gamma | 0.25 | 0.033383460735694 |

The alpha/gamma tangency holds to the supplied fit's numerical precision.
Theta/gamma retains the small reference mismatch; its exact numerical common
tangent is consequently very slightly displaced from the prescribed contacts.

## Scope and limitations

- The supplied gamma-containing phase boundaries are stable above Te. Their
  linear use 10 K below Te is a metastable extrapolation.
- Four boundary slopes with fixed curvature constrain the alpha/gamma and
  theta/gamma branches. They do not independently fit the third alpha/theta
  branch. At the run temperature, the construction predicts alpha/theta contacts
  near 0.00080123788793405 and 0.2500010059864; the supplied stable diagram gives
  interpolated alpha 0.00080176029934515 and stoichiometric theta 0.25.
- Some shifted coexistence contacts lie outside the narrow reference free-energy
  fit windows. The parabolas are a model extrapolation there.
- Material normalization, long-time cooperative growth and velocity matching are
  not established by this thermodynamic reconstruction or by short tests.

The actual frozen MicroSim campaign used a different direct grand-potential
interpolation between two temperature endpoints. The recovered older MicroSim
source and generic pair-boundary implementation support the requested
reference-plus-phase-diagram method. The new named law preserves both histories
rather than silently changing the old `linear_gp` interpretation.

## Calibration and long-run inputs

The former Te-only 0.015 µs calibration completed with 5,475 accepted steps and
755 rejected attempts. That is historical reference-temperature evidence and
must not be treated as the new undercooled run's timestep statistics. Current
local/parallel validation and actual runtime observations belong in
[VALIDATION.md](VALIDATION.md) and the per-run reports.

`Input.in` and `examples/fe_c_phase_diagram_10K.in` use a **1.5 µs example
horizon**. The smoke example uses the tested short input's parsed values and a
**0.015 µs horizon**. The requested Sol resource plan is one node, **16 MPI ranks
/ 4 × 4**, one CPU per rank, HCOLL disabled, and a **72-hour wall-time limit**
(maximum 1,152 allocated core-hours). Production uses a separately frozen input
and a horizon selected from its benchmark and domain/step constraints. Slurm
wall time is not simulated physical time. This document does not assert a
particular submission or completion status.

The executable still does not save a per-step accepted-dt history or a numerical
true-CG-residual history. Aggregate counts and output-time `dt_next` cannot be
relabeled as those missing measurements.

Project evidence: `research/pearl_v04_phase_diagram_20260921/audit/REPORT.md`,
`AUDIT.json`, `SLOPE_FIT_POINTS.csv`, and `MICROSIM_LAW_REVIEW.md` record the
independent arithmetic, data provenance, extraction window and source review.
The earlier raw files, reference-only calibration and frozen historical artifacts
remain preserved.
