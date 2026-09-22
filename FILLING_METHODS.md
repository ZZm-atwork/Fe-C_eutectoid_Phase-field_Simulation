# Filling the 2D Pearl domain

The filling section of an input file describes **the initial sharp phase map**. Thermodynamic parameters, diffusion, interface kinetics, smoothing and time integration belong to the other input sections. Filling does not consume physical time. The existing `sharp_smooth` initialization smooths the labels before initializing carbon and chemical potential.

Pearl has three fields: **alpha/ferrite = 0, theta/cementite = 1, gamma/austenite = 2**. All methods below act on the same two-dimensional grid. They do not create crystallographic grain orientations or extra gamma fields.

## Units and ordering

- Positions are global integer **grid indices**, not meters or cell-center coordinates. A coordinate distance of one represents one grid spacing `dx`. Fractional junction coordinates are allowed only for adhered-arm and pearlite-ellipse methods.
- Widths, lengths and radii are in **grid cells**. Fractions/aspect ratios are dimensionless. Angles are degrees counterclockwise from +X.
- Cube bounds are **inclusive**: x=0..16 occupies 17 columns. Box bounds in center-box and Y methods are **half-open**: `[0,137)` includes columns 0..136.
- The map starts as gamma everywhere. Commands run in the order written. Overlay methods preserve cells outside their shapes; later overlapping shapes overwrite earlier labels. The center-box and Y methods replace the entire map.
- Geometry is clipped at domain edges; it does not automatically wrap across the periodic X seam. For a shape that crosses a seam, explicitly describe both pieces.
- A command is `METHOD = {comma-separated values};`; the final semicolon is optional. Do not put words such as `alpha` inside numeric lists. Malformed lists, invalid phase IDs, noninteger index arguments, out-of-range boxes and unsupported methods produce errors.
- Whole-map Y and center-box commands are best placed first, followed by overlay commands.

## Original G137 single-period filling

The completed MicroSim G137 campaign used this sharp map: theta has 17 columns and alpha 120 columns, both covering y=0..20, with gamma above. At `dx=2.79e-9 m`, the domain is `137 × 250` grid points. This example uses the correct parser order, **phase,xlo,ylo,zlo,xhi,yhi,zhi**:

```text
FILLCUBE = {1, 0, 0, 0, 16, 20, 0};
FILLCUBE = {0, 17, 0, 0, 136, 20, 0};
```

The older fill first included x=17 in theta, then overwrote it with alpha. The explicit disjoint ranges above produce the same labels. The current input's existing `seed_height`/`theta_columns` recipe remains available; use either that recipe or an explicit filling recipe as described in the input documentation.

## Implemented methods

### FILLCUBE — rectangular phase region

```text
FILLCUBE = {phase, xlo, ylo, zlo, xhi, yhi, zhi};
```

All bounds inclusive; `zlo=zhi=0` is required. Paint alpha or theta inside the rectangle. As in ordinary MicroSim `FILLCUBE`, requesting the last phase (`2`, gamma) only recomputes the residual background; it **does not erase existing alpha/theta**. Use a whole-map initializer where a separate inside/outside background is needed.

### FILLCENTERBOX2D — one box and a specified background

```text
FILLCENTERBOX2D = {xlo, xhi, ylo, yhi, inside_phase, outside_phase};
```

Replace the entire map: inside the half-open box receives `inside_phase`, everywhere else receives `outside_phase`. Phase IDs must differ. Bounds must fit within the grid.

### FILLCYLINDER — circular 2D cross-section

```text
FILLCYLINDER = {phase, center_x, center_y, zlo, zhi, radius};
```

Require `zlo=zhi=0`; the mask is `(x-center_x)^2+(y-center_y)^2 <= radius^2`. Radius is positive, center is an integer grid location, and the circle may be clipped at an edge. This is a disk in the 2D model, not a three-dimensional cylinder simulation. Ordinary last-phase behavior is the same as `FILLCUBE`.

### FILLYJUNCTION2D — three ferrite arms in austenite

```text
FILLYJUNCTION2D = {xlo, xhi, ylo, yhi,
                  junction_x, junction_y, angle0, angle1, angle2,
                  arm_length, arm_half_width, ferrite_phase, austenite_phase};
```

Write each command on **one input line**; the layout here is for readability. The box is half-open. Three straight arms start at the junction; each spans `0 <= distance_along_arm <= arm_length` and `abs(distance_across_arm) <= arm_half_width`. Ferrite occupies the arm union; austenite occupies all other cells, including outside the box. This replaces the entire map.

### FILLYJUNCTIONLAMELLAE2D — alternating phases along three arms

```text
FILLYJUNCTIONLAMELLAE2D = {xlo, xhi, ylo, yhi,
                         junction_x, junction_y, angle0, angle1, angle2,
                         arm_length, arm_half_width, ferrite_width, cementite_width,
                         ferrite_phase, cementite_phase, austenite_phase};
```

Use one physical input line. Box and arm geometry match the simple Y. Lamellae alternate **along each arm**, with period `ferrite_width+cementite_width`: ferrite first, then theta. First listed arm wins where arms overlap. Outside the arms/box is austenite; this replaces the entire map. A three-phase recipe normally ends with `0,1,2`.

This method exists in the supplied older MicroSim Y-lamella files. Their file presence is not evidence that this geometry was used in the completed G137/G274/G411 campaign.

### FILLADHEREDPEARLITEARM2D — an attached alpha/theta chain

```text
FILLADHEREDPEARLITEARM2D = {junction_x, junction_y, angle,
                          arm_length, pair_count, cementite_fraction,
                          peak_half_width, contact_half_width, bridge_radius,
                          ferrite_phase, cementite_phase};
```

Use one physical input line. Overlay one chain on the current background. There are `pair_count` integer pairs; period=`arm_length/pair_count`. Each pair places alpha over `(1-cementite_fraction)*period`, then theta over the remaining span. The transverse lobe profile is sinusoidal, with specified peak and contact half-widths. Outer endpoints taper to zero. A positive `bridge_radius` adds a ferrite disk at the start and a theta disk at the end; zero disables bridges. The end disk takes precedence if both overlap.

Require positive length/count/peak, `0 < cementite_fraction < 1`, `0 <= contact_half_width <= peak_half_width`, and nonnegative bridge radius. Historical Case I used this overlay on a three-gamma-grain map. Pearl preserves the **alpha/theta overlay geometry on one gamma background**, not the missing grain-boundary physics.

### FILLPEARLITEELLIPSE2D — an elliptical pearlite envelope

```text
FILLPEARLITEELLIPSE2D = {junction_x, junction_y, angle, s_start,
                       period_count, period, aspect_ratio, cementite_fraction,
                       bridge_radius, ferrite_phase, cementite_phase};
```

Use one physical input line. Overlay one ellipse whose long axis follows the given arm direction, starting `s_start` cells from the junction. Total length=`period_count*period`; semimajor axis=`length/2`, semiminor axis=`aspect_ratio*semimajor`. **Aspect ratio means b/a**, not mathematical eccentricity. Theta occupies the centered portion of each longitudinal period, with fraction `cementite_fraction`; alpha occupies the rest. Optional ferrite disks at both ends make connected arrays. These disks override theta within their masks.

The number of periods is a positive integer; lengths/aspect ratio are positive; `s_start` and bridge radius are nonnegative. Historical Case II used these overlays on three distinct gamma grains. Their alpha/theta masks are available here, but the three-grain background is not.

## Historical methods that remain unsupported

These names are documented for migration. They are **not accepted as working Pearl filling commands**. The solver gives an error rather than silently ignoring or remapping them.

| MicroSim method | Why it is not enabled in this 2D/three-phase Pearl model |
|---|---|
| `FILLTHREEGRAINYJUNCTION2D` | Historical inputs assign three different gamma grain fields (IDs 2,3,4). Pearl has one gamma field. Mapping them all to 2 removes grain boundaries; mapping them to alpha/theta/gamma changes phase physics. |
| `FILLYJUNCTIONPEARLITE2D` | Legacy irregular seeded packet method requires at least five fields and overlays the three-grain background. A single-gamma adaptation would need its own explicit definition and validation. |
| `FILLELLIPSE` | Generic legacy routine differs from the dedicated pearlite ellipse. Its parser truncates fractional geometry and its rotation reads an uninitialized angle; it is not copied as a validated method. |
| `FILLCYLINDERNEXLP` | Variant permitting explicit last-phase painting; outside the selected deterministic method set. Ordinary `FILLCYLINDER` is supported with its original residual-background meaning. |
| `FILLSPHERE`, `FILLSPHERENEXLP` | Three-dimensional ball geometry. Pearl remains 2D; use a supported disk when a 2D cross-section is intended. |
| `FILLCYLINDERRANDOM`, `FILLCYLINDERRANDOMNEXLP` | Random particle-placement routines with nominal rather than exact area fraction, legacy exclusion semantics and no explicit input seed. No new random-initialization contract is introduced here. |
| `FILLSPHERERANDOM`, `FILLSPHERERANDOMNEXLP` | Random three-dimensional spheres; both the 3D model and a new reproducible random-placement contract are outside this change. |
| `FILLCUBERANDOM`, `FILLCUBEPATTERN` | Legacy random variant assignments, time-seeded RNG and matrix assumptions are not equivalent to the fixed three physical phase fields. |
| `FILLVORONOI2D`, `FILLVORONOI3D` | Arbitrary grain variants are not provided by Pearl; the legacy implementation also has uninitialized exclusion-state and label-selection defects. 3D is additionally unavailable. |
| `FILLCUBEVELOCITY` | Initializes lattice-Boltzmann fluid velocity/density, neither of which exists in this diffusion/phase-field solver. |

## Provenance and validation scope

The completed q0 and G137/G274/G411 geometry campaigns used `FILLCUBE`. Supplied older pearlite files contain Y lamellae, three-grain Y maps, attached arms and elliptical envelopes. Other listed methods appear in the MicroSim library or unrelated examples; they are not established campaign-tested pearlite inputs.

Source details, argument orders, differing bounds and known legacy defects are recorded in [the filling inventory](../research/pearl_v04_readability_20260920/FILLING_INVENTORY.md). Supporting geometry checks compare compatible masks against the original MicroSim functions. Passing a mask comparison does not establish stability or cooperative growth for a new scientific geometry. No long run is authorized or performed merely by adding these input methods.
