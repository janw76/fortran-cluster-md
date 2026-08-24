# CODE_MAP.md — Detailed File-by-File Reference

## Include Files (Common Blocks)

### `const.inc` — Constants and Array Limits
**Common block:** (none — all `PARAMETER` declarations)

| Symbol | Value | Meaning |
|--------|-------|---------|
| `mxa` | 6000 | Max atoms |
| `mxnlist` | `mxa*mxa/4` (9,000,000) | Max neighbor list entries |
| `mxsorts` | 10 | Max atomic species |
| `cokb` | 1.38066×10⁻²³ J/K | Boltzmann constant (SI) |
| `cikb` | 1.38066×10⁻⁶ | Boltzmann constant (internal units: 10⁻¹⁷ J/K) |
| `comu` | 1.66056×10⁻²⁷ kg | Atomic mass unit (SI) |
| `cimu` | 1.66056 | Atomic mass unit (internal: 10⁻²⁷ kg) |
| `cona` | 6.022045×10²³ | Avogadro's number |
| `copl` | 6.62618×10⁻³⁴ Js | Planck's constant |
| `pi` | 3.14159265 | π |

### `atomc.inc` — Atomic Coordinates and Derivatives
**Common block:** `/ATOMC/`

| Variable | Type | Description |
|----------|------|-------------|
| `x0(mxa), y0(mxa), z0(mxa)` | double | Positions (Å) |
| `x1..x5, y1..y5, z1..z5` | double | 1st through 5th time derivatives |
| `x ≡ x0, vx ≡ x1, ax ≡ x2` | equivalence | Aliases for readability |

The higher derivatives (3rd–5th) are used by Gear integrators. The Verlet integrators only use positions (x0), velocities (x1), and accelerations (x2).

### `atomp.inc` — Atomic Parameters
**Common block:** `/ATOMP/`

| Variable | Type | Description |
|----------|------|-------------|
| `natom` | int | Current total atom count |
| `atsorts` | int | Number of species |
| `api(mxa)` | int | Species index per atom |
| `apnum(mxsorts)` | int | Atom count per species |
| `apmas(mxsorts)` | double | Mass per species (internal units) |
| `apsig(mxsorts)` | double | LJ σ per species (Å) |
| `apeps(mxsorts)` | double | LJ ε per species (internal energy) |
| `appr(mxsorts)` | double | Repulsive LJ exponent |
| `appa(mxsorts)` | double | Attractive LJ exponent |
| `apsym(mxsorts)` | char*2 | Atomic symbol per species |
| `apsymi(mxa)` | char*2 | Atomic symbol per atom |
| `apth(mxsorts)` | char*1 | Thermostat flag per species |
| `att` | double | System temperature (K) |
| `atskin` | double | Stillinger cluster criterion radius (Å) |
| `atrcut` | double | Potential cutoff radius (Å) |
| `atrskin` | double | Neighbor list skin radius (Å) |
| `atnlist(mxnlist)` | int | Neighbor list (Verlet list) |
| `atnidx(mxa)` | int | Per-atom index into neighbor list |
| `atnenum` | int | Neighbor count for Frenkel criterion |
| `Q` | double | Nosé-Hoover coupling mass |
| `PSI` | double | Nosé-Hoover friction coefficient |
| `s` | double | Nosé-Hoover auxiliary coordinate |
| `rndctrl` | char*1 | Random velocity control |

**Equivalences:** `ats ≡ apsym(1)`, `atm ≡ apmas(1)`, `atr ≡ apsig(1)`, `ate ≡ apeps(1)`, etc. — provide backward-compatible single-species access.

### `boxpp.inc` — Box Parameters
**Common block:** `/BOXPP/`

| Variable | Type | Description |
|----------|------|-------------|
| `boxx, boxy, boxz` | double | Box dimensions (Å) |
| `boxxn, boxyn, boxzn` | double | Target box dimensions for resize |
| `boxxs, boxys, boxzs` | double | Max resize step per axis (Å) |
| `beng` | double | Wall energy parameter |
| `bpot` | double | Wall potential exponent |
| `bcor` | char*1 | Show box corners in XYZ output |
| `bshr` | char*1 | Box resize active (X/-) |
| `bperx, bpery, bperz` | char*1 | PBC per axis (X/-) |
| `boxt` | int | Steps between box resizes |

### `cvtpp.inc` — Cluster Statistics
**Common block:** `/CVTPP/`

| Variable | Type | Description |
|----------|------|-------------|
| `cl(mxa)` | int | Stillinger cluster label per atom |
| `clf(mxa)` | int | ten Wolde/Frenkel cluster label per atom |
| `cs(mxa)` | int | Stillinger cluster sizes |
| `csf(0:mxa)` | int | Frenkel cluster sizes |
| `cc(mxa), ccf(mxa)` | int | "Visited" markers for cluster search |
| `csi(0:mxa)` | int | Cluster size distribution (Stillinger) |
| `csif(0:mxa)` | int | Cluster size distribution (Frenkel) |
| `cek(mxa), cekf(mxa)` | double | Kinetic energy per cluster |
| `ceks(0:mxa), ceksf(0:mxa)` | double | KE summed by cluster size |
| `cgs(mxa), cgsf(0:mxa)` | int | Carrier gas atoms per cluster |
| `cgsa(mxa), cgsfa(0:mxa)` | int | Carrier gas summed by size |
| `cnum, cnumf` | int | Number of clusters found |
| `cbiggest, cbiggestf` | int | Largest cluster size |
| `cfpt(mxa)` | int | First passage time per size (Stillinger) |
| `cfptf(mxa)` | int | First passage time per size (Frenkel) |
| `csover` | int | Min cluster size for nucleation abort |
| `clwait` | int | Persistence time before abort |
| `clbrkt` | int | Timestep when large cluster first found |
| `clnuctrl` | char*1 | Nucleation checking enabled (X/-) |
| `ccn, cen` | int | Condensation/evaporation counters |
| `ctrat` | int | Rate statistics interval |
| `mcts, omcts` | double | Cluster temperature averages |

### `energ.inc` — Energy
**Common block:** `/ENERG/`

| Variable | Type | Description |
|----------|------|-------------|
| `ekin` | double | Total kinetic energy |
| `epot` | double | Total potential energy |
| `H` | double | Hamiltonian (E_kin + E_pot + NH terms) |
| `sumv2` | double | Sum of v² (for Nosé-Hoover) |
| `ekinA(mxsorts)` | double | Kinetic energy per species |
| `ewall, ewallsum` | double | Wall energy (instantaneous, summed) |
| `fwall, fwallsum` | double | Wall force (instantaneous, summed) |
| `ethst` | char*1 | Thermostat type (X/N/-) |
| `eclth` | char*1 | Cluster thermostat flag |
| `ecltim` | int | Cluster thermostat timer |

### `filep.inc` — File Parameters
**Common block:** `/FILEP/`

All file names (short and full-path versions) and output intervals for every output file type. See the README for the full list.

### `timep.inc` — Time Parameters
**Common block:** `/TIMEP/`

| Variable | Type | Description |
|----------|------|-------------|
| `dt` | double | Timestep (fs) |
| `ntim` | int | Total timesteps |
| `btim` | int | Starting timestep |
| `sctim` | int | Screen output interval |
| `nbtim` | int | Neighbor update interval |
| `chgNtim` | int | Atom count change interval |

### `dmmpp.inc` — Displacement Statistics
**Common block:** `/DMMPP/`

| Variable | Type | Description |
|----------|------|-------------|
| `dmmh(0:255)` | double | Ring buffer of max displacements |
| `dmms` | double | Cumulative displacement since last neighbor update |
| `dmmax` | double | Historical maximum displacement |
| `dmmav` | double | Moving average displacement |
| `dmmsu` | double | Running sum for moving average |
| `dmmst` | int | Step counter since last neighbor update |
| `dmmi` | int | Ring buffer index |

### `connp.inc` — Connectivity (Alpha/Unused)
**Common block:** (partially declared, never used in simulation)

Declares bond, angle, and torsion connectivity arrays. Was never completed — the code uses non-bonded LJ interactions only.

---

## Source Files — Detailed Reference

### Main Program

#### `main.f` — `program cluster`
**Purpose:** Top-level simulation driver.
**Author:** Stephan Wonczak (1997), extended by Jan Wedekind (2005–2007)
**Includes:** `const.inc`, `atomc.inc`, `filep.inc`, `timep.inc`, `energ.inc`, `boxpp.inc`, `atomp.inc`, `cvtpp.inc`

**Flow:**
1. `kpsInit(0)` → initialize RNG
2. `readdata()` → parse `.3d` input
3. Open output files conditionally
4. Print simulation summary
5. Initialize: restart file OR `initconf()`
6. Set Nosé-Hoover initial values
7. `stoddard(-1)` → initial cluster search
8. **Main loop** (label 10): output → thermostat → resize → nucleation check → `integrator(i)` → increment
9. Loop exits on completion or nucleation abort (goto 20)

**Calls:** `kpsInit`, `readdata`, `readRST`, `initconf`, `stopdead`, `stoddard`, `wrtXYZ`, `wrtVEL`, `wrtACC`, `wrtRST`, `wrtCSI`, `wrtAVE`, `thermostat`, `resize`, `chgNatom`, `chk_nucleation`, `integrator`, `rates`

---

### Initialization Files

#### `init-cube1.f` — `subroutine initconf()`
**Purpose:** Place atoms on a cubic lattice, assign random velocities, warm up.
- Computes cube side length from N^(1/3)
- Places atoms at regular intervals within the box
- For multi-species: randomly swaps positions to mix species
- Calls `stopdead(-1)` → `neighbor(0)` → `AccAtom()` → `AccWall()`
- Performs 100 warmup integration steps
- Contains helper `subroutine swap(a,b)` for position exchange

**Includes:** `const.inc`, `atomc.inc`, `atomp.inc`, `boxpp.inc`, `energ.inc`

#### `init-center.f` — `subroutine initconf()`
**Purpose:** Similar to init-cube1 but centers the initial cluster in the box. Used for diffusion studies.

#### `init-ball1.f` — `subroutine initconf()`
**Purpose:** Spherical starting configuration.

#### `init-2atom.f` — `subroutine initconf()`
**Purpose:** Two-atom test system.

#### `init-2cubes.f` — `subroutine initconf()`
**Purpose:** Two separate cubes of atoms (collision studies).

---

### Integrators

#### `int-verlet-mav.f` — `subroutine integrator(step)`
**Purpose:** Velocity-Verlet with adaptive neighbor list (main production integrator).
**Algorithm:** Swope-Andersen-Berens-Wilson variant:
1. `x(t+dt) = x(t) + [v(t) + ½a(t)·dt]·dt`
2. `v(t+½dt) = v(t) + ½a(t)·dt`
3. Track max displacement → trigger neighbor list rebuild if skin exceeded
4. Reset accelerations, call `AccAtom()` + `AccWall()`
5. `v(t+dt) = v(t+½dt) + ½a(t+dt)·dt`
6. Compute `ekin`

**Includes:** `const.inc`, `atomp.inc`, `atomc.inc`, `timep.inc`, `energ.inc`, `dmmpp.inc`
**Calls:** `neighbor(step)`, `AccAtom()`, `AccWall()`

#### `int-verlet-nh-mav.f` — `subroutine integrator(step)`
**Purpose:** Velocity-Verlet with Nosé-Hoover thermostat coupling.
- Modified equations of motion: includes `-PSI·v` friction terms
- Updates Nosé-Hoover variables `s` and `PSI`
- **Iterative velocity/PSI convergence** (up to 1000 iterations, tolerance 10⁻¹⁰)
- Only works with single atom species

**Includes:** Same as above
**Calls:** `neighbor(step)`, `AccAtom()`, `AccWall()`

#### `int-verlet.f` — `subroutine integrator(step)`
**Purpose:** Basic Verlet without displacement tracking. For simple systems.

#### `int-verlet-dmm.f` — `subroutine integrator(step)`
**Purpose:** Verlet with displacement tracking but without moving average optimization.

#### `int-euler1.f` — `subroutine integrator(step)`
**Purpose:** Simple Euler method (first-order). Teaching/testing only.

#### `int-euler2.f` — `subroutine integrator(step)`
**Purpose:** Euler with acceleration correction.

#### `int-gear3.f` through `int-gear5.f` — `subroutine integrator(step)`
**Purpose:** Gear predictor-corrector algorithms, 3rd through 5th order. Use higher derivatives from `atomc.inc`.

---

### Force Calculations — Atom-Atom

#### `for-a-lj-n-pbc-m.f` — `subroutine AccAtom()`
**Purpose:** **Primary production force routine.** LJ 12-6 with cutoff, PBC (minimum image convention), neighbor list, multi-species.
- Handles species 1-1 (full LJ), 2-2 (LJ), and 1-2 (cross, geometric mixing) interactions
- Uses precomputed σ⁶ and σ¹² constants
- Minimum image convention via `dx = dx + boxx*(bf1 - int(bf2 + dx*ibx))`
- Accumulates forces in local arrays, then converts to accelerations using per-atom mass
- O(N·neighbors) complexity

**Includes:** `const.inc`, `atomc.inc`, `atomp.inc`, `energ.inc`, `boxpp.inc`

#### `for-a-lj-n-pbc.f` — `subroutine AccAtom()`
**Purpose:** Single-species LJ with PBC and neighbor list.

#### `for-a-lj-n-pbcsw.f` — `subroutine AccAtom()`
**Purpose:** LJ with shifted/truncated potential and PBC.

#### `for-a-lj-n.f` — `subroutine AccAtom()`
**Purpose:** LJ with neighbor list, no PBC.

#### `for-a-lj1.f` through `for-a-lj4.f` — `subroutine AccAtom()`
**Purpose:** Progressively optimized LJ without neighbor list (O(N²)). Historical development sequence.

#### `for-a-dummy.f` — `subroutine AccAtom()`
**Purpose:** No-op. Used when atom-atom forces are not needed.

---

### Force Calculations — Atom-Wall

#### `for-w-i9.f` — `subroutine AccWall()`
**Purpose:** Repulsive r⁻⁹ wall potential on all box faces. For confined (non-PBC) simulations.

#### `for-w-i12.f` — `subroutine AccWall()`
**Purpose:** Repulsive r⁻¹² wall potential.

#### `for-w-iX.f` — `subroutine AccWall()`
**Purpose:** Repulsive r⁻ˣ wall, exponent read from `.3d` file.

#### `for-w-i9-dbl.f` — `subroutine AccWall()`
**Purpose:** Double-wall variant with exponential decay.

#### `for-w-i9g.f` — `subroutine AccWall()`
**Purpose:** r⁻⁹ wall with ghost particles.

#### `for-w-lj.f` — `subroutine AccWall()`
**Purpose:** Full LJ 12-6 wall potential.

#### `for-w-dummy.f` — `subroutine AccWall()`
**Purpose:** No-op. **Used for all PBC simulations** (no walls needed).

---

### Neighbor Lists

#### `neighbor-pbc-mav.f` — `subroutine neighbor(step)`
**Purpose:** **Primary production neighbor list.** PBC, adaptive skin radius via moving average displacement tracking.
- Wraps all positions into primary box via `dmod`
- Builds Verlet neighbor list: for each atom i, lists atoms j>i within `atrskin`
- Encoded as: `[-i, j1, j2, ..., -i+1, j1, ...]` in `atnlist`
- Dynamically adjusts `atrskin` based on observed displacement rate
- Fatal error if skin exceeded before rebuild
- Prints diagnostic: step, list size, displacement stats

**Includes:** `const.inc`, `atomc.inc`, `atomp.inc`, `boxpp.inc`, `timep.inc`, `dmmpp.inc`

Silent variant: the same file compiled with `-fpp -DSILENT` (object still named `neighbor-pbc-mav-silent.o`) suppresses the per-step diagnostic write. For production runs.

#### `neighbor-pbc.f`
**Purpose:** PBC neighbor list without adaptive skin.

#### `neighbor-pbcsw-mav.f`
**Purpose:** PBC with shifted potential cutoff and moving average.

#### `neighbor.f`
**Purpose:** Basic neighbor list, no PBC.

#### `neighbor-dummy.f`
**Purpose:** No-op. For O(N²) force routines that don't use neighbor lists.

---

### Cluster Identification (Stoddard Algorithm)

#### `stoddard-pbc-n-m-fre.f` — `subroutine stoddard(step)`
**Purpose:** **Most complete cluster finder.** Implements both Stillinger and ten Wolde/Frenkel cluster definitions with PBC and multi-species support.

**Stillinger criterion:** Two atoms belong to the same cluster if their distance < `atskin` (typically 1.8σ).

**ten Wolde/Frenkel criterion:** An atom is a "liquid monomer" if it has ≥ `atnenum` (typically 5) neighbors of species 1 within `atskin`. Frenkel clusters are then built only from liquid monomers.

**Algorithm:**
1. Build full neighbor list with back-references using PBC minimum image
2. Identify liquid monomers (Frenkel): atoms with ≥ atnenum species-1 neighbors
3. Stillinger cluster search: flood-fill from each unvisited atom through neighbor connections
4. Frenkel cluster search: flood-fill among liquid monomers only
5. Compute statistics: size distribution, kinetic energy per cluster, first passage times

**Includes:** `const.inc`, `atomp.inc`, `atomc.inc`, `cvtpp.inc`, `boxpp.inc`

#### `stoddard-pbc-n-m-fre-mov.f`
**Purpose:** Same as above plus tracking of cluster movement for visualization.

#### `stoddard-pbc-n-m.f`
**Purpose:** Multi-species with PBC, Stillinger only (no Frenkel).

#### `stoddard-pbc-n.f`
**Purpose:** Single-species with PBC and neighbor list optimization.

#### `stoddard-pbcsw-n.f`
**Purpose:** Shifted potential variant.

#### `stoddard-pbc.f`
**Purpose:** PBC without neighbor list (O(N²) distance checks).

#### `stoddard.f`
**Purpose:** Basic cluster search, no PBC.

#### `stoddard-dummy.f`
**Purpose:** No-op.

---

### Thermostat

#### `scalev-m.f` — `subroutine thermostat(step)`
**Purpose:** **Primary thermostat.** Isokinetic velocity scaling for multi-species.
- Per species: compute kinetic energy, compute scaling factor `kf = sqrt(E_target / E_actual)`
- Only scale species marked with `apth(i)='X'`
- Always scale at step -1 (initialization)
- In carrier-gas simulations: only the carrier gas (He) is thermostated

#### `scalev.f` — `subroutine thermostat(step)`
**Purpose:** Single-species velocity scaling.

#### `scalev-nh-m.f` — `subroutine thermostat(step)`
**Purpose:** Combined Nosé-Hoover + velocity scaling. Scales Nosé-Hoover variables when called at initialization.

#### `scalev-clu.f` — `subroutine thermostat(step)`
**Purpose:** Cluster thermostat — scales velocities of cluster atoms separately from vapor.

---

### I/O — Input

#### `read-3d.f` — `subroutine readdata()`
**Purpose:** Parse the `.3d` configuration file. Handles version compatibility (v1.02 through v1.30).
- Reads from command-line argument (first arg)
- Sequential fixed-format reads with version branching
- Converts input units (nm→Å, u→internal mass, K→internal energy)
- Initializes displacement statistics, cluster variables, and zeroes all coordinate arrays
- Contains helper `integer function fspc(a)` to find first space in a string

#### `read-rst.f` — `subroutine readRST(step, fnum)`
**Purpose:** Read binary restart file. Restores positions, velocities, accelerations (all 6 derivative levels), and energies.

#### `read-xyz.f` — `subroutine readXYZ(step, fnum)`
**Purpose:** Read one frame from an XYZ trajectory file. Used by helper programs.

#### `read-xyz-m.f` — `subroutine readXYZ(step, fnum)`
**Purpose:** Read XYZ for multi-species systems.

#### `read-vel.f` — `subroutine readVEL(step, fnum)`
**Purpose:** Read velocity data file.

---

### I/O — Output

#### `wrt-xyz-m-vmd.f` — `subroutine wrtXYZ(step, fnum)`
**Purpose:** **Primary coordinate writer.** VMD-compatible XYZ format for multi-species. Also writes separate Stillinger and Frenkel cluster visualization files.

#### `wrt-xyz.f` — `subroutine wrtXYZ(step, fnum)`
**Purpose:** Basic XYZ output, single species.

#### `wrt-xyz1.f` — Color-coded by cluster size (monomers=H, dimers=N, trimers=O, etc.)

#### `wrt-xyz2.f` — Outputs only the largest cluster.

#### `wrt-xyz3.f` — Outputs clusters larger than a threshold.

#### `wrt-xyz-m.f` / `wrt-xyz-m1.f` — Multi-species XYZ variants.

#### `wrt-xyz-m-ats1.f` — Outputs only species 1 atoms.

#### `wrt-csi.f` — `subroutine wrtCSI(step, fnum)`
**Purpose:** Write Stillinger cluster size distribution.

#### `wrt-csi-fre.f` — `subroutine wrtCSI(step, fnum)`
**Purpose:** Write both Stillinger and ten Wolde/Frenkel cluster size distributions to separate files.

#### `wrt-ave.f` — `subroutine wrtAVE(step, fnum)`
**Purpose:** Write energy/temperature/pressure averages (11-column format).

#### `wrt-rst.f` — `subroutine wrtRST(step, fnum)`
**Purpose:** Write binary restart file (all 6 derivative levels + energies).

#### `wrt-rst-multi.f`
**Purpose:** Restart with unique filenames per output (no overwrite).

#### `wrt-vel.f` / `wrt-vel-m-ats1.f` — `subroutine wrtVEL(step, fnum)`
**Purpose:** Write velocity data.

#### `wrt-acc.f` — `subroutine wrtACC(step, fnum)`
**Purpose:** Write acceleration data.

#### `wrt-fpt.f` — `subroutine wrtFPT()`
**Purpose:** Write first-passage times for all observed cluster sizes, in nanoseconds. Both Stillinger and Frenkel.

---

### Utilities

#### `kps.f` — Random Number Generator
**Purpose:** Kirkpatrick-Stoll lagged-Fibonacci generator (250 integers, XOR-based).
**Subroutines:**
- `kpsInit(flag)` — Seed from system time (flag=0) or deterministic (flag=1)
- `kpsRND()` → double in [0,1)
- `kpsINT()` → positive 31-bit integer
- `kpsINTf()` → full 32-bit integer (may be negative)
- `kpsMAX(maximum)` → integer in [0, maximum)
- `kpsMAX2()` → scaled integer using preset maximum
- `kpsSetMax(max)` — Set the scaling factor for kpsMAX2

Uses its own COMMON block `/kps/`.

#### `stopdead.f` — `subroutine stopdead(step)`
**Purpose:** Zero total system momentum (single-species version).

#### `stopdead-m.f` — `subroutine stopdead(step)`
**Purpose:** Zero total system momentum (multi-species). Runs 10 iterations of momentum subtraction + thermostat for convergence.

#### `resize.f` — `subroutine resize(step)`
**Purpose:** Gradually shrink or expand box dimensions toward target. Shifts atom positions accordingly. Triggers neighbor list rebuild.

#### `resize2.f`
**Purpose:** Alternative resize implementation.

#### `rates.f` — `subroutine rates(step)`
**Purpose:** Track condensation/evaporation events on the largest cluster. Compares current and previous Stoddard results to count atoms gained/lost.

#### `chg-natom.f` — `subroutine chgNatom(step)`
**Purpose:** Dynamically add or remove atoms. Removal: picks a random monomer. Addition: places at random non-overlapping position. Adjusts neighbor lists and momentum.

#### `chk-nucleation.f` — `subroutine chk_nucleation(step)`
**Purpose:** Auto-abort: if a cluster exceeds `csover` atoms and persists for `clwait` timesteps, terminate the simulation. Saves computation when nucleation has clearly occurred.

---

### Helper/Analysis Programs (`helpers/`)

All helpers are standalone programs that read simulation output files for post-processing.

| File | Program | Purpose |
|------|---------|---------|
| `sizer.f` | `sizer` | Analyze .csi file: extract time series of largest cluster size, temperature, and probability distributions. Supports Stillinger, Frenkel, and old formats. |
| `rdf.f` | `rdf` | Compute radial distribution function g(r) and structure factor S(q) from .xyz trajectory |
| `density.f` | `density` | Compute density profile of largest cluster (radial mass distribution) |
| `densitp.f` | `densitp` | Same as density but with PBC |
| `selclu.f` | `selclu` | Color-code atoms by cluster membership for visualization |
| `selats.f` | `selats` | Select atoms by species for visualization |
| `bars.f` | `bars` | Time-dependent bar chart of cluster size distribution from .csi file |
| `colgas.f` | `colgas` | Track carrier gas exchange between cluster and vapor |
| `bigcl.f` | `bigcl` | Extract and output only the largest cluster |
| `focusclp.f` | `focusclp` | Focus view on cluster with PBC unwrapping |
| `marker.f` | `marker` | Mark specific atoms with different symbols |
| `mesher.f` | `mesher` | Cluster size vs time mesh plot data |
| `press.f` | `press` | Average pressure from .ave file |
| `pressw.f` | `pressw` | Pressure from wall forces (non-PBC only) |
| `thinout.f` | `thinout` | Reduce XYZ file by skipping frames |
| `reverse.f` | `reverse` | Reverse time ordering of XYZ file |
| `veltest.f` | `veltest` | Velocity distribution statistics |
| `sysdens.f` | `sysdens` | System density on a spatial grid |
| `timefile.f` | `timefile` | Generate SciAn time file from .3d |
| `makelog.f` | `makelog` | Generate log files |
| `stotest.f` | `stotest` | Test Stoddard cluster routine on XYZ data |
| `rate-xyz.f` | `rate-xyz` | Compute rates from XYZ trajectories |
| `ratstat.f` | `ratstat` | Statistical analysis of rate data |
| `ratstat-t.f` | `ratstat-t` | Time-resolved rate statistics |
| `passage-time.f` | `passage-time` | Mean first-passage time analysis |
| `mpt.f` | `mpt` | MPT helper utility |
| `tempavg.f` / `tempavg2.f` | `tempavg` | Cluster temperature averaging |
| `meantemp.f` | `meantemp` | Mean temperature computation |
| `t-hist.f` | `t-hist` | Temperature histograms |
| `diffusivity.f` | `diffusivity` | Barrier diffusivity from sizer output |
| `probability-distrib.f` | `probability` | Cluster size probability distributions |

---

## Call Graph

### Main Simulation Path (production: `make frenkel`)

```
program cluster (main.f)
├── kpsInit()                              [kps.f]
├── readdata()                             [read-3d.f]
│   └── fspc()
├── readRST() ─── OR ─── initconf()       [read-rst.f / init-cube1.f]
│                          ├── kpsInit()
│                          ├── kpsRND()
│                          ├── swap()
│                          ├── stopdead()  [stopdead-m.f]
│                          │   └── thermostat()  [scalev-m.f]
│                          ├── neighbor()  [neighbor-pbc-mav.f -DSILENT]
│                          ├── AccAtom()   [for-a-lj-n-pbc-m.f]
│                          ├── AccWall()   [for-w-dummy.f]
│                          └── integrator() × 100 warmup
│
├── stoddard()                             [stoddard-pbc-n-m-fre.f]
│
└── ─── MAIN LOOP ───
    ├── rates()                            [rates.f]
    │   └── stoddard()
    ├── wrtXYZ()                           [wrt-xyz-m-vmd.f]
    ├── wrtVEL()                           [wrt-vel.f]
    ├── wrtACC()                           [wrt-acc.f]
    ├── wrtRST()                           [wrt-rst.f]
    ├── wrtCSI()                           [wrt-csi-fre.f]
    │   └── stoddard()
    ├── wrtAVE()                           [wrt-ave.f]
    ├── thermostat()                       [scalev-m.f]
    ├── resize()                           [resize.f]
    │   └── neighbor()
    ├── chgNatom()                         [chg-natom.f]
    │   ├── stoddard()
    │   ├── kpsMAX()
    │   ├── stopdead()
    │   └── neighbor()
    ├── chk_nucleation()                   [chk-nucleation.f]
    └── integrator()                       [int-verlet-mav.f]
        ├── neighbor()     (conditionally)
        ├── AccAtom()
        └── AccWall()
```

### Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                        .3d INPUT FILE                        │
│  (box size, species, temperature, cutoffs, file paths, ...)  │
└──────────────────────────────┬──────────────────────────────┘
                               │ readdata()
                               ▼
┌──────────────────────────────────────────────────────────────┐
│                     COMMON BLOCKS (Shared State)             │
│                                                              │
│  /ATOMP/  ←→  /ATOMC/  ←→  /ENERG/  ←→  /BOXPP/           │
│  species     positions    energies     box dims              │
│  params      velocities   temperature  PBC flags             │
│  neighbor    accelerations                                   │
│  list                                                        │
│              /CVTPP/  ←→  /TIMEP/  ←→  /DMMPP/             │
│              cluster      timestep     displacement          │
│              stats        params       tracking              │
│                                                              │
│              /FILEP/                                          │
│              file names & intervals                           │
└───────────┬──────────┬──────────┬──────────┬────────────────┘
            │          │          │          │
    ┌───────▼──┐ ┌─────▼────┐ ┌──▼───┐ ┌───▼────┐
    │ AccAtom()│ │integrator│ │stodd-│ │neighbor│
    │ AccWall()│ │    ()    │ │ard() │ │   ()   │
    │          │ │          │ │      │ │        │
    │ reads:   │ │ updates: │ │reads:│ │reads:  │
    │ positions│ │ pos,vel, │ │ pos, │ │ pos,   │
    │ neighbor │ │ acc,ekin │ │ nlist│ │ box    │
    │ list     │ │          │ │      │ │        │
    │ writes:  │ │ calls:   │ │write:│ │writes: │
    │ acc,epot │ │ AccAtom  │ │ cl,  │ │ nlist, │
    └──────────┘ │ AccWall  │ │ cs,  │ │ nidx   │
                 │ neighbor │ │ csi  │ └────────┘
                 └──────────┘ └──────┘
                                 │
                     ┌───────────┼───────────┐
                     ▼           ▼           ▼
              ┌──────────┐ ┌─────────┐ ┌─────────┐
              │ wrtCSI() │ │ rates() │ │chk_nuc()│
              │ .csi/.fsi│ │ .rat    │ │ abort?  │
              └──────────┘ └─────────┘ └─────────┘

                     OUTPUT FILES
    ┌────────┬────────┬────────┬────────┬────────┐
    │  .xyz  │  .vel  │  .rst  │  .ave  │  .fpt  │
    │ coords │ veloc. │restart │ E,T,P  │ FPT    │
    └────────┴────────┴────────┴────────┴────────┘
```

---

## Miscellaneous Files

### Data Files (helpers/)
- `nlist.dat` — Sample neighbor list data for testing
