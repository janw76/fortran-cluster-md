# CLUSTER — Molecular Dynamics Simulation of Vapor-Liquid Nucleation

A Fortran 77 molecular dynamics (MD) code for simulating homogeneous vapor-liquid nucleation of argon and argon/carrier-gas mixtures using Lennard-Jones potentials with periodic boundary conditions.

## Scientific Background

This code was used in the PhD thesis **"Nano-Droplets at Birth — Computer Experiments on Gas Phase Nucleation"** by Jan Wedekind (Universität zu Köln, 2006). The thesis investigated homogeneous vapor-liquid nucleation of argon below the triple point, producing six nucleation rate isotherms based on more than 7500 individual MD simulations in the temperature range 45 K ≤ T ≤ 70 K.

Key scientific contributions enabled by this code:

- **Mean First-Passage Time (MFPT) method** — a new technique for determining nucleation rates, critical cluster sizes, and Zeldovich factors directly from simulation kinetics
- **Thermostat comparison** — demonstrated that simple isokinetic velocity scaling introduces negligible error compared to Nosé-Hoover for argon nucleation
- **Finite-size effects** — showed that previous simulations used systems up to 100× larger than necessary, enabling dramatic computational savings
- **Nucleation rate isotherms** — first-ever rate isotherms from MD simulations, enabling application of the nucleation theorem

### The Physics

Supersaturated argon vapor spontaneously forms liquid droplets through nucleation — the stochastic formation of a critical-sized cluster that can grow into a stable droplet. The simulation:

1. Places argon atoms (and optionally helium carrier gas) in a periodic box
2. Evolves the system via Newtonian mechanics with Lennard-Jones 12-6 potentials
3. Identifies clusters using the **Stillinger criterion** (atoms within a distance cutoff) and optionally the **ten Wolde/Frenkel** definition (based on local coordination number)
4. Monitors cluster size distributions over time to detect nucleation events
5. Records first-passage times for each cluster size to enable MFPT analysis

### References

- Wonczak, S. (2001). *Molekulardynamische Simulationen von Argon-Clustern*. PhD thesis, Universität zu Köln. **(Original CLUSTER code)**
- Wedekind, J. (2006). *Nano-Droplets at Birth — Computer Experiments on Gas Phase Nucleation*. PhD thesis, Universität zu Köln.
- Wedekind, J., Strey, R., & Reguera, D. (2007). New method to analyze simulations of activated processes. *J. Chem. Phys.* **126**, 134103.
- Wedekind, J., Reguera, D., & Strey, R. (2006). Finite-size effects in simulations of nucleation. *J. Chem. Phys.* **125**, 214505.

---

## Code Structure Overview

The code is organized as a modular "construction kit" (Bausatz) where different implementations of each component (forces, integrators, neighbor lists, etc.) can be swapped at compile time via the Makefile.

### Directory Layout

```
cluster/
├── main.f                    # Main program
├── Makefile                  # Build system with ~15 build targets
├── sample.3d                 # Sample input configuration file (v1.30)
├── sample-v1.12.3d           # Older format input file
├── sample-v1.2.3d            # Another input file variant
├── README                    # Original German documentation
│
├── *.inc                     # Common block include files (shared state)
├── init-*.f                  # Initialization routines
├── int-*.f                   # Time integrators
├── for-a-*.f                 # Atom-atom force calculations
├── for-w-*.f                 # Atom-wall force calculations
├── neighbor-*.f              # Neighbor list construction
├── stoddard-*.f              # Cluster identification (Stoddard algorithm)
├── scalev-*.f                # Thermostat (velocity scaling)
├── read-*.f                  # Input routines
├── wrt-*.f                   # Output routines
├── kps.f                     # Random number generator
├── resize.f                  # Box resizing
├── stopdead*.f               # Zero total momentum
├── rates.f                   # Condensation/evaporation rates
├── chg-natom.f               # Dynamic atom count changes
├── chk-nucleation.f          # Nucleation detection & auto-abort
│
└── helpers/                  # Post-processing & analysis tools
    ├── sizer.f               # Cluster size time-series analysis
    ├── rdf.f                 # Radial distribution function
    ├── density.f             # Cluster density profiles
    ├── selclu.f              # Cluster visualization (color-coding)
    ├── diffusivity.f         # Barrier diffusivity analysis
    ├── bars.f                # Cluster size histograms
    ├── colgas.f              # Gas exchange tracking
    ├── ratstat.f             # Rate statistics
    ├── passage-time.f        # Mean first-passage time analysis
    ├── mpt.f                 # MPT helper
    ├── ... (20+ helper programs)
    └── *.inc                 # Copies of include files for helpers
```

### File Categories

| Category | Files | Purpose |
|----------|-------|---------|
| **Main Program** | `main.f` | Simulation loop: init → integrate → output → repeat |
| **Include Files** | `const.inc`, `atomc.inc`, `atomp.inc`, `boxpp.inc`, `cvtpp.inc`, `energ.inc`, `filep.inc`, `timep.inc`, `dmmpp.inc`, `connp.inc` | Shared state via COMMON blocks |
| **Initialization** | `init-cube1.f`, `init-ball1.f`, `init-2atom.f`, `init-2cubes.f`, `init-center.f` | Starting configurations |
| **Integrators** | `int-verlet*.f`, `int-euler*.f`, `int-gear*.f` | Velocity-Verlet, Euler, Gear (3rd–5th order) |
| **Atom Forces** | `for-a-lj*.f`, `for-a-dummy.f` | LJ 12-6 with cutoff, PBC, multi-species |
| **Wall Forces** | `for-w-i9.f`, `for-w-i12.f`, `for-w-iX.f`, `for-w-lj.f`, `for-w-dummy.f` | Repulsive wall potentials |
| **Neighbor Lists** | `neighbor*.f` | Verlet neighbor lists with PBC & adaptive skin |
| **Cluster ID** | `stoddard*.f` | Stillinger + ten Wolde/Frenkel cluster definitions |
| **Thermostat** | `scalev*.f`, `scalev-nh-m.f` | Velocity scaling, Nosé-Hoover |
| **I/O** | `read-3d.f`, `read-rst.f`, `read-xyz.f`, `wrt-*.f` | Input/output routines |
| **Analysis** | `rates.f`, `chk-nucleation.f`, `wrt-fpt.f` | Runtime nucleation analysis |
| **Utilities** | `kps.f`, `resize.f`, `stopdead*.f`, `chg-natom.f` | RNG, box resize, momentum zeroing |

---

## Common Blocks (Shared State)

All inter-module communication happens via Fortran 77 COMMON blocks defined in `.inc` files:

### `const.inc` — Constants & Array Limits
- `mxa = 6000` — Maximum number of atoms
- `mxnlist = mxa²/4` — Maximum neighbor list entries (9,000,000)
- `mxsorts = 10` — Maximum atom species
- Physical constants: Boltzmann (`cikb`), atomic mass unit (`cimu`), Avogadro, Planck, π
- Internal units: length in Å, time in fs, mass in 10⁻²⁷ kg, energy in 10⁻¹⁷ J

### `atomc.inc` — Atomic Coordinates & Derivatives
- Positions: `x0(mxa), y0(mxa), z0(mxa)` (aliased as `x, y, z`)
- Velocities: `x1/y1/z1` (aliased as `vx, vy, vz`)
- Accelerations: `x2/y2/z2` (aliased as `ax, ay, az`)
- Higher derivatives up to 5th order (for Gear integrators): `x3..x5`, `y3..y5`, `z3..z5`

### `atomp.inc` — Atomic Parameters
- `natom` — Total particle count
- `atsorts` — Number of species
- Per-species: symbol, count, mass, σ (LJ diameter), ε (LJ well depth), exponents
- `atnlist(mxnlist)` — Neighbor list array
- `atnidx(mxa)` — Per-atom index into neighbor list
- `att` — System temperature
- `atrcut, atrskin, atskin` — Cutoff, neighbor skin, and Stillinger cluster radii
- `Q, PSI, s` — Nosé-Hoover thermostat variables

### `boxpp.inc` — Box Parameters
- `boxx, boxy, boxz` — Box dimensions
- `bperx, bpery, bperz` — Periodic boundary condition flags per axis
- `bshr` — Box shrink/expand control
- `beng, bpot` — Wall potential parameters

### `cvtpp.inc` — Cluster Statistics
- `cl(mxa)` — Cluster label per atom (Stillinger)
- `clf(mxa)` — Cluster label per atom (ten Wolde/Frenkel)
- `cs(mxa), csf(mxa)` — Cluster sizes
- `csi(mxa), csif(mxa)` — Cluster size distributions
- `cek, cekf, ceks, ceksf` — Cluster kinetic energies
- `cbiggest, cbiggestf` — Largest cluster size
- `cfpt(mxa), cfptf(mxa)` — First passage times per cluster size
- `csover, clwait, clnuctrl` — Nucleation detection parameters

### `energ.inc` — Energies
- `ekin, epot, H` — Kinetic energy, potential energy, Hamiltonian
- `ekinA(mxsorts)` — Kinetic energy per species
- `ewall, ewallsum, fwall, fwallsum` — Wall interaction energies/forces
- `ethst, eclth` — Thermostat control flags

### `filep.inc` — File Parameters
- Filenames and paths for all output files (`.xyz`, `.vel`, `.acc`, `.rst`, `.csi`, `.ave`, `.rat`, `.fpt`)
- Output intervals (`ftxyz`, `ftvel`, etc.)
- `frctrl` — Restart control flag

### `timep.inc` — Time Parameters
- `dt` — Timestep in femtoseconds
- `ntim` — Total number of timesteps
- `btim` — Starting timestep number
- `sctim` — Screen output interval
- `nbtim` — Neighbor list update interval

### `dmmpp.inc` — Displacement/Optimization Statistics
- Moving average of maximum atomic displacements for adaptive neighbor list updates
- `dmmh(0:255)` — History ring buffer
- `dmms, dmmax, dmmav` — Current, maximum, and average displacement sums

---

## Program Flow

### 1. Startup
```
main.f: program cluster
  ├── kpsInit(0)           — Initialize random number generator
  └── readdata()           — Parse .3d input file (read-3d.f)
       ├── Read file paths, output intervals
       ├── Read timestep, box dimensions, PBC flags
       ├── Read thermostat settings, Nosé-Hoover Q
       ├── Read cluster detection parameters
       └── Read per-species: symbol, N, mass, σ, ε, exponents
```

### 2. Initialization
```
  ├── Open output files (xyz, vel, acc, csi, ave, rat, fpt)
  ├── Print simulation parameters to screen
  │
  ├── IF restart file:
  │     └── readRST()      — Read positions/velocities from binary restart
  │
  ├── ELSE (fresh start):
  │     └── initconf()     — Create initial configuration (init-cube1.f)
  │          ├── Place atoms on cubic lattice
  │          ├── Assign random velocities
  │          ├── Swap atom types randomly (for mixtures)
  │          ├── stopdead()    — Zero center-of-mass momentum
  │          ├── neighbor(0)   — Build initial neighbor list
  │          ├── AccAtom()     — Compute atom-atom forces
  │          ├── AccWall()     — Compute atom-wall forces
  │          └── 100 warmup integration steps
  │
  ├── Initialize Nosé-Hoover: PSI, s, H
  └── stoddard(-1)         — Initial cluster identification
```

### 3. Main Simulation Loop
```
  DO i = btim to btim+ntim:
  │
  ├── rates(i)             — Track condensation/evaporation (if enabled)
  ├── wrtXYZ(i)            — Write coordinates (at intervals)
  ├── wrtVEL(i)            — Write velocities
  ├── wrtACC(i)            — Write accelerations
  ├── wrtRST(i)            — Write restart file
  ├── wrtCSI(i)            — Write cluster statistics
  ├── wrtAVE(i)            — Write averages (T, E, pressure)
  │
  ├── Screen output: T_kinetic, E_kin, E_pot, H
  │
  ├── thermostat(i)        — Scale velocities (if enabled)
  ├── resize(i)            — Adjust box size (if enabled)
  ├── chgNatom(i)          — Add/remove atoms (if enabled)
  │
  ├── chk_nucleation(i)    — Check for critical cluster (if enabled)
  │     └── If cluster > csover persists for clwait steps → ABORT
  │
  └── integrator(i)        — Advance one timestep
       ├── Update positions: x += (v + ½a·dt)·dt
       ├── Update velocities (half-step): v += ½a·dt
       ├── Track max displacement for adaptive neighbor list
       ├── IF displacement exceeds skin: neighbor(step)
       ├── Reset accelerations to zero
       ├── AccAtom()        — Recompute atom-atom forces
       ├── AccWall()        — Recompute atom-wall forces
       └── Complete velocity update: v += ½a(new)·dt
```

### 4. Termination
- Simulation ends when all timesteps complete or nucleation detection triggers abort
- Final first-passage times can be written via `wrtFPT()`

---

## Input Files — The `.3d` Configuration Format

The `.3d` files are fixed-format text files controlling all simulation parameters. See `sample.3d` for the latest version (v1.30). Key sections:

| Section | Parameters |
|---------|-----------|
| **Files** | Output directory, filenames for .xyz/.vel/.acc/.rst/.csi/.ave/.rat, output intervals |
| **Time** | Timestep (fs), total steps, start step, screen output interval |
| **Box** | Dimensions (nm), periodicity per axis (X/-), wall potential parameters |
| **Box Resize** | Enable/disable, target size, step size, interval |
| **Restart** | Use restart file (X/-), restart filename |
| **Thermostat** | Type (X=velocity-scaling, N=Nosé-Hoover, -=none), Nosé-Hoover coupling Q |
| **Physics** | Temperature (K), cluster criterion (σ), Frenkel neighbor count, potential cutoff (σ), skin radius (σ), neighbor update interval |
| **Nucleation** | Enable auto-abort (X/-), minimum cluster size, persistence time |
| **Species** | Per species: symbol, atom count, mass (u), σ (nm), ε (K), LJ exponents, thermostat flag, atom-number-change control |

Example (argon + helium carrier gas):
```
atom symbol                            : Ar
number of atoms                        :      125
atom mass                           [u]:       40.0
Lennard-Jones Sigma                [nm]:        0.3405
Lennard-Jones Epsilon               [K]:      120.0
...
atom symbol                            : He
number of atoms                        :      875
atom mass                           [u]:        4.0
Lennard-Jones Sigma                [nm]:        0.258
Lennard-Jones Epsilon               [K]:       10.22
```

---

## Output Files

| Extension | Content | Format |
|-----------|---------|--------|
| `.xyz` | Atomic coordinates at each output step | XYZ format (readable by VMD, SciAn) |
| `.vel` | Atomic velocities | Text, same structure as .xyz |
| `.acc` | Atomic accelerations | Text |
| `.rst` | Restart checkpoint | Fortran unformatted binary |
| `.csi` | Cluster size distribution vs time (Stillinger) | Text: size, count, carrier-gas fraction, kinetic energy, temperature |
| `fsi` | Cluster size distribution vs time (ten Wolde/Frenkel) | Same format as .csi |
| `.ave` | Time series of E_kin, E_pot, E_total, wall energy, temperature, pressure | Text, 11 columns |
| `.rat` | Condensation/evaporation rates of largest cluster | Text |
| `fpt` | First-passage times per cluster size (Stillinger & tW/F) | Text: size, FPT_Stillinger (ns), FPT_Frenkel (ns) |

---

## How to Compile

### Prerequisites
- A Fortran compiler: `gfortran` (recommended), or Intel `ifort`
- GNU `make`

### Quick Start with gfortran

The original Makefile uses `ifort`. To compile with modern `gfortran`:

```bash
cd cluster/

# Edit Makefile: change compiler and flags
# FF= gfortran
# FLAGS= -O3 -std=legacy -w -fallow-argument-mismatch

# Build the main nucleation simulation (PBC, multi-species, carrier gas, Frenkel clusters)
make frenkel

# Or build other variants:
make cluster    # Basic cluster simulation (no PBC)
make pbc        # With periodic boundary conditions
make nuc        # Multi-species + PBC
make silent     # Silent neighbor list updates
make nose       # With Nosé-Hoover thermostat
make diffusion  # Diffusion study variant

# Build analysis tools:
make sizer      # Cluster size time-series analyzer
make rdf        # Radial distribution function
make density    # Density profiles
```

### Build Targets (Simulation)

| Target | Binary | Description |
|--------|--------|-------------|
| `cluster` | `cluster` | Basic: no PBC, single species |
| `pbc` | `cluster` | With PBC, single species |
| `mav` | `cluster` | PBC + adaptive neighbor list |
| `nuc` | `clu-nuc` | PBC + multi-species |
| `silent` | `clu-silent` | Like nuc, silent neighbor updates |
| `frenkel` | `clu-fre` | Full: PBC + multi-species + Frenkel clusters + FPT |
| `nose` | `clu-nose` | With Nosé-Hoover thermostat |
| `diffusion` | `diffusion` | Diffusion coefficient study |

### Running

```bash
./clu-fre sample.3d
```

The program reads the `.3d` file as its sole command-line argument and creates output files in the directory specified within the `.3d` file.

---

## Internal Units

| Quantity | Unit | Symbol |
|----------|------|--------|
| Length | Ångström (10⁻¹⁰ m) | Å |
| Time | Femtosecond (10⁻¹⁵ s) | fs |
| Mass | 10⁻²⁷ kg | — |
| Energy | 10⁻¹⁷ J | — |
| Temperature | Kelvin | K |

Input is in nm and atomic mass units; conversion happens in `read-3d.f`.

---

## Authors & Attribution

The core simulation package **CLUSTER** was developed by **Stefan Wonczak** for his PhD thesis at the Universität zu Köln (2001). **Jan Wedekind** substantially extended the code for his own PhD work (2003–2006), adding the MFPT analysis framework, new thermostats, finite-size optimization, and additional analysis tools. Both worked at the Institute for Physical Chemistry in the group of Prof. Reinhard Strey.

A complete documentation of the original CLUSTER package can be found in Wonczak's thesis:

> Wonczak, S. (2001). *Molekulardynamische Simulationen von Argon-Clustern*. PhD thesis, Universität zu Köln.

The code reflects the practices of its era:
- **Fortran 77** fixed-format (columns 1–72, continuation in column 6)
- **COMMON blocks** for all inter-module data sharing
- **Static arrays** with hardcoded maximum sizes (`mxa=6000`)
- **No dynamic memory allocation**
- **Makefile-based build** with different link combinations for different simulation variants
- Originally compiled with SGI, Cray, Sun, and Intel Fortran compilers
- German comments throughout (the code was written in Cologne)

The simulation campaigns for the thesis ran on clusters at the Universität zu Köln computing center (RRZK), with individual simulations taking up to 2000 ns of simulated time.

---

## See Also

- [`docs/CODE_MAP.md`](docs/CODE_MAP.md) — Detailed file-by-file reference with call graphs
- [`docs/MODERNIZATION.md`](docs/MODERNIZATION.md) — Roadmap for compiling and modernizing the code
