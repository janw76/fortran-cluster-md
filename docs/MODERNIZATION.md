# MODERNIZATION.md — Revival and Modernization Roadmap

## Current State Assessment

### What We Have
- **~100 Fortran 77 source files** (fixed-format, columns 1–72)
- **10 COMMON block include files** for shared state
- **Static arrays** with hardcoded `mxa=6000` limit (9M entries in neighbor list)
- **No dynamic memory allocation**
- **Makefile-based build** with Intel `ifort` as default compiler
- **German comments** mixed with English
- **No version control integration** (archive target creates `.tar.gz`)
- **No unit tests** (just `test.f` = hello world)
- **No documentation** beyond the original German README (now rectified)

### What Works Well
- **Modular "construction kit" design** — swapping integrators, force routines, thermostats at link time is clever
- **Complete, self-contained** — no external library dependencies
- **Well-commented common blocks** — each variable documented
- **Production-tested** — 7500+ simulations for the thesis ran successfully
- **Physically correct** — results match theoretical predictions

### What's Dated
- COMMON blocks instead of modules
- `GOTO`-based control flow (labels 10, 20, 30...)
- `EQUIVALENCE` statements for aliasing
- Implicit typing still partially relied upon (despite `implicit none`)
- Fixed-format source (column 6 continuations, comment in column 1)
- No error handling beyond `write(*,*) ... stop`
- Hardcoded array sizes limit system size to 6000 atoms

---

## Phase 1: Compile on Modern gfortran (Minimal Changes)

### Goal: Get it running with zero physics changes

**Estimated effort: 1–2 hours**

#### Step 1: Modify the Makefile

```makefile
# Replace:
FF= ifort
FLAGS= -O3 -warn -traceback

# With:
FF= gfortran
FLAGS= -O3 -std=legacy -w -fallow-argument-mismatch -ffixed-form
```

Key flags:
- `-std=legacy` — Accept all F77 constructs without warnings
- `-fallow-argument-mismatch` — Needed for the `EQUIVALENCE` and implicit type coercions
- `-ffixed-form` — Force fixed-format parsing (should be default for `.f` files)
- `-w` — Suppress warnings (there will be hundreds)

#### Step 2: Fix Known Issues

1. **`fdate()` in `kps.f`** — Non-standard intrinsic. Replace with:
   ```fortran
   character*24 idate
   call fdate(idate)
   ```
   gfortran supports `fdate` as an extension. If not available, replace with `date_and_time()`:
   ```fortran
   character*8 date
   character*10 time
   call date_and_time(date, time)
   read(time, '(I2)') jh
   read(time(3:4), '(I2)') jm
   read(time(5:6), '(I2)') js
   ```

2. **`ieor()` and `iand()` in `kps.f`** — Standard intrinsics, should work fine.

3. **`getarg()` in `read-3d.f`** — gfortran supports this as an extension. For strict compliance, use `get_command_argument()`.

4. **`connp.inc`** — Has incomplete array declarations (`cnvibd(), cnbena(), cntora()`) that will cause compilation errors. Since this file is only included by unused helpers, either:
   - Fix the declarations
   - Remove the include from any files that use it

5. **Integer overflow in `mxnlist`** — `mxa*mxa/4 = 6000*6000/4 = 9,000,000`. This is fine for 32-bit integers but uses ~36 MB for the neighbor list alone. On modern systems this is nothing.

#### Step 3: Test Compilation

```bash
cd cluster/
# Edit Makefile as above, then:
make clean
make frenkel    # The most feature-complete target
./clu-fre sample.3d
```

#### Expected Issues
- Some warnings about unused variables (harmless)
- Possible type mismatch warnings in `write` format statements
- The `connp.inc` issue if any file includes it

---

## Phase 2: Modernize to Fortran 90/95 (Functional Upgrade)

### Goal: Modern Fortran while preserving physics

**Estimated effort: 2–4 weeks**

### 2a. Free-Form Source

Convert all `.f` files to `.f90` (free-form):
- Remove column restrictions
- Replace `C` and `*` comments with `!`
- Remove column-6 continuation characters, use `&`
- Tool: [`f77_to_f90`](https://github.com/pjb7687/f77_to_f90) or manual with regex

```bash
# Automated first pass (still needs manual cleanup):
for f in *.f; do
  sed 's/^[cC*]/!/; s/     \$/\&/' "$f" > "${f%.f}.f90"
done
```

### 2b. Replace COMMON Blocks with Modules

Create Fortran 90 modules for each `.inc` file:

```fortran
! module_atomc.f90 — replaces atomc.inc
module atomc_mod
  use const_mod, only: mxa
  implicit none
  
  double precision :: x(mxa), y(mxa), z(mxa)
  double precision :: vx(mxa), vy(mxa), vz(mxa)
  double precision :: ax(mxa), ay(mxa), az(mxa)
  ! ... higher derivatives for Gear integrators
  
end module atomc_mod
```

Replace `include 'atomc.inc'` with `use atomc_mod`.

### 2c. Replace Static Arrays with Allocatable

```fortran
! Before (const.inc):
parameter (mxa=6000)

! After (module_const.f90):
module const_mod
  implicit none
  integer :: mxa = 6000          ! Now adjustable at runtime
  integer :: mxnlist             ! Computed from mxa
end module

! In atomc_mod:
double precision, allocatable :: x(:), y(:), z(:)
! Allocated after reading natom from .3d file
```

This removes the 6000-atom limit and reduces memory for small systems.

### 2d. Replace GOTO with Structured Control

```fortran
! Before:
10  continue
    ...
    if (i.le.ntim+btim) goto 10

! After:
do while (i <= ntim + btim)
    ...
    i = i + 1
end do
```

### 2e. Remove EQUIVALENCE

Replace `equivalence (x, x0)` with just using one set of names, or explicit assignment.

### 2f. Encapsulate Subroutines in Modules

Group related routines into modules:
- `forces_mod` — AccAtom, AccWall
- `integrator_mod` — integrator
- `neighbor_mod` — neighbor list construction
- `cluster_mod` — Stoddard algorithm
- `io_mod` — all read/write routines
- `thermostat_mod` — velocity scaling, Nosé-Hoover

### 2g. Add a Proper Build System

Replace the Makefile with CMake:

```cmake
cmake_minimum_required(VERSION 3.12)
project(cluster LANGUAGES Fortran)

# Select simulation variant via CMake option
option(USE_PBC "Enable periodic boundary conditions" ON)
option(USE_FRENKEL "Enable ten Wolde/Frenkel clusters" ON)
option(USE_NOSE_HOOVER "Use Nosé-Hoover thermostat" OFF)

# Sources based on options
set(SOURCES main.f90 kps.f90 read_3d.f90 ...)
if(USE_FRENKEL)
  list(APPEND SOURCES stoddard_pbc_n_m_fre.f90)
else()
  list(APPEND SOURCES stoddard_pbc_n.f90)
endif()

add_executable(cluster ${SOURCES})
target_compile_options(cluster PRIVATE -O3)
```

---

## Phase 3: Parallelization

### 3a. OpenMP — CPU Multi-Threading

**Best opportunities (in order of impact):**

1. **Force calculation (`AccAtom`)** — The inner loop over neighbor pairs is embarrassingly parallel:
   ```fortran
   !$omp parallel do private(j,dx,dy,dz,rij2,rij6,fij) &
   !$omp& reduction(+:esum) schedule(dynamic)
   do k = 1, num_pairs
     ! ... force calculation ...
   end do
   !$omp end parallel do
   ```
   **Caveat:** The current neighbor list encoding (negative indices as atom markers) complicates parallelization. Restructure to a pair list: `pairs(2, npairs)`.

2. **Neighbor list construction** — O(N²) pair distance checks are parallelizable:
   ```fortran
   !$omp parallel do private(j,dx,dy,dz,r2) schedule(dynamic)
   do i = 1, natom
     ! build neighbors for atom i
   end do
   ```
   **Caveat:** Need thread-local list buffers merged afterwards.

3. **Stoddard cluster search** — The initial distance/neighbor computation is parallelizable; the flood-fill is inherently serial but fast (O(N)).

4. **Kinetic energy computation** — Simple `!$omp parallel do reduction(+:ekin)`.

**Expected speedup:** 4–8× on 8 cores for N > 1000 atoms.

### 3b. GPU Acceleration

#### OpenACC (Easiest Path)

```fortran
!$acc parallel loop gang vector reduction(+:esum)
do k = 1, num_pairs
  ! force calculation
end do
```

**Requires:** NVIDIA GPU + PGI/NVHPC compiler. **Not directly usable on Apple M1 GPU.**

#### Metal / Apple GPU

Apple M1 Max has 32 GPU cores and unified memory. Options:

1. **Metal Compute Shaders** — Write force kernel in Metal Shading Language, call from Fortran via C bridge:
   ```
   Fortran → C wrapper → Objective-C → Metal API → GPU kernel
   ```
   Significant effort but full access to M1 GPU.

2. **OpenCL** (deprecated on macOS but still functional):
   ```
   Fortran → C wrapper → OpenCL API → GPU kernel
   ```

3. **Accelerate framework (vDSP/BLAS)** — For vectorized math operations, useful for the force accumulation.

**Realistically:** GPU acceleration is overkill for 6000 atoms. GPUs shine at N > 100,000. For this code's typical system sizes (100–3000 atoms), CPU optimization with OpenMP is more practical.

#### For Larger Systems (Future)

If scaling beyond 6000 atoms:
- **CUDA** (via NVHPC) for NVIDIA GPUs
- **HIP** for AMD GPUs
- **SYCL** for portable GPU code

### 3c. MPI — Domain Decomposition

**Assessment: Probably overkill.** The thesis optimized for small systems (as few as 32 argon atoms). MPI domain decomposition introduces significant complexity and communication overhead that only pays off for N > 10,000 atoms with large boxes.

If needed: standard spatial decomposition — divide box into subdomains, each MPI rank handles one, exchange boundary atoms via halo regions.

---

## Phase 4: Apple M1 Max — Specific Considerations

### Architecture
- **CPU:** 10-core (8 performance + 2 efficiency), ARM64 (aarch64)
- **GPU:** 32-core Apple GPU (Metal only, no CUDA/OpenCL)
- **Memory:** Unified 32/64 GB (no CPU↔GPU transfer overhead!)
- **Neural Engine:** 16-core (not useful for MD)

### Compiler Options

```bash
# Homebrew gfortran (GCC):
brew install gcc
gfortran-14 -O3 -march=native -mcpu=apple-m1 -std=legacy ...

# Flang (LLVM Fortran, ARM-native):
brew install llvm
flang-new -O3 -march=armv8.5-a ...
```

### Performance vs 2006 Hardware

The thesis computations ran on:
- Intel Pentium 4 / Xeon clusters (~2–3 GHz, single core)
- RRZK computing center (exact specs unknown, likely 2004-era Xeons)

**Single-core comparison:**
| Metric | ~2005 Xeon | Apple M1 Max P-core | Speedup |
|--------|-----------|---------------------|---------|
| Clock | 3.0 GHz | 3.2 GHz | 1.1× |
| IPC | ~1.0 | ~8.0 | 8× |
| SIMD width | SSE2 (128-bit) | NEON (128-bit) | 1× |
| Memory BW | ~6 GB/s | ~200 GB/s | 33× |
| L1 cache | 32 KB | 192 KB | 6× |
| **Est. single-core** | | | **5–10×** |

**Multi-core:**
| | 2005 cluster node | Apple M1 Max | Speedup |
|--|-------------------|-------------|---------|
| Cores used | 1–2 (no OpenMP in code) | 8 performance | 4–8× |
| **Est. total** | | | **20–80×** |

**Practical impact:** A simulation that took 48 hours in 2006 would complete in **~1–4 hours** on M1 Max. The entire thesis campaign (7500 simulations) that likely required months of cluster time could be completed in **days to a week** on a single M1 Max.

### Unified Memory Advantage

The M1's unified memory architecture eliminates CPU↔GPU data transfer — the dominant bottleneck in GPU-accelerated MD on discrete GPUs. If Metal GPU kernels were implemented, the force calculation could access the same memory as the CPU code with zero copy overhead. This makes GPU acceleration more attractive for smaller N than on traditional GPU architectures.

---

## Phase 5: Alternative Approaches

### 5a. Python + C/Fortran Kernels

Rewrite the simulation control in Python, keep performance-critical loops in Fortran/C:

```python
# simulation.py
import numpy as np
from cluster_forces import compute_lj_forces  # Compiled Fortran/C
from cluster_neighbor import build_neighbor_list

class MDSimulation:
    def __init__(self, config_file):
        self.config = parse_3d(config_file)
        self.positions = np.zeros((self.config.natom, 3))
        # ...
    
    def step(self):
        forces = compute_lj_forces(self.positions, self.neighbor_list, ...)
        # Velocity-Verlet in numpy (fast enough for N < 10000)
        self.velocities += 0.5 * forces / self.masses * self.dt
        self.positions += self.velocities * self.dt
        # ...
```

**Pros:** Modern tooling, easy visualization (matplotlib, nglview), simple parallelization, easy to extend
**Cons:** 10–100× slower for inner loops vs compiled Fortran; needs ctypes/f2py for kernels

### 5b. Julia Port

Julia is increasingly popular in computational physics:

```julia
# cluster.jl
struct LJParams
    σ::Float64
    ε::Float64
    rcut::Float64
end

function compute_forces!(forces, positions, neighbors, params::LJParams)
    @inbounds @simd for k in 1:length(neighbors)
        # LJ force calculation
    end
end
```

**Pros:** Near-Fortran speed, modern syntax, built-in parallelism, excellent scientific ecosystem
**Cons:** Smaller community than Python, JIT compilation warmup

### 5c. Use an Existing MD Framework

For production-quality nucleation simulations today:
- **LAMMPS** — Industry-standard MD code, supports LJ, PBC, all thermostats, GPU acceleration, MPI. Would need custom nucleation analysis.
- **GROMACS** — Primarily for biomolecules but handles LJ systems.
- **OpenMM** — Python-friendly, GPU-accelerated.

The unique value of *this* code is the integrated nucleation analysis (Stoddard, Frenkel clusters, MFPT, first-passage times). These would need to be reimplemented as post-processing or plugins.

---

## Step-by-Step Modernization Plan

### Phase 1: Just Make It Run (1–2 hours)
- [ ] Modify Makefile for gfortran
- [ ] Fix any compilation errors (connp.inc, fdate)
- [ ] Compile `frenkel` target
- [ ] Run with `sample.3d`, verify output matches expectations
- [ ] Git commit: "Build with modern gfortran"

### Phase 2: Free-Form Conversion (1 week)
- [ ] Convert all `.f` → `.f90` (automated + manual cleanup)
- [ ] Verify compilation and output unchanged
- [ ] Git commit: "Convert to free-form Fortran 90"

### Phase 3: Module Architecture (2 weeks)
- [ ] Create modules from include files
- [ ] Replace COMMON blocks with `use` statements
- [ ] Replace EQUIVALENCE with explicit naming
- [ ] Replace GOTO with structured control flow
- [ ] Replace static arrays with allocatable
- [ ] Verify output bit-for-bit identical
- [ ] Git commit: "Replace COMMON blocks with F90 modules"

### Phase 4: Build System (1 day)
- [ ] Write CMakeLists.txt
- [ ] Configure build variants via CMake options
- [ ] Add install targets for helper programs
- [ ] Git commit: "CMake build system"

### Phase 5: OpenMP Parallelization (1 week)
- [ ] Restructure neighbor list for parallel access
- [ ] Add OpenMP to force calculation
- [ ] Add OpenMP to neighbor list build
- [ ] Add OpenMP to energy computation
- [ ] Benchmark: measure scaling 1–8 cores
- [ ] Git commit: "OpenMP parallelization"

### Phase 6: Testing & Validation (ongoing)
- [ ] Create reference output from original code
- [ ] Write validation script comparing new vs reference
- [ ] Add CI (GitHub Actions with gfortran)
- [ ] Document any numerical differences and their causes

### Phase 7: Enhanced Analysis (optional)
- [ ] Python wrapper for output file parsing
- [ ] Jupyter notebook for MFPT analysis
- [ ] Modern visualization pipeline (VMD scripts or nglview)

---

## Appendix: Compiler Compatibility Matrix

| Compiler | Platform | Notes |
|----------|----------|-------|
| gfortran 12+ | Linux, macOS (Homebrew), Windows (MinGW) | Recommended. Use `-std=legacy` |
| Intel ifort | Linux, macOS (x86 only) | Original compiler. Being replaced by ifx |
| Intel ifx | Linux | oneAPI Fortran compiler, successor to ifort |
| NVIDIA nvfortran | Linux | Needed for OpenACC GPU offloading |
| flang-new (LLVM) | Linux, macOS (ARM native) | Maturing rapidly, good for M1 |
| NAG Fortran | Linux, macOS | Strictest standard compliance checking |
