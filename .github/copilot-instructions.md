# Copilot Instructions: Fortran Cluster MD

## Project Overview
This is a Fortran 77/90 molecular dynamics simulation code for studying Lennard-Jones atomic clusters. The code uses a modular architecture with specialized components for different simulation variants and analysis tools.

## Architecture & Key Components

### Core Simulation Flow
- **Main program**: `cluster/main.f` - Simulation loop with configurable I/O intervals
- **Input files**: `.3d` files (e.g., `sample.3d`) define all simulation parameters in fixed-format sections
- **Integration algorithms**: Multiple time-stepping methods via `int-*.f` files (Verlet, Gear 3rd-5th order, Euler)
- **Force calculations**: Modular force routines (`for-a-*.f` for atom-atom, `for-w-*.f` for walls)
- **Boundary conditions**: PBC (periodic) and non-PBC variants with neighbor lists for optimization

### Critical File Naming Conventions
The code uses systematic naming patterns that determine functionality:
- `*-pbc*`: Periodic boundary conditions enabled
- `*-n*`: Neighbor list optimization 
- `*-m*`: Multi-species support
- `*-mav*`: Modified neighbor list with average displacement tracking
- `*-silent*`: Reduced output for production runs
- `*-fre*`: Frenkel-type analysis features
- `*-nh*`: Nose-Hoover thermostat integration

### Build System Architecture
The `Makefile` defines multiple simulation variants as separate targets:
- `cluster`: Basic simulation
- `pbc`, `mav`: Boundary condition variants  
- `silent`, `frenkel`, `nose`: Specialized analysis modes
- Helper utilities built to `$(HOME)/bin/` with dependencies

### Data Structures (via Common Blocks)
- `const.inc`: Physical constants and array size limits (`mxa=6000` max atoms)
- `atomp.inc`: Particle parameters (masses, LJ parameters, species data)
- `atomc.inc`: Coordinate arrays with time derivatives for high-order integrators
- `cvtpp.inc`: Cluster analysis data (size distributions, connectivity)
- `boxpp.inc`: Simulation box and boundary conditions
- `energ.inc`: Energy components and thermodynamic quantities

## Development Workflows

### Building Simulations
```bash
cd cluster/
make clean
make <variant>  # e.g., make silent, make frenkel

# Apple Silicon optimized builds with OpenMP
./setup_apple_silicon.sh  # Configure environment
make apple-silent          # Parallel silent mode
make apple-frenkel         # Parallel Frenkel analysis
make benchmark            # Build both versions for comparison
```

### Performance Optimization for Apple Silicon
The codebase includes OpenMP-parallelized versions optimized for modern multi-core processors:
- `*-omp.f` files: OpenMP parallelized force calculations, neighbor lists, and integrators
- Automatic thread configuration based on performance/efficiency core detection
- SIMD vectorization hints for Apple Silicon's wide vector units
- Memory-efficient neighbor list construction with reduced cache misses

Run `./perf_test.sh` to compare serial vs parallel performance.

### Running Simulations
```bash
./cluster < sample.3d
# Or specific variants:
./clu-silent < input.3d
./clu-frenkel < nucleation.3d
```

### Analysis Pipeline
1. Main simulation outputs: `.xyz` (trajectories), `.csi` (cluster sizes), `.ave` (averages)
2. Build analysis tools: `make sizer`, `make rdf`, `make density`, etc.
3. Process with helpers: `~/bin/sizer input.3d > analysis.dat`

### Parameter File Structure (`.3d`)
Files use fixed-column format with version-dependent sections:
- File paths and output intervals
- Time parameters (timestep, total steps)  
- Box geometry and boundary conditions
- Thermostat and restart settings
- Species definitions (LJ parameters, masses)

## Project-Specific Conventions

### Integration Method Selection
Choose integrator by linking appropriate `int-*.f` file in Makefile target:
- `int-verlet-mav.f`: Standard production choice (neighbor lists + displacement tracking)
- `int-gear*.f`: Higher-order accuracy for sensitive systems
- `int-verlet-nh-mav.f`: Nose-Hoover thermostat (single species only)

### Cluster Analysis Workflow
The code specializes in nucleation studies:
1. Use `stoddard*.f` for real-time cluster identification
2. Configure `clnuctrl='X'` to auto-terminate when clusters exceed threshold
3. Post-process with `sizer`, `passage-time`, `tempavg` for statistics

### Helper Utility Patterns
Analysis programs follow consistent argument patterns:
```bash
~/bin/sizer input.3d switch > output.dat
~/bin/tempavg sizer.dat output.dat startsize endsize  
~/bin/mpt passage-times.dat output.dat startsize endsize
```

## Key Integration Points

### Force Calculation Modules
Forces are computed via function calls: `AccAtom()` for inter-particle forces, `AccWall()` for boundaries. These call the linked `for-*.f` implementations.

### Neighbor List Management  
Neighbor lists (`neighbor*.f`) are rebuilt based on displacement criteria in `dmmpp.inc`. The `*-mav` variants use average displacement for better performance.

### Thermostat Integration
Velocity scaling and Nose-Hoover coupling requires careful integration with the time-stepping scheme. Only single-species systems support Nose-Hoover.

## Common Patterns

When adding new analysis: Follow helper program structure with `readdata()` call, argument parsing via `getarg()`, and fixed-format output for downstream processing.

When modifying integrators: Maintain the `integrator(step)` interface and update both position/velocity arrays and energy calculations.

When debugging simulations: Check `.ave` output file for energy conservation and use `ftave` parameter to control output frequency.

When optimizing for modern hardware: Use the OpenMP-parallelized variants (`*-omp.f`) which include vectorization hints and cache-friendly memory access patterns optimized for Apple Silicon processors.
