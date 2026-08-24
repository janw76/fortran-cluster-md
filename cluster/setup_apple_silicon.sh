#!/bin/bash
#===============================================================================
# Apple Silicon Optimization Script for Fortran Cluster MD
#===============================================================================
# This script sets up the optimal environment for running on Apple Silicon

echo "Setting up Apple Silicon optimizations for Fortran Cluster MD..."

# Check if we're on Apple Silicon
if [[ $(uname -m) != "arm64" ]]; then
    echo "Warning: This script is optimized for Apple Silicon (arm64) processors"
    echo "Current architecture: $(uname -m)"
fi

# Set optimal thread counts for Apple Silicon
# M1/M2/M3 have performance and efficiency cores
get_optimal_threads() {
    # Get total logical cores
    total_cores=$(sysctl -n hw.logicalcpu)
    # Get performance cores (approximate)
    perf_cores=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || echo $((total_cores / 2)))
    
    echo "Detected cores: $total_cores total, $perf_cores performance"
    
    # For MD simulations, use performance cores + some efficiency cores
    # but leave some headroom for system
    optimal_threads=$((perf_cores + perf_cores / 2))
    if [ $optimal_threads -gt $((total_cores - 2)) ]; then
        optimal_threads=$((total_cores - 2))
    fi
    
    echo $optimal_threads
}

optimal_threads=$(get_optimal_threads)

# Set OpenMP environment variables
export OMP_NUM_THREADS=$optimal_threads
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_SCHEDULE=guided,32

# Apple Silicon specific optimizations
export OMP_THREAD_LIMIT=$optimal_threads
export OMP_MAX_ACTIVE_LEVELS=1

# Memory optimizations
export OMP_STACKSIZE=64M

# Set compiler paths if needed
if command -v gfortran >/dev/null 2>&1; then
    export APPLE_FF=gfortran
elif command -v gfortran-13 >/dev/null 2>&1; then
    export APPLE_FF=gfortran-13
elif command -v gfortran-12 >/dev/null 2>&1; then
    export APPLE_FF=gfortran-12
else
    echo "Warning: No gfortran found. Please install via: brew install gcc"
fi

echo "OpenMP Configuration:"
echo "  Threads: $OMP_NUM_THREADS"
echo "  Binding: $OMP_PROC_BIND"
echo "  Places: $OMP_PLACES"
echo "  Schedule: $OMP_SCHEDULE"

# Create a simple performance test
create_perf_test() {
    cat > perf_test.sh << 'EOF'
#!/bin/bash
echo "Performance comparison test"
echo "=========================="

if [ ! -f "sample.3d" ]; then
    echo "Error: sample.3d not found. Please ensure you're in the cluster directory."
    exit 1
fi

# Create a short test input
cp sample.3d test_perf.3d
sed -i '' 's/total number of timesteps.*/total number of timesteps              :     1000/' test_perf.3d
sed -i '' 's/screen output interval.*/screen output interval                 :      100/' test_perf.3d

echo "Building optimized version..."
make clean > /dev/null 2>&1
make apple-silent > /dev/null 2>&1

if [ -f "clu-silent-omp" ]; then
    echo "Running OpenMP optimized version..."
    time ./clu-silent-omp < test_perf.3d > /dev/null 2>&1
    echo ""
fi

echo "Building standard version..."
make clean > /dev/null 2>&1  
make silent > /dev/null 2>&1

if [ -f "clu-silent" ]; then
    echo "Running standard version..."
    time ./clu-silent < test_perf.3d > /dev/null 2>&1
fi

rm -f test_perf.3d
echo ""
echo "Test complete. Compare the 'real' times above."
EOF
    chmod +x perf_test.sh
}

create_perf_test

echo ""
echo "Setup complete! To build optimized versions:"
echo "  make apple-silent    # OpenMP parallelized silent mode"
echo "  make apple-frenkel   # OpenMP parallelized Frenkel analysis"
echo "  make benchmark       # Build both standard and optimized versions"
echo ""
echo "To run performance comparison:"
echo "  ./perf_test.sh"
echo ""
echo "For best performance, ensure your input uses large enough systems"
echo "(>1000 atoms) to benefit from parallelization."
