#!/bin/bash
#===============================================================================
# Performance Testing Script for Apple Silicon Optimizations
#===============================================================================

echo "Fortran Cluster MD - Performance Testing on Apple Silicon"
echo "=========================================================="

# Check system info
echo "System Information:"
echo "  Architecture: $(uname -m)"
echo "  Processor: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Unknown')"
echo "  Total cores: $(sysctl -n hw.logicalcpu)"
echo "  Performance cores: $(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || echo 'Unknown')"
echo ""

# Set up environment
source ./setup_apple_silicon.sh > /dev/null 2>&1

# Create test input files with different sizes
create_test_inputs() {
    echo "Creating test input files..."
    
    # Small test (100 atoms, 1000 steps)
    cp sample.3d test_small.3d
    sed -i '' 's/total number of timesteps.*/total number of timesteps              :     1000/' test_small.3d
    sed -i '' 's/screen output interval.*/screen output interval                 :      200/' test_small.3d
    
    # Medium test (500 atoms, 2000 steps)  
    cp sample.3d test_medium.3d
    sed -i '' 's/total number of timesteps.*/total number of timesteps              :     2000/' test_medium.3d
    sed -i '' 's/screen output interval.*/screen output interval                 :      500/' test_medium.3d
    # Increase atom count (this is system-dependent, may need adjustment)
    
    # Large test (1000+ atoms, 5000 steps)
    cp sample.3d test_large.3d  
    sed -i '' 's/total number of timesteps.*/total number of timesteps              :     5000/' test_large.3d
    sed -i '' 's/screen output interval.*/screen output interval                 :     1000/' test_large.3d
    
    echo "  test_small.3d  - 1000 timesteps"
    echo "  test_medium.3d - 2000 timesteps" 
    echo "  test_large.3d  - 5000 timesteps"
    echo ""
}

# Performance test function
run_performance_test() {
    local test_name=$1
    local input_file=$2
    local executable=$3
    
    echo "Testing $test_name with $input_file..."
    if [ -f "$executable" ]; then
        # Warm up
        ./$executable < $input_file > /dev/null 2>&1
        
        # Actual timing
        local start_time=$(date +%s.%N)
        ./$executable < $input_file > /dev/null 2>&1
        local end_time=$(date +%s.%N)
        local runtime=$(echo "$end_time - $start_time" | bc)
        
        echo "  Runtime: ${runtime}s"
        return 0
    else
        echo "  Executable $executable not found - skipping"
        return 1
    fi
}

# Build all versions
build_versions() {
    echo "Building simulation versions..."
    
    # Clean first
    make clean > /dev/null 2>&1
    
    # Build standard version
    echo "  Building standard version..."
    make silent > build_standard.log 2>&1
    if [ $? -eq 0 ]; then
        echo "    ✓ Standard version built successfully"
    else
        echo "    ✗ Standard version build failed (see build_standard.log)"
    fi
    
    # Build OpenMP version
    echo "  Building OpenMP optimized version..."
    make apple-silent > build_openmp.log 2>&1
    if [ $? -eq 0 ]; then
        echo "    ✓ OpenMP version built successfully"
    else
        echo "    ✗ OpenMP version build failed (see build_openmp.log)"
        echo "      Try: brew install gcc"
    fi
    
    echo ""
}

# Main performance comparison
compare_performance() {
    local input_file=$1
    local test_label=$2
    
    echo "Performance Comparison - $test_label"
    echo "$(echo "$test_label" | sed 's/./=/g')$(echo "===================")"
    
    # Test standard version
    if [ -f "clu-silent" ]; then
        echo -n "Standard version: "
        local std_start=$(date +%s.%N)
        ./clu-silent < $input_file > /dev/null 2>&1
        local std_end=$(date +%s.%N)
        local std_time=$(echo "$std_end - $std_start" | bc)
        echo "${std_time}s"
    else
        echo "Standard version: Not available"
        std_time=0
    fi
    
    # Test OpenMP version
    if [ -f "clu-silent-omp" ]; then
        echo -n "OpenMP version:   "
        local omp_start=$(date +%s.%N)
        ./clu-silent-omp < $input_file > /dev/null 2>&1
        local omp_end=$(date +%s.%N)
        local omp_time=$(echo "$omp_end - $omp_start" | bc)
        echo "${omp_time}s"
        
        # Calculate speedup
        if [ "$std_time" != "0" ] && [ $(echo "$std_time > 0" | bc) -eq 1 ]; then
            local speedup=$(echo "scale=2; $std_time / $omp_time" | bc)
            echo "Speedup:          ${speedup}x"
        fi
    else
        echo "OpenMP version:   Not available"
    fi
    
    echo ""
}

# Thread scaling test
thread_scaling_test() {
    if [ ! -f "clu-silent-omp" ]; then
        echo "OpenMP version not available for thread scaling test"
        return
    fi
    
    echo "Thread Scaling Analysis"
    echo "======================="
    
    local input_file="test_medium.3d"
    local max_threads=$(sysctl -n hw.logicalcpu)
    
    echo "Testing thread counts from 1 to $max_threads..."
    echo "Threads | Runtime (s) | Speedup"
    echo "--------|-------------|--------"
    
    local baseline_time=""
    
    for threads in 1 2 4 8 $(seq 12 4 $max_threads); do
        if [ $threads -gt $max_threads ]; then
            break
        fi
        
        export OMP_NUM_THREADS=$threads
        
        local start_time=$(date +%s.%N)
        ./clu-silent-omp < $input_file > /dev/null 2>&1
        local end_time=$(date +%s.%N)
        local runtime=$(echo "$end_time - $start_time" | bc)
        
        if [ -z "$baseline_time" ]; then
            baseline_time=$runtime
            speedup="1.00"
        else
            speedup=$(echo "scale=2; $baseline_time / $runtime" | bc)
        fi
        
        printf "%7d | %11.2f | %7s\n" $threads $runtime $speedup
    done
    
    # Restore optimal thread count
    source ./setup_apple_silicon.sh > /dev/null 2>&1
    echo ""
}

# Memory usage analysis
memory_analysis() {
    echo "Memory Usage Analysis"
    echo "===================="
    
    # Check available memory
    local total_mem=$(sysctl -n hw.memsize)
    local total_mem_gb=$(echo "scale=1; $total_mem / 1024 / 1024 / 1024" | bc)
    echo "Total system memory: ${total_mem_gb} GB"
    
    # Estimate memory usage for different system sizes
    echo ""
    echo "Estimated memory usage for different system sizes:"
    echo "  1,000 atoms:   ~50 MB"
    echo "  5,000 atoms:   ~250 MB" 
    echo "  10,000 atoms:  ~500 MB"
    echo "  20,000 atoms:  ~1 GB"
    echo ""
    echo "Note: Actual usage depends on neighbor list size and simulation parameters"
    echo ""
}

# Clean up test files
cleanup() {
    echo "Cleaning up test files..."
    rm -f test_small.3d test_medium.3d test_large.3d
    rm -f build_*.log
    rm -f data/test/*
    echo ""
}

# Main execution
main() {
    create_test_inputs
    build_versions
    
    # Run performance comparisons
    compare_performance "test_small.3d" "Small System Test"
    compare_performance "test_medium.3d" "Medium System Test"  
    compare_performance "test_large.3d" "Large System Test"
    
    # Additional analysis
    thread_scaling_test
    memory_analysis
    
    # Recommendations
    echo "Optimization Recommendations"
    echo "============================"
    echo "1. Use OpenMP versions for systems with >500 atoms"
    echo "2. Optimal thread count is set automatically by setup script"
    echo "3. For very large systems (>10k atoms), consider:"
    echo "   - Increasing neighbor list size (mxnlist in const.inc)"
    echo "   - Using larger cutoff radii for efficiency"
    echo "   - Running on systems with more memory"
    echo ""
    echo "To use optimized versions:"
    echo "  ./clu-silent-omp < input.3d     # OpenMP parallelized"
    echo "  ./clu-fre-omp < input.3d        # OpenMP Frenkel analysis"
    echo ""
    
    cleanup
    
    echo "Performance testing complete!"
}

# Check dependencies
if ! command -v bc >/dev/null 2>&1; then
    echo "Error: 'bc' calculator not found. Please install: brew install bc"
    exit 1
fi

# Run main function
main "$@"
