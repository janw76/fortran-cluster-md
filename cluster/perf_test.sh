#!/bin/bash
# Performance comparison: clu-silent vs clu-silent-omp across 3 input sizes
set -u

export OMP_NUM_THREADS=$(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || sysctl -n hw.logicalcpu)
export OMP_PROC_BIND=true
export OMP_PLACES=cores

cd "$(dirname "$0")"
trap 'rm -f test_small.3d test_medium.3d test_large.3d' EXIT

for size in small:1000:200 medium:2000:500 large:5000:1000; do
    IFS=: read -r name steps interval <<< "$size"
    f="test_${name}.3d"
    cp sample.3d "$f"
    sed -i '' "s/total number of timesteps.*/total number of timesteps              :     ${steps}/" "$f"
    sed -i '' "s/screen output interval.*/screen output interval                 :      ${interval}/" "$f"
done

printf "%-8s %-16s %-10s\n" "Size" "Executable" "Time(s)"
for size in small medium large; do
    f="test_${size}.3d"
    for exe in clu-silent clu-silent-omp; do
        if [ ! -x "./$exe" ]; then
            echo "$size $exe: skipped (binary not found)"
            continue
        fi
        t=$( { /usr/bin/time -p ./"$exe" < "$f" > /dev/null; } 2>&1 | awk '/real/{print $2}' )
        printf "%-8s %-16s %-10s\n" "$size" "$exe" "$t"
    done
done

exit 0
