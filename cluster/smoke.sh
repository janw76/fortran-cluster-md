#!/bin/bash
#-------------------------------------------------------------------------------
# smoke.sh -- regression gate: build + run every main sim target with a tiny
# 50-timestep input and assert exit 0 + at least one finite energy line.
# bash 3.2 safe (macOS /bin/bash): no associative arrays, no `<<<` here-string
# reliance beyond what 3.2 already supports.
#-------------------------------------------------------------------------------
set -u
cd "$(dirname "$0")" || exit 1

export PATH=/opt/homebrew/bin:$PATH
export OMP_NUM_THREADS=1

TARGETS="cluster:cluster pbc:cluster mav:cluster mul:cluster-m nuc:clu-nuc \
silent:clu-silent frenkel:clu-fre diffusion:diffusion \
moviefrenkel:clu-fre-mov nose:clu-nose rst:clu-rst test:clutest \
apple-silent:clu-silent-omp apple-frenkel:clu-fre-omp"

WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"; make clean >/dev/null 2>&1' EXIT

RESULTS=""
OVERALL=0

for pair in $TARGETS; do
  target=${pair%%:*}
  bin=${pair##*:}

  make "$target" >"$WORKDIR/build-$target.log" 2>&1
  if [ $? -ne 0 ] || [ ! -x "./$bin" ]; then
    echo "FAIL  $target  (build failed; see $WORKDIR/build-$target.log)"
    RESULTS="$RESULTS
FAIL  $target"
    OVERALL=1
    continue
  fi

  rundir="$WORKDIR/$target"
  mkdir -p "$rundir/data/test"
  input="$rundir/tiny.3d"
  cp sample.3d "$input"
  sed -i '' "s/^total number of timesteps.*/total number of timesteps              :       50/" "$input"
  sed -i '' "s/^screen output interval.*/screen output interval                 :       10/" "$input"
  sed -i '' "s/^pure random velocities after init.*/pure random velocities after init      : -/" "$input"
  # rst (wrt-rst-multi) requires a '####' sequence placeholder in the restart
  # filename, or it aborts (bare `stop`, exit 0) before any screen output.
  if [ "$target" = "rst" ]; then
    sed -i '' "s|^filename for restart-file.*|filename for restart-file              : test-####.rst|" "$input"
  fi

  abs_bin="$(pwd)/$bin"
  out="$rundir/out.log"
  ( cd "$rundir" && OMP_NUM_THREADS=1 "$abs_bin" tiny.3d >"$out" 2>&1 )
  run_rc=$?

  status="FAIL"
  reason="no finite energy line"
  if grep -qF '*** FATAL' "$out"; then
    reason="FATAL error in output"
  elif [ $run_rc -ne 0 ]; then
    reason="nonzero exit ($run_rc)"
  elif grep -Eq '^[[:space:]]*[0-9]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+' "$out" \
     && ! grep -Eiq 'nan|infinity|\*\*\* (FATAL|ERROR)' "$out"; then
    status="PASS"
  fi

  if [ "$status" = "PASS" ]; then
    echo "PASS  $target"
  else
    echo "FAIL  $target  ($reason; see $out)"
    OVERALL=1
  fi
  RESULTS="$RESULTS
$status  $target"

  # Multi-thread sanity case (bead 0ae.4): re-run apple-silent's own binary
  # once more with OMP_NUM_THREADS=2 to catch thread-count-dependent
  # regressions (e.g. the 3zk neighbor-merge bug), using the same PASS/FAIL
  # predicate as the 1-thread run above.
  if [ "$target" = "apple-silent" ]; then
    mt_label="apple-silent (OMP_NUM_THREADS=2)"
    mt_rundir="$WORKDIR/$target-omp2"
    mkdir -p "$mt_rundir/data/test"
    mt_input="$mt_rundir/tiny.3d"
    cp "$input" "$mt_input"
    mt_out="$mt_rundir/out.log"
    ( cd "$mt_rundir" && OMP_NUM_THREADS=2 "$abs_bin" tiny.3d >"$mt_out" 2>&1 )
    mt_rc=$?

    mt_status="FAIL"
    mt_reason="no finite energy line"
    if grep -qF '*** FATAL' "$mt_out"; then
      mt_reason="FATAL error in output"
    elif [ $mt_rc -ne 0 ]; then
      mt_reason="nonzero exit ($mt_rc)"
    elif grep -Eq '^[[:space:]]*[0-9]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+[[:space:]]+[-0-9.]+' "$mt_out" \
       && ! grep -Eiq 'nan|infinity|\*\*\* (FATAL|ERROR)' "$mt_out"; then
      mt_status="PASS"
    fi

    if [ "$mt_status" = "PASS" ]; then
      echo "PASS  $mt_label"
    else
      echo "FAIL  $mt_label  ($mt_reason; see $mt_out)"
      OVERALL=1
    fi
    RESULTS="$RESULTS
$mt_status  $mt_label"
  fi
done

echo
echo "=== summary ==="
echo "$RESULTS" | grep -v '^$'

exit $OVERALL
