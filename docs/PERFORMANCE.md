# Performance

## Multi-thread validated: yes (bit-exact below 2500 atoms; ~1e-15 relative agreement above)

**Configuration validated:** `apple-silent` target (binary `clu-silent-omp`),
linking `neighbor-pbc-mav-omp.o` + `for-a-lj-n-pbc-m-omp.o` +
`int-verlet-mav-omp.o` (see `cluster/Makefile` `APPLE_SILENT` object list),
built with `APPLE_FLAGS` (`-fopenmp -march=native -ffast-math`). Validated at
`OMP_NUM_THREADS=1` vs `4` on a 12-core Apple Silicon host (6 performance
cores), `natom=1000` (sample.3d, below the `NATOM_PAR_THRESHOLD=2500` force-
loop guard in `for-a-lj-n-pbc-m-omp.f`), 200 timesteps.

**Determinism arm** (fixed `CLUSTER_SEED`, 1 vs 4 threads, same binary/input):
screen energy output and the binary `test.rst` restart file were **byte-
identical** between thread counts (not just approximately equal).

**Ensemble arm** (16 replicas, distinct seeds, same seed set reused across
both thread counts, `scripts/campaign.py --threads`): every one of the 16
seed-matched replica pairs produced byte-identical `out.log` between
1-thread and 4-thread runs; mean/std of final epot and T were therefore
identical to all reported digits (T: mean 150.0013392873, std 0.0345130300,
n=16, both arms; epot: mean -0.0002570991, std 0.0001373603, n=16, both
arms).

**Why agreement is exact rather than merely approximate:** at `natom=1000`
(below the 2500-atom threshold), `AccAtom`'s parallel branch is skipped
entirely in favor of the serial, term-identical force loop (see
`for-a-lj-n-pbc-m-omp.f` lines 106-152), so pairwise force accumulation order
is independent of `OMP_NUM_THREADS`. The neighbor list build
(`neighbor-pbc-mav-omp.f`) writes per-atom-disjoint scratch slots
(`local_counts(i)`, `local_neighbors(:,i)`) regardless of how OpenMP's
`SCHEDULE(DYNAMIC,32)` interleaves iterations across threads, then merges
them into `atnlist`/`atnidx` in a single, fixed `i=1..natom` order inside an
`!$OMP SINGLE` block — so list order and content are thread-count-invariant
too. The one thread-count-sensitive floating-point reduction in this
configuration is `ekin` (`REDUCTION(+:ekin)` in `int-verlet-mav-omp.f`); at
these thread counts and this problem size it did not perturb any observed
digit, but reduction-order-driven float divergence remains a latent
possibility at larger `natom` (where the `for-a-lj-n-pbc-m-omp.f` atomic
branch activates) or on different core counts/compilers — treat "exact"
agreement as this configuration's observed behavior, not a guarantee.

**Caveats:**
- RNG (`kpsRND`/`kpsINT` in `kps.f`, `COMMON /kps/` state, no locking) is
  only called from `kpsInit` (pre-run) and `chgNatom` (growth path, called
  from the serial step loop in `main.f`, never from inside an `!$OMP
  PARALLEL` region in the objects this target links). It is not called
  concurrently by multiple threads in this configuration. If a future target
  calls RNG from inside a parallel region, determinism would need to be
  re-validated.
- Above `NATOM_PAR_THRESHOLD=2500` atoms, `for-a-lj-n-pbc-m-omp.f` switches
  to an atomic-update parallel branch. Measured (review run, natom=2700 —
  500 Ar + 2200 He derived from sample.3d, 100 steps, fixed seed, 1 vs 4
  threads): screen energies identical to all 10 printed digits; binary
  `test.rst` NOT byte-identical — max abs diff 2.3e-18, max relative 3.0e-15
  across 2589 differing doubles. 1-thread runs are bit-reproducible
  run-to-run; 4-thread runs are NOT bit-reproducible even against themselves
  (rerun max relative 2.3e-16) — expected consequence of `!$OMP ATOMIC`
  accumulation-order interleaving. Above the threshold, expect ~1e-15
  relative agreement, not bit-exactness, and no run-to-run bit
  reproducibility at >1 thread.
- Verified statically (not exercised at runtime, since sample.3d has no
  zero-neighbor atoms): an atom with zero neighbors cannot cause a read of
  uninitialized `istart`/`iend` in `for-a-lj-n-pbc-m-omp.f`. `istart(i) =
  atnidx(i)+1`, `iend(i) = atnidx(i+1)-1` (lines 82-85); for a zero-neighbor
  atom `atnidx(i+1) = atnidx(i)+1`, giving `istart(i) = iend(i)+1`, i.e. an
  empty `do k=istart(i),iend(i)` range, which Fortran executes zero times.
  `atnidx` itself is fully populated for `i=1..natom` inside the `!$OMP
  SINGLE` block of `neighbor-pbc-mav-omp.f` (lines 123-147) before the
  implicit barrier at `!$OMP END SINGLE`, so by the time `AccAtom`'s
  parallel loop (`for-a-lj-n-pbc-m-omp.f` lines 170-219) reads `istart`/
  `iend`, every entry is initialized and read-only.

## Multi-thread smoke sanity case

`cluster/smoke.sh` includes an `apple-silent` run with `OMP_NUM_THREADS=2`
in addition to the existing single-thread pass, using the same FATAL/
Infinity/exit-code PASS criteria as the rest of the suite.
