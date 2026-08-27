# OPTIONS.md — Build Targets and `.3d` Input Reference

Reverse-engineered from `cluster/Makefile`, `cluster/read-3d.f`, `cluster/README`
(German, translated below) and `docs/CODE_MAP.md`. Where the German README and
the code disagreed, the code wins (noted inline). Machine-readable mirror:
`scripts/options.json` (same data, keyed the same way; validates with
`python3 -m json.tool`).

---

## (a) Make targets — MD simulation binaries

All targets below build with the default `FLAGS` (`gfortran -O3 -std=legacy
-ffixed-line-length-132`) unless noted. `apple-silent`/`apple-frenkel`
override `FLAGS` to `APPLE_FLAGS` (`-fopenmp -march=native -ffast-math
-funroll-loops -ftree-vectorize -std=legacy -ffixed-line-length-132`) for
**both** compiling and linking (GNU make prerequisite propagation), and the
`.built-with-{default,apple}-flags` stamp files force a full object rebuild
whenever a build switches between the two flag sets (mixing `-ffast-math`
and non-fast-math objects is a silent correctness hazard, not a link error).

| target | binary | species | thermostat | boundary | integrator | force | neighbor | init | stoddard |
|---|---|---|---|---|---|---|---|---|---|
| `cluster` | `cluster` | single | scalev | walls (`for-w-i9`) | int-verlet-dmm | for-a-lj-n | neighbor | init-cube1 | stoddard |
| `pbc` | `cluster` | single | scalev | pbc | int-verlet-dmm | for-a-lj-n-pbc | neighbor-pbc | init-cube1 | stoddard-pbc-n |
| `mav` | `cluster` | single | scalev | pbc | int-verlet-mav | for-a-lj-n-pbc | neighbor-pbc-mav | init-cube1 | stoddard-pbc-n |
| `mul` | `cluster-m` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav | init-cube1 | stoddard-pbc-n-m |
| `nuc` | `clu-nuc` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav (verbose) | init-cube1 | stoddard-pbc-n-m |
| `silent` | `clu-silent` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | init-cube1 | stoddard-pbc-n-m |
| `frenkel` | `clu-fre` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | init-cube1 | stoddard-pbc-n-m-fre |
| `diffusion` | `diffusion` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | **init-center** | stoddard-pbc-n-m-fre |
| `moviefrenkel` | `clu-fre-mov` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | init-cube1 | stoddard-pbc-n-m-**fre-mov** |
| `nose` | `clu-nose` | multi | **scalev-nh-m** | pbc | **int-verlet-nh-mav** | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | init-cube1 | stoddard-pbc-n-m |
| `apple-silent` | `clu-silent-omp` | multi | scalev-m | pbc | int-verlet-mav-**omp** | for-a-lj-n-pbc-m-**omp** | neighbor-pbc-mav-**omp** | init-cube1 | stoddard-pbc-n-m |
| `apple-frenkel` | `clu-fre-omp` | multi | scalev-m | pbc | int-verlet-mav-**omp** | for-a-lj-n-pbc-m-**omp** | neighbor-pbc-mav-**omp** | init-cube1 | stoddard-pbc-n-m-fre |
| `rst` | `clu-rst` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav-silent | init-cube1 | stoddard-pbc-n-m |
| `test` | `clutest` | multi | scalev-m | pbc | int-verlet-mav | for-a-lj-n-pbc-m | neighbor-pbc-mav (verbose) | init-cube1 | stoddard-pbc-n-m |
| `benchmark` | *(meta)* | — | — | — | — | — | — | — | — |

Writer/misc differences not captured in the table above (see full `objs`
list per target in `scripts/options.json`):

- **cluster** only: adds `wrt-fpt` (first-passage-time output), single-species `stopdead`.
- **pbc**, **mav**: no `wrt-fpt`; `pbc` uses `stopdead-m`, `mav` uses `stopdead` (inconsistent with every other PBC target, which uses `stopdead-m`).
- **mul**: `wrt-xyz-m` (plain multi-species XYZ, not VMD).
- **nuc**, **silent**, **frenkel**, **diffusion**, **moviefrenkel**, **nose**, **rst**, **test**, **apple-silent**, **apple-frenkel**: all use `stopdead-m`.
- **frenkel**, **diffusion**, **moviefrenkel**, **apple-frenkel**: `wrt-csi-fre` (Stillinger + Frenkel distributions) instead of `wrt-csi`.
- **silent**, **apple-silent**: `wrt-fpt` present; **frenkel**/**apple-frenkel** also have `wrt-fpt`; **diffusion**/**moviefrenkel**/**nuc**/**nose**/**rst**/**test** do not.
- **rst**: `wrt-rst-multi` (uniquely-named restart file per checkpoint) instead of `wrt-rst`; `wrt-xyz-m` not VMD.
- **test**: `wrt-vel-m-ats1` / `wrt-xyz-m-ats1` (species-1-only output).
- **nuc**, **nose**, **silent**, **apple-silent**, **frenkel**, **diffusion**, **moviefrenkel**, **apple-frenkel**: `wrt-xyz-m-vmd` (VMD-compatible, drives cluster-visualization side files).

**Discrepancies / hazards found (code wins):**
- `cluster`, `pbc`, and `mav` all link to an output binary literally named
  `cluster` (`-o cluster` in the Makefile) — building one after another
  silently overwrites the previous binary. Not documented in the README.
- `mav`'s object list uses `stopdead.o` (single-species) while `pbc` uses
  `stopdead-m.o` (multi-species-capable) even though both are otherwise
  PBC/single-species targets — looks like an inconsistency rather than an
  intentional choice; flagged here rather than silently normalized.
- `neighbor-pbc-mav-silent.o` and `neighbor-pbc-mav.o` are **the same
  source file** (`neighbor-pbc-mav.f`) compiled twice with different `-cpp`
  defines (`-DSILENT` vs none) — this is *not* a separate `.f` file, unlike
  every other `-silent`/verbose pairing in this table.
- The README documents only `neighbor.f` / `neighbor-pbc.f` /
  `neighbor-pbc-mav.f` and `int-verlet.f` / `int-verlet-dmm.f` /
  `int-verlet-mav.f` as of its last edit — it predates
  `neighbor-pbc-mav-omp.f`, `int-verlet-mav-omp.f`, `int-verlet-nh-mav.f`,
  `for-a-lj-n-pbc-m*.f`, `stoddard-pbc-n-m*.f`, `scalev-m.f`,
  `scalev-nh-m.f`, and the whole multi-species/Nose-Hoover/OpenMP feature
  set. `docs/CODE_MAP.md` is current and was treated as authoritative
  alongside the source for these.

## Helper targets (`cluster/helpers/*.f` — not MD binaries)

Post-processing / analysis programs, installed to `$(BINDIR)` (default
`$(HOME)/bin`, overridable via `make BINDIR=./bin-test <target>`). Several
share an output binary name between a box (non-PBC) and PBC variant
(`selclu`/`selclup`, `density`/`densityp`, `colgas`/`colgasp`,
`bigcl`/`bigclp` all install to the same filename — last one built wins).
Full purpose/objs list: `scripts/options.json` → `helper_targets`.

`thinyout, selats, selclu, selclup, sizer, mesher, veltest, pressw, marker,
press, density, densityp, sysdens, colgas, colgasp, rdf, timefile, makelog,
bars, reverse, rate-xyz, stotest, ratstat, ratstatt, bigcl, bigclp,
focusclp, passagetime, mpt, tempavg, meantemp, t-hist, tempavg2,
diffusivity, probability` (35 targets).

`makelog`, `reverse`, and `mpt` are undocumented in the README (marked
`??????` there) — purpose inferred from source file names / call sites only.

## Housekeeping targets

`clean` (remove `.o`/binaries/stamps), `archive` (tar+gzip the source tree),
`bdir` (create `$(BINDIR)` if missing; prerequisite of every helper target).

---

## (b) Module trait glossary

**Thermostat**
- `scalev` — single-species isokinetic velocity scaling.
- `scalev-m` — multi-species isokinetic velocity scaling, gated per species by the `apth(i)` flag from the `.3d` file. Primary production thermostat.
- `scalev-nh-m` — Nosé-Hoover coupling + velocity scaling at init; pairs only with `int-verlet-nh-mav`.
- `scalev-clu` — cluster-vs-vapor separate thermostat; **source exists but is not linked into any Makefile target** — dead code in every shipped binary (its only inputs, the `.3d` fields `eclth`/`ecltim`, are consequently unused too — see field table below).

**Boundary**
- walls — `for-w-i9.o` (or another `for-w-*.f`) supplies real repulsive forces on box faces; only the `cluster` target uses this.
- pbc — `for-w-dummy.o` (no-op) is linked instead; periodicity is controlled per-axis by `bperx`/`bpery`/`bperz` in the `.3d` file, enforced by the paired `for-a-*-pbc*.f` minimum-image force routine and `neighbor-pbc*.f` list. Every target except `cluster` uses this.

**Integrator**
- `int-verlet-dmm` — Verlet + raw cumulative displacement tracking (no smoothing) to trigger neighbor-list rebuilds.
- `int-verlet-mav` — Verlet + moving-average displacement tracking, adaptive skin radius. Primary production integrator.
- `int-verlet-mav-omp` — OpenMP-parallel `int-verlet-mav`.
- `int-verlet-nh-mav` — `int-verlet-mav` + Nosé-Hoover friction terms, iterative velocity/PSI convergence (≤1000 iters, tol 1e-10). Single active species only per `docs/CODE_MAP.md`.
- **Orphaned** (present in `cluster/`, documented in README, no Makefile target): `int-verlet.f`, `int-euler1.f`, `int-euler2.f`, `int-gear3.f`, `int-gear4.f`, `int-gear5.f`.

**Force variant** (`for-a-*.f`, subroutine `AccAtom()`)
- `for-a-lj-n` — single-species, cutoff, neighbor list, no PBC.
- `for-a-lj-n-pbc` — single-species, cutoff, neighbor list, minimum-image PBC.
- `for-a-lj-n-pbc-m` — multi-species (1-1/2-2/1-2 geometric mixing), cutoff, neighbor list, PBC. Primary production force routine.
- `for-a-lj-n-pbc-m-omp` — OpenMP-parallel `for-a-lj-n-pbc-m`.
- **Orphaned**: `for-a-lj1.f`..`for-a-lj4.f` (O(N²) historical sequence), `for-a-lj-n-pbcsw.f` (shifted potential), `for-a-dummy.f` (no-op).

**Neighbor module** (`neighbor*.f`, subroutine `neighbor(step)`)
- `neighbor` — basic list, no PBC.
- `neighbor-pbc` — PBC, fixed skin, no adaptive rebuild.
- `neighbor-pbc-mav` — PBC, adaptive skin via moving-average displacement; per-step diagnostic printed.
- `neighbor-pbc-mav-silent` — **same source** as `neighbor-pbc-mav.f`, compiled with `gfortran -cpp -DSILENT` into a differently-named object; suppresses the diagnostic. Used by every "production" multi-species target.
- `neighbor-pbc-mav-omp` — OpenMP-parallel adaptive PBC list.
- **Orphaned**: `neighbor-dummy.f`, `neighbor-pbcsw-mav.f`.

**Cluster search (`stoddard*.f`, subroutine `stoddard(step)`)**
- `stoddard` / `stoddard-pbc` — O(N²) Stillinger search, no-PBC / PBC. (`stoddard-pbc` is used only by helper programs, never an MD binary target.)
- `stoddard-pbc-n` — Stillinger, PBC, neighbor-list-optimized, single-species.
- `stoddard-pbc-n-m` — same, multi-species.
- `stoddard-pbc-n-m-fre` — adds the ten Wolde/Frenkel liquid-monomer criterion.
- `stoddard-pbc-n-m-fre-mov` — adds cluster-movement tracking for visualization.
- **Orphaned**: `stoddard-pbcsw-n.f`, `stoddard-dummy.f`.

**Writers** — see the per-target table above for which combination each target uses; `wrt-csi` vs `wrt-csi-fre` (Frenkel distribution too), `wrt-rst` vs `wrt-rst-multi` (numbered files), and the `wrt-xyz*` family (`wrt-xyz`, `wrt-xyz-m`, `wrt-xyz-m-vmd`, `wrt-xyz-m-ats1`/`wrt-vel-m-ats1`, plus `wrt-xyz1`/`2`/`3`/`wrt-xyz-m1` which are helper-program-only, never an MD binary target).

---

## (c) `.3d` input field reference

All fields are read with `implicit none`, fixed-column `read (1,'(...)')`
statements in `subroutine readdata()` (`cluster/read-3d.f`). Most string and
float fields use column offset `41X` (i.e. the value starts at column 42, to
the right of a `: `-terminated 41-character label); the file version gate
(`fver`, line 58) is the sole exception at `41X,F8.6`, and **`ntim`
(read-3d.f:114) uses `40X,I10` — one column earlier than every other
`I10` field** (see discrepancy note in the table). `sample.3d` line numbers
below are for the shipped `fver=1.30` sample; version-gated fields absent
from earlier-version files fall back to hardcoded defaults set in the
`else` branches (see `read-3d.f` for exact fallback values).

Full machine-readable table (all 56 top-level fields + 10 per-species
fields + 12 legacy pre-1.10 fields = **78 fields total**, one entry per
`READ` statement in `read-3d.f`) is in `scripts/options.json` →
`input_3d_fields` / `input_3d_per_species_fields` /
`input_3d_legacy_pre_1_10_fields`. Summary:

### Header / output file table (read-3d.f:57–110, `.3d` lines 1–18)
| field | source | `.3d` line | meaning |
|---|---|---|---|
| `fver` | read-3d.f:58 | 2 | format version; must be ≥1.02 |
| `fpath` | read-3d.f:66 | 4 | directory prepended to every filename below |
| `fxyz`/`ftxyz` | read-3d.f:68,71 | 5/6 | trajectory filename / output interval (0=off) |
| `fvel`/`ftvel` | read-3d.f:73,76 | 7/8 | velocity filename / interval |
| `facc`/`ftacc` | read-3d.f:77,80 | 9/10 | acceleration filename / interval |
| `frst`/`ftrst` | read-3d.f:81,84 | 11/12 | restart filename / interval |
| `fcsi`/`ftcsi` | read-3d.f:85,95 | 13/14 | cluster-size filename / interval |
| `fave`/`ftave` | read-3d.f:96,99 | 15/16 | averages filename / interval |
| `frat`/`ftrat` | read-3d.f:101,104 | 17/18 | rate filename / interval (fver≥1.06) |

`ffpt`/`fxfpt` (first-passage-time output) and `fcsif`/`fxcsif` (Frenkel
`.fsi` output) are **not read from the file** — they're derived at
read-3d.f:88-93 from hardcoded literals `'fpt'`/`'fsi'` plus `fpath`.

### Time / box / thermostat-control table (read-3d.f:112–207, `.3d` lines 19–63)
| field | source | `.3d` line | meaning |
|---|---|---|---|
| `dt` | read-3d.f:113 | 20 | timestep [fs] |
| `ntim` | read-3d.f:114 | 21 | total timesteps (⚠ `40X` column offset, not `41X`) |
| `btim` | read-3d.f:115 | 22 | starting timestep number |
| `sctim` | read-3d.f:116 | 23 | screen-output interval |
| `boxx/y/z` | read-3d.f:118-120 | 25-27 | box size [nm]→[Å] (×10) |
| `bperx/y/z` | read-3d.f:122-124 | 28-30 | PBC per axis, X/- (fver≥1.08) |
| `bpot`, `beng` | read-3d.f:130-131 | 31-32 | wall potential exponent / energy [K]→internal |
| `bcor` | read-3d.f:132 | 33 | show box corners in XYZ, X/- |
| `bshr`, `boxxn/yn/zn`, `boxxs/ys/zs`, `boxt` | read-3d.f:134-141 | 34-41 | box resize target/step/interval (fver≥1.05) |
| `atsorts` (legacy) | read-3d.f:163 | *(none — pre-1.10 only)* | species count, superseded by line 63 field for fver≥1.10 |
| `frctrl`, `finp` | read-3d.f:164-165 | 43-44 | use restart file? / restart filename |
| `ethst` | read-3d.f:168 | 45 | thermostat type: X/N/- |
| `Q` | read-3d.f:171 | 46 | Nosé-Hoover coupling constant (fver≥1.2) |
| `rndctrl` | read-3d.f:174 | 47 | random velocities after init, X/- |
| `eclth`, `ecltim` | read-3d.f:176-177 | 48-49 | cluster thermostat flag/timer (fver≥1.07) — **unused**: only consumer `scalev-clu.f` is unlinked |
| `att`, `atskin` | read-3d.f:184-185 | 51-52 | system temp [K] / Stillinger criterion [σ] (fver≥1.10) |
| `atnenum` | read-3d.f:187 | 53 | Frenkel neighbor-count threshold (fver≥1.30) |
| `atrcut`, `atrskin` | read-3d.f:189-190 | 54-55 | LJ cutoff / neighbor skin [σ] |
| `nbtim` | read-3d.f:191 | 56 | neighbor-list update interval |
| `chgNtim` | read-3d.f:194 | 57 | atom-count-change interval (fver≥1.11) |
| `clnuctrl`, `csover`, `clwait` | read-3d.f:200-202 | 59-61 | nucleation abort flag / min size / persistence steps (fver≥1.12) |
| `atsorts` | read-3d.f:205 | 63 | number of atom species (fver≥1.10 position — the one that matters for the shipped sample) |

### Per-species block (read-3d.f:229-254, `.3d` lines 65–75 per species, repeats `atsorts` times)
| field | source | meaning |
|---|---|---|
| `apsym(i)` | read-3d.f:230 | symbol |
| `apnum(i)` | read-3d.f:231 | atom count |
| `apmas(i)` | read-3d.f:232 | mass [u]→internal |
| `apsig(i)` | read-3d.f:233 | LJ σ [nm]→[Å] |
| `apeps(i)` | read-3d.f:234 | LJ ε [K]→internal |
| `appr(i)`, `appa(i)` | read-3d.f:235-236 | repulsive/attractive LJ exponents |
| `apth(i)` | read-3d.f:237 | thermostat flag, X/- |
| `apNch(i)`, `apnnum(i)` | read-3d.f:239-240 | atom-count-change flag / target count (fver≥1.11) |

### Legacy pre-1.10 single-species block (read-3d.f:208-225, no lines in `fver≥1.10` files)
`ats, natom, atm, atr, ate, atpr, atpa, att, atskin, atrcut, atrskin, nbtim`
— only executed for `fver<1.10`; not exercised by the shipped `sample.3d`
(`fver=1.30`). Field meanings mirror the per-species/system-parameter
fields above (single implicit species). See `scripts/options.json` →
`input_3d_legacy_pre_1_10_fields`.

---

## Sources

- `cluster/Makefile` (all targets/objs)
- `cluster/read-3d.f` (every field, cited by line)
- `cluster/sample.3d` (line-number ground truth for the current format)
- `cluster/README` (German; module purpose translations)
- `docs/CODE_MAP.md` (current, post-README module documentation)
