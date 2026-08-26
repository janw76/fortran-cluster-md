#!/usr/bin/env python3
"""campaign.py -- run N replicas of a cluster-md binary in parallel, each with
its own CLUSTER_SEED and its own isolated output directory.

Output isolation: the .3d input's "path for following files (relative)" line
names a directory (e.g. data/test/) that all cluster output files are written
into, *relative to the process's cwd* (see cluster/read-3d.f / sample.3d).
So each replica is run with cwd=rep-NNN/ and that relative dir pre-created
there (cluster/smoke.sh does the same: `mkdir -p data/test` per run).

Seed spacing: kpsInit (cluster/kps.f) coerces CLUSTER_SEED via
ibm=ior(abs(seed),1) -- only the low bit is forced, so two seeds differing by
1 can collide (abs(2)|1 == abs(3)|1 == 3). Stepping by 2 per replica keeps all
seed_base+2*i values the same parity, so their |1 coercions stay distinct.
"""
import argparse
import concurrent.futures
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone


def perf_core_count():
    try:
        out = subprocess.run(
            ["sysctl", "-n", "hw.perflevel0.logicalcpu"],
            capture_output=True, text=True, timeout=5,
        )
        if out.returncode == 0 and out.stdout.strip():
            return int(out.stdout.strip())
    except (OSError, ValueError):
        pass
    return os.cpu_count() or 1


def output_subdir(input_path):
    """Parse the "path for following files (relative)" line of a .3d input."""
    with open(input_path) as f:
        for line in f:
            if line.startswith("path for following files"):
                return line.split(":", 1)[1].strip()
    return None


def run_replica(index, binary, input_path, campaign_dir, seed, subdir, threads):
    repdir = os.path.join(campaign_dir, f"rep-{index:03d}")
    os.makedirs(repdir, exist_ok=True)
    if subdir:
        os.makedirs(os.path.join(repdir, subdir), exist_ok=True)

    input_name = os.path.basename(input_path)
    shutil.copy(input_path, os.path.join(repdir, input_name))

    env = dict(os.environ)
    env["CLUSTER_SEED"] = str(seed)
    env["OMP_NUM_THREADS"] = str(threads)

    log_path = os.path.join(repdir, "out.log")
    start = time.time()
    with open(log_path, "w") as log:
        proc = subprocess.run(
            [binary, input_name], cwd=repdir, env=env,
            stdout=log, stderr=subprocess.STDOUT,
        )
    duration = time.time() - start

    log_text = open(log_path).read()
    fatal_detected = "*** FATAL" in log_text
    failed = proc.returncode != 0 or fatal_detected

    return {
        "index": index,
        "dir": repdir,
        "seed_env": seed,
        "expected_seed": abs(seed) | 1,
        "exit_code": proc.returncode,
        "duration_s": duration,
        "failed": failed,
        "fatal_detected": fatal_detected,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Run N replicas of a cluster-md binary concurrently, "
                     "each with a distinct CLUSTER_SEED and isolated cwd.")
    ap.add_argument("binary", help="path to the cluster-md binary")
    ap.add_argument("input", help="path to the .3d input file")
    ap.add_argument("-n", "--replicas", type=int, required=True,
                     help="number of replicas to run")
    ap.add_argument("--jobs", type=int, default=perf_core_count(),
                     help="max concurrent replicas "
                          "(default: performance-core count)")
    ap.add_argument("--seed-base", type=int, default=1,
                     help="base seed; replica i uses seed_base + 2*i "
                          "(default: 1)")
    ap.add_argument("--threads", type=int, default=1,
                     help="OMP_NUM_THREADS per replica (default: 1)")
    args = ap.parse_args()

    binary = os.path.abspath(args.binary)
    input_path = os.path.abspath(args.input)
    if not os.path.isfile(binary):
        sys.exit(f"error: binary not found: {binary}")
    if not os.path.isfile(input_path):
        sys.exit(f"error: input not found: {input_path}")

    subdir = output_subdir(input_path)
    timestamp = datetime.now().strftime("%Y%m%dT%H%M%S")
    campaign_dir = os.path.join("runs", timestamp)
    os.makedirs(campaign_dir, exist_ok=True)

    start_ts = datetime.now(timezone.utc).isoformat()
    t0 = time.time()

    replicas = [None] * args.replicas
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(
                run_replica, i, binary, input_path, campaign_dir,
                args.seed_base + 2 * i, subdir, args.threads,
            ): i
            for i in range(args.replicas)
        }
        for fut in concurrent.futures.as_completed(futures):
            result = fut.result()
            replicas[result["index"]] = result
            status = "FAIL" if result["failed"] else "PASS"
            print(f"{status}  rep-{result['index']:03d}  "
                  f"seed={result['expected_seed']}  "
                  f"exit={result['exit_code']}  "
                  f"{result['duration_s']:.2f}s")

    wall_s = time.time() - t0
    end_ts = datetime.now(timezone.utc).isoformat()
    serial_sum_s = sum(r["duration_s"] for r in replicas)

    manifest = {
        "binary": binary,
        "input": input_path,
        "jobs": args.jobs,
        "seed_base": args.seed_base,
        "threads": args.threads,
        "start": start_ts,
        "end": end_ts,
        "wall_s": wall_s,
        "serial_sum_s": serial_sum_s,
        "replicas": replicas,
    }
    with open(os.path.join(campaign_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    n_failed = sum(1 for r in replicas if r["failed"])
    print(f"\n{args.replicas - n_failed}/{args.replicas} replicas passed  "
          f"wall={wall_s:.2f}s  serial_sum={serial_sum_s:.2f}s")
    print(f"manifest: {os.path.join(campaign_dir, 'manifest.json')}")

    sys.exit(1 if n_failed else 0)


if __name__ == "__main__":
    main()
