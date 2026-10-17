#!/usr/bin/env python3
"""Reproduce the benchmark in sequential isolated R processes (Linux)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CASES = {
    "spca-small": ("lowrank-small", ["spca", "rspca"]),
    "spca-medium": ("lowrank-medium", ["spca", "rspca"]),
    "spca-large": ("lowrank-large", ["spca", "rspca"]),
    "pca-lowrank": ("lowrank-large", ["prcomp", "rspca"]),
    "pca-flat": ("flat-small", ["prcomp", "rspca"]),
    "sparse-input": ("sparse-large", ["rspca-dense", "rspca-sparse"]),
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, default=ROOT / "R",
                        help="R source directory; can resume an existing frozen snapshot")
    parser.add_argument("--results-dir", type=Path, default=HERE / "results")
    parser.add_argument("--cases", nargs="+", choices=list(CASES), default=list(CASES))
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    args.data_dir.mkdir(parents=True, exist_ok=True)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[variable] = "1"
    env["TZ"] = "UTC"
    files = sorted(args.source_dir.glob("*.R"))
    hashes = {file.name: hashlib.sha256(file.read_bytes()).hexdigest() for file in files}
    fingerprint = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    snapshot = args.data_dir / ("source-" + fingerprint[:16])
    snapshot.mkdir(exist_ok=True)
    for file in files:
        if file.resolve() != (snapshot / file.name).resolve():
            shutil.copy2(file, snapshot / file.name)
    manifest_path = args.results_dir / "source-manifest.json"
    if manifest_path.exists():
        saved = json.loads(manifest_path.read_text())
        if saved["source_fingerprint"] != fingerprint:
            raise SystemExit("Source changed; use a new results directory to avoid mixing revisions.")
    else:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        manifest_path.write_text(json.dumps({
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
            "git_head_before_uncommitted_candidate_changes": commit,
            "source_fingerprint": fingerprint,
            "r_source_sha256": hashes,
            "note": "Working-source benchmark. Verify these implementations against the final release before transferring evidence.",
        }, indent=2) + "\n")
        with (args.results_dir / "hardware.txt").open("w") as out:
            for command in (["uname", "-a"], ["lscpu"], ["free", "-h"],
                            ["cat", "/sys/fs/cgroup/memory.max"],
                            ["cat", "/sys/fs/cgroup/cpu.max"]):
                out.write("$ " + " ".join(command) + "\n")
                out.flush()
                subprocess.run(command, stdout=out, stderr=subprocess.STDOUT, check=False)

    for input_id in dict.fromkeys(CASES[case][0] for case in args.cases):
        path = args.data_dir / f"{input_id}.rds"
        if not path.exists():
            subprocess.run([args.rscript, "--vanilla", str(HERE / "prepare-inputs.R"),
                            input_id, str(path)], env=env, check=True)

    for case in args.cases:
        input_id, methods = CASES[case]
        for repetition in range(1, args.repetitions + 1):
            # Alternate order to reduce systematic drift between methods.
            order = methods if repetition % 2 else list(reversed(methods))
            for method in order:
                name = f"{case}-{method}-{repetition}"
                record = args.results_dir / f"{name}.csv"
                timing = args.results_dir / f"{name}.time.txt"
                if record.exists():
                    continue
                command = ["/usr/bin/time", "-v", "-o", str(timing),
                           args.rscript, "--vanilla", str(HERE / "worker.R"),
                           str(snapshot), str(args.data_dir / f"{input_id}.rds"),
                           method, case, str(repetition), str(record), fingerprint,
                           str(args.results_dir / "r-environment.txt")]
                print("START", name, flush=True)
                with (args.results_dir / f"{name}.log").open("w") as log:
                    process = subprocess.run(command, env=env, stdout=log,
                                             stderr=subprocess.STDOUT, check=False)
                print((args.results_dir / f"{name}.log").read_text().strip(), flush=True)
                if process.returncode:
                    raise SystemExit(f"Benchmark failed ({process.returncode}): {name}")
                with record.open(newline="") as stream:
                    row = next(csv.DictReader(stream))
                for line in timing.read_text().splitlines():
                    if "Maximum resident set size (kbytes):" in line:
                        row["process_peak_rss_mib"] = str(float(line.split(":")[-1])/1024)
                with record.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)

    records = []
    for path in sorted(args.results_dir.glob("*.csv")):
        if path.name in ("measurements.csv", "summary.csv"):
            continue
        with path.open(newline="") as stream:
            records.extend(csv.DictReader(stream))
    if records:
        with (args.results_dir / "measurements.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    print(f"Saved {len(records)} raw measurements to {args.results_dir}", flush=True)


if __name__ == "__main__":
    main()
