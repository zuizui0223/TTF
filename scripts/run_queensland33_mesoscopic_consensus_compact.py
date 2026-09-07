#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


SIGNALS = (
    ("0", "0.5", "s0-a0.5"),
    ("0", "1", "s0-a1"),
    ("0", "2", "s0-a2"),
    ("0", "3", "s0-a3"),
    ("0", "4", "s0-a4"),
    ("0.25", "2", "s0.25-a2"),
    ("0.5", "2", "s0.5-a2"),
    ("1", "2", "s1-a2"),
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the frozen Queensland v0.3b consensus Gate in one process slot.")
    parser.add_argument("--layouts", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--aggregate-output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--bootstrap", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260911)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for layout in range(12):
        for shared, amplitude, slug in SIGNALS:
            output = args.output_dir / f"layout-{layout}-{slug}.json"
            cmd = [
                sys.executable,
                "scripts/run_queensland33_mesoscopic_consensus_calibration.py",
                "--layouts", str(args.layouts),
                "--layout-index", str(layout),
                "--shared-fraction", shared,
                "--amplitude", amplitude,
                "--replicates", str(args.replicates),
                "--bootstrap", str(args.bootstrap),
                "--alpha", "0.05",
                "--k", "2",
                "--noise-sd", "0.8",
                "--transition-width", "0.20",
                "--min-shared-crossing-species", "14",
                "--prior-strength", "1.0",
                "--evidence-temperature", "1.0",
                "--split-seed", "20260907",
                "--seed", str(args.seed),
                "--output", str(output),
            ]
            subprocess.run(cmd, check=True)

    subprocess.run(
        [
            sys.executable,
            "scripts/aggregate_queensland33_mesoscopic_consensus.py",
            "--input-dir", str(args.output_dir),
            "--type1-ceiling", "0.10",
            "--power-floor", "0.80",
            "--output", str(args.aggregate_output),
        ],
        check=True,
    )
    print(args.aggregate_output.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
