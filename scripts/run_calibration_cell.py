#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ttf.calibration import run_calibration


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one TTF calibration cell.")
    parser.add_argument("--shared-fraction", type=float, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--permutations", type=int, default=199)
    parser.add_argument("--species", type=int, default=40)
    parser.add_argument("--records-per-species", type=int, default=60)
    parser.add_argument("--bandwidth", type=float, default=0.2)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cells = run_calibration(
        shared_fractions=[args.shared_fraction],
        amplitudes=[args.amplitude],
        n_replicates=args.replicates,
        n_permutations=args.permutations,
        n_species=args.species,
        records_per_species=args.records_per_species,
        bandwidth=args.bandwidth,
        alpha=args.alpha,
        seed=args.seed,
    )
    if len(cells) != 1:
        raise RuntimeError("single-cell runner produced an unexpected number of cells")

    payload = {
        "schema": "ttf_calibration_cell_v0.1",
        "config": {
            "shared_fraction": args.shared_fraction,
            "amplitude": args.amplitude,
            "replicates": args.replicates,
            "permutations": args.permutations,
            "species": args.species,
            "records_per_species": args.records_per_species,
            "bandwidth": args.bandwidth,
            "alpha": args.alpha,
            "seed": args.seed,
        },
        "cell": cells[0].to_dict(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["cell"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
