#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def blind(text: str) -> str:
    text = re.sub(
        r"## Data and code for peer review\n[\s\S]*?\n## Keywords",
        "## Data and code for peer review\n\n"
        "An anonymized peer-review archive containing the implementation, "
        "machine-readable method contracts, known-truth benchmark rules, "
        "frozen benchmark receipts, response-blind climate and prospective-stop "
        "receipts, tests and figure-generation code is supplied with the submission. The "
        "archive contains no version-control history or author-identifying "
        "repository metadata. The manuscript uses no genetic-response quantity "
        "from the permanently closed relational v0.3 B/C program.\n\n"
        "## Keywords",
        text,
    )
    text = re.sub(
        r"## Acknowledgements\n[\s\S]*?\n## Author contributions",
        "## Acknowledgements\n\nOmitted for double-anonymous review.\n\n"
        "## Author contributions",
        text,
    )
    text = re.sub(
        r"## Author contributions\n[\s\S]*?\n## Conflict of interest",
        "## Author contributions\n\nOmitted for double-anonymous review.\n\n"
        "## Conflict of interest",
        text,
    )
    text = re.sub(
        r"## Data availability\n[\s\S]*?\n## References",
        "## Data availability\n\n"
        "An anonymized peer-review archive containing all non-standard "
        "analysis and simulation code, frozen method contracts, scalar result "
        "handoffs and tests is supplied with the submission. Public ecological "
        "source datasets are identified in the Methods and reference list. "
        "Repository identity and permanent archival DOI will be restored after "
        "peer review/acceptance in accordance with journal policy.\n\n"
        "## References",
        text,
    )
    return text


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument(
        "--source",
        type=Path,
        default=ROOT/"manuscript/ttf_q_methods_v0.1.md",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=ROOT/"manuscript/ttf_q_methods_blinded_v0.1.md",
    )
    args=ap.parse_args()
    result=blind(args.source.read_text())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(result)
    print(args.output)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
