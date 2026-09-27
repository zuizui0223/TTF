#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path


RUNNING_TITLE = "Dimensions of butterfly specialization"

ANON_DATA_CODE = """## Data and Code Availability

For double-anonymous review, the analysis code, frozen scientific protocols, provenance receipts and figure-generation workflow will be supplied through an anonymized review repository. The large external datasets are obtained from their original providers (LepTraits, HOSTS, WCVP, GBIF, CHELSA and WGSRPD) using the versions, identifiers and query rules described in the Methods.

A permanent public archival snapshot and DOI will replace this anonymized review statement in the final public version.

"""


def render_blinded(text: str) -> str:
    text = re.sub(
        r"^\*\*Working manuscript[^\n]*\*\*\s*\n+",
        "",
        text,
        count=1,
        flags=re.MULTILINE,
    )

    lines = text.splitlines()
    if not lines or not lines[0].startswith("# "):
        raise ValueError("expected manuscript title as first Markdown heading")
    if len(lines) < 2 or not lines[1].startswith("**Running title:**"):
        lines.insert(1, f"**Running title:** {RUNNING_TITLE}")
        lines.insert(2, "")
    text = "\n".join(lines).rstrip() + "\n"

    text = text.replace(
        "Genetic responses from the original TTF program were not used in these analyses.",
        "Genetic responses from the precursor transferability program were not used in these analyses.",
    )

    data_start = text.find("## Data and Code Availability")
    if data_start < 0:
        raise ValueError("could not locate Data and Code Availability section")
    next_section = text.find("\n## ", data_start + len("## Data and Code Availability"))
    if next_section < 0:
        data_end = len(text)
    else:
        data_end = next_section + 1
    text = text[:data_start] + ANON_DATA_CODE + text[data_end:]

    provenance_start = text.find("## Repository provenance")
    if provenance_start >= 0:
        text = text[:provenance_start].rstrip() + "\n"

    forbidden = (
        "TTF repository",
        ".github/workflows/",
        "Repository provenance",
        "zuizui0223",
    )
    for token in forbidden:
        if token.lower() in text.lower():
            raise ValueError(f"blinded manuscript still contains identifying/internal token: {token}")

    return text


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    rendered = render_blinded(args.input.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
