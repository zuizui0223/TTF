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

    refs_start = text.find("## References")
    data_start = text.find("## Data and Code Availability")
    figures_start = text.find("## Figure legends")
    provenance_start = text.find("## Repository provenance")
    if min(refs_start, data_start, figures_start) < 0:
        raise ValueError(
            "could not locate References / Data and Code Availability / Figure legends"
        )

    end = len(text) if provenance_start < 0 else provenance_start
    body_end = min(refs_start, data_start, figures_start)
    refs_block = text[refs_start:data_start].strip()
    figures_block = text[figures_start:end].strip()

    text = (
        text[:body_end].rstrip()
        + "\n\n"
        + refs_block
        + "\n\n"
        + ANON_DATA_CODE.strip()
        + "\n\n"
        + figures_block
        + "\n"
    )

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
