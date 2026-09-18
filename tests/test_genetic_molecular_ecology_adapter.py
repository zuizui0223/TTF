from __future__ import annotations

import re
from pathlib import Path


ADAPTER = Path("docs/GENETIC_TTF_MOLECULAR_ECOLOGY_ADAPTER.md")


def test_molecular_ecology_adapter_is_format_only_and_under_abstract_limit() -> None:
    text = ADAPTER.read_text()

    match = re.search(
        r"## Candidate abstract — \d+ words\n\n(.+?)\n\n## Candidate keywords",
        text,
        flags=re.S,
    )
    assert match is not None
    abstract = match.group(1).strip()
    words = abstract.split()
    assert len(words) == 213
    assert len(words) <= 250

    assert "211 survived" in abstract
    assert "103-training/108-evaluation" in abstract
    assert "T = 0.0326" in abstract
    assert "profiled-private p = 0.7393" in abstract
    assert "S = -0.2261" in abstract
    assert "upper-tail p = 0.0090" in abstract
    assert "consistent with lineage-conditioned spatial structure" in abstract

    keywords = re.findall(r"^- (.+)$", text.split("## Candidate keywords", 1)[1].split("## Cover-letter", 1)[0], flags=re.M)
    assert 4 <= len(keywords) <= 6

    assert "candidate first-shot adapter" in text
    assert "journal choice remains a human editorial decision" in text
    assert "formatting only" in text
    assert "proof of zero transfer" in text
    assert "specific historical mechanism" in text
