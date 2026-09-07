#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import shutil
import urllib.parse
import urllib.request

PMCID = "PMC13126619"
DOI = "10.1111/mec.70355"
ARTICLE_URL = f"https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/"
EXPECTED = {
    "MEC-35-e70355-s001.zip": "Figures S1-S4 and Table S1",
    "MEC-35-e70355-s002.zip": "Appendix S1 NEXUS RADseq SNP data",
    "MEC-35-e70355-s003.zip": "Appendix S2 STRUCTURE RADseq SNP data",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def request_bytes(url: str, *, timeout: int = 120) -> bytes:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "TTF independent-biogeography benchmark acquisition/0.2",
            "Accept": "text/html,application/xhtml+xml,application/zip,*/*;q=0.8",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


class LinkCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def supplement_urls() -> dict[str, str]:
    html = request_bytes(ARTICLE_URL, timeout=60).decode("utf-8", errors="replace")
    parser = LinkCollector()
    parser.feed(html)
    resolved = [urllib.parse.urljoin(ARTICLE_URL, href) for href in parser.hrefs]

    found: dict[str, str] = {}
    for expected_name in EXPECTED:
        matches = []
        for url in resolved:
            path_name = Path(urllib.parse.urlparse(url).path).name
            if path_name.lower() == expected_name.lower():
                matches.append(url)
        matches = sorted(set(matches))
        if len(matches) != 1:
            zip_names = sorted(
                {
                    Path(urllib.parse.urlparse(url).path).name
                    for url in resolved
                    if urllib.parse.urlparse(url).path.lower().endswith(".zip")
                }
            )
            raise RuntimeError(
                f"expected exactly one article link for {expected_name}; "
                f"found {matches}; article zip inventory={zip_names}"
            )
        found[expected_name] = matches[0]
    return found


def download(url: str, destination: Path) -> None:
    payload = request_bytes(url)
    if len(payload) < 1024:
        raise RuntimeError(f"download for {destination.name} was unexpectedly small: {len(payload)} bytes")
    destination.write_bytes(payload)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Acquire and freeze the open-access supplements for the Iranian Plateau TTF benchmark."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    urls = supplement_urls()
    frozen = []
    for name, description in EXPECTED.items():
        destination = output / name
        download(urls[name], destination)
        frozen.append(
            {
                "name": name,
                "description": description,
                "url": urls[name],
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        )

    manifest = {
        "schema": "ttf_iran_plateau_acquisition_v0.2",
        "source": {
            "pmcid": PMCID,
            "doi": DOI,
            "article_url": ARTICLE_URL,
            "acquisition_mode": "supplement links parsed from the public PMC article HTML",
        },
        "supplements": frozen,
        "selection": {
            "expected_files": list(EXPECTED),
            "published_roles": EXPECTED,
            "outcome_values_used": False,
            "boundary_labels_used_for_selection": False,
            "fallback_source_used": False,
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
