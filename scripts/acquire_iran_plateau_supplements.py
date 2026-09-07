#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import urllib.request
import zipfile

PMCID = "PMC13126619"
DOI = "10.1111/mec.70355"
SUPPLEMENT_ARCHIVE_URL = (
    f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
)
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


def download(url: str, destination: Path, *, timeout: int = 180) -> None:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "TTF independent-biogeography benchmark acquisition/0.3",
            "Accept": "application/zip,application/octet-stream,*/*;q=0.8",
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)
    if destination.stat().st_size < 1024:
        raise RuntimeError(
            f"supplement archive was unexpectedly small: {destination.stat().st_size} bytes"
        )


def select_members(archive: Path) -> dict[str, str]:
    with zipfile.ZipFile(archive) as zf:
        inventory = [name for name in zf.namelist() if not name.endswith("/")]
    selected: dict[str, str] = {}
    for expected in EXPECTED:
        matches = [name for name in inventory if Path(name).name.lower() == expected.lower()]
        if len(matches) != 1:
            raise RuntimeError(
                f"expected exactly one {expected} in Europe PMC supplementary archive; "
                f"found {matches}; inventory={inventory}"
            )
        selected[expected] = matches[0]
    return selected


def safe_extract_member(archive: Path, member: str, destination: Path) -> None:
    with zipfile.ZipFile(archive) as zf:
        info = zf.getinfo(member)
        if info.is_dir():
            raise RuntimeError(f"expected file member, got directory: {member}")
        with zf.open(info) as src, destination.open("wb") as out:
            shutil.copyfileobj(src, out)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Acquire and freeze open-access supplements for the Iranian Plateau TTF benchmark."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="ttf_iran_epmc_") as tmp:
        outer = Path(tmp) / "europepmc-supplementary-files.zip"
        download(SUPPLEMENT_ARCHIVE_URL, outer)
        selected = select_members(outer)

        frozen = []
        for name, description in EXPECTED.items():
            destination = output / name
            safe_extract_member(outer, selected[name], destination)
            if not zipfile.is_zipfile(destination):
                raise RuntimeError(f"published supplement is not a valid ZIP: {name}")
            frozen.append(
                {
                    "name": name,
                    "description": description,
                    "archive_member": selected[name],
                    "bytes": destination.stat().st_size,
                    "sha256": sha256(destination),
                }
            )

        manifest = {
            "schema": "ttf_iran_plateau_acquisition_v0.3",
            "source": {
                "pmcid": PMCID,
                "doi": DOI,
                "supplement_archive_url": SUPPLEMENT_ARCHIVE_URL,
                "supplement_archive_sha256": sha256(outer),
                "acquisition_mode": "Europe PMC REST supplementaryFiles endpoint",
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
