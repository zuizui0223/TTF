#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import urllib.request
import xml.etree.ElementTree as ET

PMCID = "PMC13126619"
DOI = "10.1111/mec.70355"
OA_API = f"https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi?id={PMCID}"
EXPECTED_SUFFIXES = ("s001.zip", "s002.zip", "s003.zip")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "TTF benchmark acquisition/0.1"},
    )
    with urllib.request.urlopen(request, timeout=120) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)


def oa_package_url() -> str:
    request = urllib.request.Request(
        OA_API,
        headers={"User-Agent": "TTF benchmark acquisition/0.1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        xml = response.read()
    root = ET.fromstring(xml)
    links = root.findall(".//link")
    candidates = [link.attrib.get("href", "") for link in links if link.attrib.get("format") == "tgz"]
    candidates = [url for url in candidates if url]
    if len(candidates) != 1:
        raise RuntimeError(f"expected one PMC OA tgz link, found {candidates}")
    url = candidates[0]
    if url.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        url = "https://ftp.ncbi.nlm.nih.gov/" + url[len("ftp://ftp.ncbi.nlm.nih.gov/"):]
    if not url.startswith("https://"):
        raise RuntimeError(f"unexpected OA package URL scheme: {url}")
    return url


def safe_extract(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            if destination not in target.parents and target != destination:
                raise RuntimeError(f"unsafe archive path: {member.name}")
        tar.extractall(destination)


def find_supplements(root: Path) -> list[Path]:
    zips = sorted(path for path in root.rglob("*.zip") if path.is_file())
    selected: list[Path] = []
    for suffix in EXPECTED_SUFFIXES:
        matches = [path for path in zips if path.name.lower().endswith(suffix)]
        if len(matches) != 1:
            inventory = [str(path.relative_to(root)) for path in zips]
            raise RuntimeError(
                f"expected exactly one *{suffix}; found {len(matches)}; zip inventory={inventory}"
            )
        selected.append(matches[0])
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Acquire and freeze the open-access supplements for the Iranian Plateau TTF benchmark."
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    package_url = oa_package_url()
    with tempfile.TemporaryDirectory(prefix="ttf_iran_oa_") as tmp:
        tmpdir = Path(tmp)
        package = tmpdir / "oa_package.tar.gz"
        extracted = tmpdir / "extracted"
        extracted.mkdir()
        download(package_url, package)
        safe_extract(package, extracted)
        supplements = find_supplements(extracted)

        frozen = []
        for source in supplements:
            destination = output / source.name
            shutil.copy2(source, destination)
            frozen.append(
                {
                    "name": destination.name,
                    "bytes": destination.stat().st_size,
                    "sha256": sha256(destination),
                }
            )

        manifest = {
            "schema": "ttf_iran_plateau_acquisition_v0.1",
            "source": {
                "pmcid": PMCID,
                "doi": DOI,
                "oa_api": OA_API,
                "oa_package_url": package_url,
                "oa_package_sha256": sha256(package),
            },
            "supplements": frozen,
            "selection": {
                "expected_suffixes": list(EXPECTED_SUFFIXES),
                "outcome_values_used": False,
                "boundary_labels_used_for_selection": False,
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
