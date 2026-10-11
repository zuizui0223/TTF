#!/usr/bin/env python3
"""Verify and safely extract the exact original phylogatR ZIP for MASK-ONLY processing.

This is an input-transport utility, not permission to open nucleotide identity.
No FASTA sequence bytes are interpreted here. The authenticated original is
checked by size and SHA-256 before any member is extracted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import stat
import zipfile
from pathlib import Path, PurePosixPath

ARCHIVE_BYTES = 274_988_692
ARCHIVE_SHA = "5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_member(info: zipfile.ZipInfo):
    name = info.filename
    if "\\" in name or not name or name.startswith("/"):
        raise RuntimeError("unsafe member name in exact source archive")
    parts = PurePosixPath(name).parts
    if ".." in parts:
        raise RuntimeError("path traversal in exact source archive")
    mode = info.external_attr >> 16
    if stat.S_ISLNK(mode):
        raise RuntimeError("symlink prohibited in exact source archive")
    if info.file_size < 0 or info.file_size > 2_000_000_000:
        raise RuntimeError("invalid file size in exact source archive")


def extract_verified_archive(
    archive: Path, output_dir: Path, *, expected_sha: str = ARCHIVE_SHA,
    expected_size: int = ARCHIVE_BYTES,
) -> Path:
    if archive.stat().st_size != expected_size or sha256(archive) != expected_sha:
        raise RuntimeError("exact frozen phylogatR archive identity mismatch")
    output_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        entries = z.infolist()
        if not entries or len(entries) > 300_000:
            raise RuntimeError("unexpected exact source archive member count")
        if sum(x.file_size for x in entries) > 8_000_000_000:
            raise RuntimeError("source extraction exceeds frozen bounded safety allowance")
        for x in entries:
            check_member(x)
        for x in entries:
            z.extract(x, output_dir)
    hits = sorted(p.parent for p in output_dir.rglob("genes.txt")
                  if (p.parent / "cite.txt").is_file())
    if len(hits) != 1:
        raise RuntimeError("exact archive must expose one phylogatR root with genes.txt and cite.txt")
    root = hits[0].resolve()
    if root != output_dir.resolve() and output_dir.resolve() not in root.parents:
        raise RuntimeError("raw source escaped verified extraction root")
    return root


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-archive", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--root-receipt", type=Path, required=True)
    args = ap.parse_args()
    root = extract_verified_archive(args.source_archive, args.output_dir)
    # Only a filesystem-relative directory is persisted, never sequence content.
    relative = root.relative_to(args.output_dir.resolve())
    args.root_receipt.parent.mkdir(parents=True, exist_ok=True)
    args.root_receipt.write_text(json.dumps({
        "schema": "ttf_historical_host_memory_exact_mask_source_intake_v0.1",
        "status": "PASS_EXACT_ARCHIVE_ROOT_FOR_MASK_ONLY",
        "archive_sha256": ARCHIVE_SHA,
        "archive_size_bytes": ARCHIVE_BYTES,
        "relative_root": str(relative),
        "sequence_identity_opened": False,
        "genetic_distance_opened": False,
    }, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS_EXACT_ARCHIVE_ROOT_FOR_MASK_ONLY",
                      "relative_root":str(relative)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
