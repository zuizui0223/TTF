#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

SCRIPT_REF = re.compile(r"scripts/[A-Za-z0-9_./-]+\.(?:py|R|sh)")


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def referenced_scripts(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return set(SCRIPT_REF.findall(text))


def build_surface(root: Path) -> dict[str, object]:
    scripts_root = root / "scripts"
    all_scripts = {
        p.relative_to(root).as_posix()
        for p in scripts_root.rglob("*")
        if p.is_file() and p.suffix in {".py", ".R", ".sh"}
    }

    seeds: set[str] = set()
    source_files: list[str] = []

    for base in (root / ".github" / "workflows", root / "tests"):
        if not base.exists():
            continue
        for path in sorted(p for p in base.rglob("*") if p.is_file()):
            rel = path.relative_to(root).as_posix()
            source_files.append(rel)
            seeds.update(ref for ref in referenced_scripts(path) if ref in all_scripts)

    reachable = set(seeds)
    frontier = list(sorted(seeds))
    missing_references: set[str] = set()

    while frontier:
        rel = frontier.pop()
        path = root / rel
        for ref in referenced_scripts(path):
            if ref not in all_scripts:
                missing_references.add(ref)
                continue
            if ref not in reachable:
                reachable.add(ref)
                frontier.append(ref)

    historical = sorted(all_scripts - reachable)
    reachable_sorted = sorted(reachable)
    return {
        "schema": "ttf_repository_script_surface_audit_v0.1",
        "source_files_scanned": len(source_files),
        "script_count": len(all_scripts),
        "active_referenced_count": len(reachable_sorted),
        "historical_candidate_count": len(historical),
        "missing_script_references": sorted(missing_references),
        "active_referenced": [
            {"path": rel, "sha256": sha256_path(root / rel)}
            for rel in reachable_sorted
        ],
        "historical_candidates": [
            {"path": rel, "sha256": sha256_path(root / rel)}
            for rel in historical
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--output", type=Path)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    result = build_surface(args.root.resolve())
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")

    if args.check and result["missing_script_references"]:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
