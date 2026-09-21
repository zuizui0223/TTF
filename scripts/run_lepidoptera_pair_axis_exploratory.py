#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import stat
import tempfile
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from ttf.lepidoptera_species_transferability import deterministic_pair_folds
from ttf.lepidoptera_trait_gradient_empirical import (
    frozen_pair_surface,
    pair_transfer_congruence,
    post_ibd_edge_response,
)
from ttf.phylogatr_confirmatory import (
    choose_one_panel_per_species,
    read_genes_rows,
    scan_phylogatr_phase1,
)
from ttf.phylogatr_empirical import extract_species_frozen_edge_distances

ALIASES = (
    "COI", "CO1", "COX1", "COXI",
    "CYTOCHROME C OXIDASE SUBUNIT I",
    "CYTOCHROME C OXIDASE SUBUNIT 1",
)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def selective_extract(archive: Path, destination: Path, survivors: set[str]) -> Path:
    with zipfile.ZipFile(archive) as z:
        infos = {info.filename: info for info in z.infolist()}
        for info in infos.values():
            p = Path(info.filename)
            if p.is_absolute() or ".." in p.parts:
                raise RuntimeError("unsafe archive member")
            mode = (int(info.external_attr) >> 16) & 0o170000
            if mode == stat.S_IFLNK:
                raise RuntimeError("symlink archive member forbidden")
        genes_members = [
            name for name in infos
            if name.endswith("genes.txt")
            and str(Path(name).parent / "cite.txt").replace("\\", "/") in infos
        ]
        if len(genes_members) != 1:
            raise RuntimeError("archive root drift")
        genes_member = genes_members[0]
        prefix = str(Path(genes_member).parent).replace("\\", "/")
        reader = csv.DictReader(
            io.StringIO(z.read(genes_member).decode("utf-8")),
            delimiter="\t",
        )
        selected = {
            genes_member,
            str(Path(prefix) / "cite.txt").replace("\\", "/"),
        }
        for row in reader:
            species = str(row.get("species", "")).strip()
            if species not in survivors:
                continue
            rel = Path(str(row.get("dir", "")))
            if rel.is_absolute() or ".." in rel.parts:
                raise RuntimeError("unsafe genes.txt dir")
            gene = str(row.get("gene", ""))
            selected.add(str(Path(prefix) / rel / f"{gene}.afa").replace("\\", "/"))
            selected.add(str(Path(prefix) / rel / "occurrences.txt").replace("\\", "/"))
        missing = sorted(name for name in selected if name not in infos)
        if missing:
            raise RuntimeError(f"survivor source files missing: {missing[:5]}")
        for name in sorted(selected):
            z.extract(infos[name], destination)
    root = destination / Path(genes_member).parent
    if not (root / "genes.txt").is_file() or not (root / "cite.txt").is_file():
        raise RuntimeError("selective extraction drift")
    return root


def response_item(item):
    name, panel = item
    distances = extract_species_frozen_edge_distances(
        panel.fasta_path,
        panel.occurrence_path,
        panel.canonical_latlon,
        panel.geometry,
        minimum_comparable_fraction=0.5,
    )
    return name, post_ibd_edge_response(
        distances.genetic_distance,
        panel.geometry,
        min_training_edges=5,
    )


def load_traits(path: Path) -> dict[str, dict[str, str]]:
    out = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            name = str(row["Species"]).strip()
            if name and name not in out:
                out[name] = row
    return out


def present(value: object) -> bool:
    return str(value).strip() not in {"", "NA"}


def wing_value(row: dict[str, str]) -> float:
    keys = (
        "WS_L", "WS_U", "FW_L", "FW_U",
        "WS_L_Fem", "WS_U_Fem", "WS_L_Mal", "WS_U_Mal",
        "FW_L_Fem", "FW_U_Fem", "FW_L_Mal", "FW_U_Mal",
    )
    values = [float(row[key]) for key in keys if present(row.get(key, ""))]
    if not values:
        raise RuntimeError("survivor lacks frozen wing-size value")
    return float(np.mean(values))


def habitat_value(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(
        str(row[key]).strip()
        for key in (
            "CanopyAffinity",
            "EdgeAffinity",
            "MoistureAffinity",
            "DisturbanceAffinity",
        )
    )


def zmap(values: dict[str, float], species: tuple[str, ...]) -> dict[str, float]:
    x = np.log1p(np.asarray([values[name] for name in species], dtype=float))
    sd = float(np.std(x))
    if sd <= np.sqrt(np.finfo(float).eps):
        z = np.zeros_like(x)
    else:
        z = (x - float(np.mean(x))) / sd
    return dict(zip(species, map(float, z)))


def pair_features(
    species: tuple[str, ...],
    source_names: np.ndarray,
    target_names: np.ndarray,
    trait_rows: dict[str, dict[str, str]],
) -> dict[str, np.ndarray]:
    wing = {name: wing_value(trait_rows[name]) for name in species}
    host = {
        name: float(trait_rows[name]["NumberOfHostplantFamilies"])
        for name in species
    }
    volt = {name: str(trait_rows[name]["Voltinism"]).strip() for name in species}
    habitat = {name: habitat_value(trait_rows[name]) for name in species}
    wz = zmap(wing, species)
    hz = zmap(host, species)

    wing_sim = np.asarray([
        np.exp(-abs(wz[source] - wz[target]))
        for source, target in zip(source_names, target_names)
    ])
    host_sim = np.asarray([
        np.exp(-abs(hz[source] - hz[target]))
        for source, target in zip(source_names, target_names)
    ])
    volt_sim = np.asarray([
        float(volt[source] == volt[target])
        for source, target in zip(source_names, target_names)
    ])
    habitat_sim = np.asarray([
        float(np.mean([
            a == b for a, b in zip(habitat[source], habitat[target])
        ]))
        for source, target in zip(source_names, target_names)
    ])
    return {
        "wing_similarity": wing_sim[:, None],
        "host_breadth_similarity": host_sim[:, None],
        "voltinism_similarity": volt_sim[:, None],
        "habitat_similarity": habitat_sim[:, None],
        "wing_plus_host": np.column_stack([wing_sim, host_sim]),
        "wing_host_joint_match": (wing_sim * host_sim)[:, None],
        "four_axes_separate": np.column_stack([
            wing_sim, host_sim, volt_sim, habitat_sim
        ]),
        "four_axes_plus_wing_host_interaction": np.column_stack([
            wing_sim, host_sim, volt_sim, habitat_sim, wing_sim * host_sim
        ]),
    }


def frozen_geometry(design, pairs: np.ndarray) -> np.ndarray:
    target = pairs[:, 0].astype(np.int64)
    source = pairs[:, 1].astype(np.int64)
    geometry = np.asarray(pairs[:, 3:7], dtype=float).copy()
    geometry[:, 1] = np.log1p(geometry[:, 1] / 500.0)
    geometry[:, 2] = np.log(geometry[:, 2])
    geometry[:, 3] = np.log(geometry[:, 3])
    return np.column_stack([
        geometry,
        np.asarray(design["geometry_kernel"], dtype=float)[target, source],
    ])


def fit_predict(x_train, y_train, x_test):
    center = np.mean(x_train, axis=0)
    scale = np.std(x_train, axis=0, ddof=0)
    scale[scale <= np.sqrt(np.finfo(float).eps)] = 1.0
    z_train = (x_train - center) / scale
    z_test = (x_test - center) / scale
    X = np.column_stack([np.ones(len(z_train)), z_train])
    Xt = np.column_stack([np.ones(len(z_test)), z_test])
    beta = np.linalg.lstsq(X, y_train, rcond=None)[0]
    return Xt @ beta


def correlation(a, b):
    aa = np.asarray(a, dtype=float) - float(np.mean(a))
    bb = np.asarray(b, dtype=float) - float(np.mean(b))
    denom = float(np.sqrt(np.dot(aa, aa) * np.dot(bb, bb)))
    return 0.0 if denom <= np.finfo(float).tiny else float(np.dot(aa, bb) / denom)


def compare_cv(
    response: np.ndarray,
    geometry: np.ndarray,
    feature: np.ndarray | None,
    folds: np.ndarray,
) -> dict:
    pred_geometry = np.empty(len(response), dtype=float)
    pred_augmented = np.empty(len(response), dtype=float)
    for fold in sorted(np.unique(folds)):
        train = folds != fold
        test = ~train
        pred_geometry[test] = fit_predict(
            geometry[train], response[train], geometry[test]
        )
        if feature is None:
            pred_augmented[test] = pred_geometry[test]
        else:
            augmented = np.column_stack([geometry, feature])
            pred_augmented[test] = fit_predict(
                augmented[train], response[train], augmented[test]
            )
    rmse_geometry = float(np.sqrt(np.mean((response - pred_geometry) ** 2)))
    rmse_augmented = float(np.sqrt(np.mean((response - pred_augmented) ** 2)))
    corr_geometry = correlation(response, pred_geometry)
    corr_augmented = correlation(response, pred_augmented)
    return {
        "geometry_only_rmse": rmse_geometry,
        "augmented_rmse": rmse_augmented,
        "rmse_improvement": rmse_geometry - rmse_augmented,
        "geometry_only_correlation": corr_geometry,
        "augmented_correlation": corr_augmented,
        "correlation_gain": corr_augmented - corr_geometry,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", type=Path, required=True)
    ap.add_argument("--survivor-design-npz", type=Path, required=True)
    ap.add_argument("--leptraits-csv", type=Path, required=True)
    ap.add_argument("--parent-primary-result", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    parent = json.loads(args.parent_primary_result.read_text())
    if parent.get("schema") != "ttf_lepidoptera_trait_gradient_empirical_primary_result_v0.1":
        raise RuntimeError("parent primary result schema drift")
    if parent.get("status") != "PRIMARY_ONE_SHOT_DECISION_COMPLETE":
        raise RuntimeError("parent primary result incomplete")
    if sha256_path(args.archive) != str(parent["source_archive_sha256"]):
        raise RuntimeError("archive SHA drift")
    if sha256_path(args.survivor_design_npz) != str(parent["survivor_design_npz_sha256"]):
        raise RuntimeError("survivor design SHA drift")

    design = np.load(args.survivor_design_npz, allow_pickle=False)
    species = tuple(map(str, design["species_order"]))
    coords = np.asarray(design["coordinates"], dtype=float)
    offsets = np.asarray(design["coordinate_offsets"], dtype=np.int64)
    expected = {
        name: coords[offsets[i]:offsets[i + 1]]
        for i, name in enumerate(species)
    }

    with tempfile.TemporaryDirectory(prefix="ttf_lepidoptera_pair_axis_") as td:
        root = selective_extract(args.archive, Path(td), set(species))
        genes = read_genes_rows(root / "genes.txt")
        universe = {str(row.get("species", "")).strip() for row in genes}
        scan = scan_phylogatr_phase1(
            root,
            aliases=ALIASES,
            excluded_species=universe - set(species),
            min_localities=12,
            min_endpoint_training_edges=5,
            neighbor_fraction=0.15,
        )
        panels = {p.species: p for p in choose_one_panel_per_species(scan.candidates)}
        if set(panels) != set(species):
            raise RuntimeError("survivor panel reconstruction drift")
        for name in species:
            if not np.array_equal(panels[name].geometry.coordinates, expected[name]):
                raise RuntimeError(f"survivor geometry drift for {name}")

        response_by_species = {}
        with ProcessPoolExecutor(max_workers=max(1, int(args.workers))) as executor:
            futures = [
                executor.submit(response_item, (name, panels[name]))
                for name in species
            ]
            for future in as_completed(futures):
                name, value = future.result()
                response_by_species[name] = value

    pairs, target_index, source_index, _, pair_id = frozen_pair_surface(design)
    pair_response = pair_transfer_congruence(
        response_by_species,
        species,
        pairs,
        pair_id,
        np.asarray(design["alignment_target"], dtype=np.int64),
        np.asarray(design["alignment_source"], dtype=np.int64),
    )
    source_names = np.asarray([species[i] for i in source_index], dtype=object)
    target_names = np.asarray([species[i] for i in target_index], dtype=object)
    geometry = frozen_geometry(design, pairs)
    trait_rows = load_traits(args.leptraits_csv)
    features = pair_features(species, source_names, target_names, trait_rows)
    folds = deterministic_pair_folds(
        source_names,
        target_names,
        folds=10,
        tag="lepidoptera-pair-axis-exploratory-v0.1",
    )

    baseline = compare_cv(pair_response, geometry, None, folds)
    models = {
        name: compare_cv(pair_response, geometry, value, folds)
        for name, value in features.items()
    }

    output = {
        "schema": "ttf_lepidoptera_pair_axis_exploratory_v0.1",
        "status": "POST_PRIMARY_EXPLORATORY_ONLY",
        "parent_primary_decision": parent["decision"],
        "pairs": int(len(pair_response)),
        "source_species": int(len(set(source_names))),
        "target_species": int(len(set(target_names))),
        "baseline": baseline,
        "models": models,
        "claim_boundary": {
            "confirmatory": False,
            "p_values_computed": False,
            "primary_result_rewritten": False,
            "purpose": "Generate candidate relational predictors for an independent future confirmatory panel.",
        },
        "serialized_sequence_identity": False,
        "serialized_edge_genetic_distance_vectors": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print(json.dumps({"baseline": baseline, "models": models}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
