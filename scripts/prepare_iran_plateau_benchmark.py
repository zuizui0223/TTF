#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re
import statistics
import xml.etree.ElementTree as ET
from zipfile import ZipFile

PMCID = "PMC13126619"
DOI = "10.1111/mec.70355"
S001 = "MEC-35-e70355-s001.zip"
S002 = "MEC-35-e70355-s002.zip"
EXPECTED_SPECIES = (
    "Allium scabriscapum",
    "Crepis heterotricha",
    "Didymophysa aucheri",
    "Dielsiocharis kotschyi",
    "Helichrysum oligocephalum",
    "Onosma microcarpa",
    "Phlomis olivieri",
    "Physoptychis gnaphalodes",
    "Tanacetum kotschyi",
)
IUPAC_DIPLOID = {
    "A": ("A", "A"),
    "C": ("C", "C"),
    "G": ("G", "G"),
    "T": ("T", "T"),
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
}
XLSX_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def digest_bytes(data: bytes) -> str:
    return sha256(data).hexdigest()


def _column_index(cell_ref: str) -> int:
    match = re.match(r"([A-Z]+)", cell_ref)
    if not match:
        raise ValueError(f"bad spreadsheet cell reference: {cell_ref}")
    out = 0
    for char in match.group(1):
        out = 26 * out + (ord(char) - 64)
    return out - 1


def read_first_sheet_xlsx(data: bytes) -> list[tuple[int, list[object | None]]]:
    """Read the simple Table S1 workbook with only the Python standard library."""
    with ZipFile(BytesIO(data)) as book:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in book.namelist():
            root = ET.fromstring(book.read("xl/sharedStrings.xml"))
            for item in root.findall("m:si", XLSX_NS):
                shared.append("".join(node.text or "" for node in item.findall(".//m:t", XLSX_NS)))
        root = ET.fromstring(book.read("xl/worksheets/sheet1.xml"))
        rows: list[tuple[int, list[object | None]]] = []
        for row in root.findall(".//m:sheetData/m:row", XLSX_NS):
            values: dict[int, object | None] = {}
            for cell in row.findall("m:c", XLSX_NS):
                index = _column_index(cell.attrib["r"])
                cell_type = cell.attrib.get("t")
                value_node = cell.find("m:v", XLSX_NS)
                if cell_type == "inlineStr":
                    value: object | None = "".join(
                        node.text or "" for node in cell.findall(".//m:t", XLSX_NS)
                    )
                elif value_node is None:
                    value = None
                else:
                    raw = value_node.text or ""
                    if cell_type == "s":
                        value = shared[int(raw)]
                    elif cell_type == "b":
                        value = raw == "1"
                    else:
                        try:
                            number = float(raw)
                            value = int(number) if number.is_integer() else number
                        except ValueError:
                            value = raw
                values[index] = value
            if values:
                dense: list[object | None] = [None] * (max(values) + 1)
                for index, value in values.items():
                    dense[index] = value
                rows.append((int(row.attrib["r"]), dense))
        return rows


def species_key(value: object) -> str:
    words = str(value).strip().split()
    if len(words) < 2:
        raise ValueError(f"cannot derive binomial species key from {value!r}")
    return f"{words[0]} {words[1]}"


def normalized_name(value: str) -> str:
    return re.sub(r"[^a-z]", "", value.lower())


def read_table_s1(bundle: Path) -> tuple[bytes, list[dict[str, object]]]:
    bundle_bytes = bundle.read_bytes()
    with ZipFile(BytesIO(bundle_bytes)) as archive:
        candidates = [name for name in archive.namelist() if name.lower().endswith("table_s1.xlsx")]
        if len(candidates) != 1:
            raise RuntimeError(f"expected one Table_S1.xlsx, found {candidates}")
        table_bytes = archive.read(candidates[0])
    rows = read_first_sheet_xlsx(table_bytes)
    if len(rows) < 3:
        raise RuntimeError("Table S1 has no data rows")
    header = rows[1][1]
    expected_header = [
        "Collection No.", "Manusc. No.", "Species", "Location", "Date", "Altitude",
        "N", "E", "Nr. of samples",
    ]
    if header[: len(expected_header)] != expected_header:
        raise RuntimeError(f"unexpected Table S1 header: {header}")

    data: list[dict[str, object]] = []
    for row_number, row in rows[2:]:
        if len(row) < 9 or row[0] is None:
            continue
        record = {
            "row": row_number,
            "collection_no": str(int(row[0])),
            "manuscript_no": None if row[1] is None else int(row[1]),
            "species": species_key(row[2]),
            "location": str(row[3]),
            "date": None if row[4] is None else str(row[4]),
            "altitude_m": float(row[5]),
            "latitude": float(row[6]),
            "longitude": float(row[7]),
            "n_samples": int(row[8]),
        }
        if not (-90 <= record["latitude"] <= 90 and -180 <= record["longitude"] <= 180):
            raise RuntimeError(f"invalid coordinate at Table S1 row {row_number}")
        data.append(record)
    return table_bytes, data


def parse_nexus(text: str) -> tuple[int, int, list[str], dict[str, str]]:
    ntax_match = re.search(r"NTAX\s*=\s*(\d+)", text, re.I)
    nchar_match = re.search(r"NCHAR\s*=\s*(\d+)", text, re.I)
    labels_match = re.search(r"TAXLABELS\s*(.*?)\s*;", text, re.I | re.S)
    matrix_match = re.search(r"\bMATRIX\b(.*?)\n\s*;\s*\n\s*END\s*;", text, re.I | re.S)
    if not all((ntax_match, nchar_match, labels_match, matrix_match)):
        raise RuntimeError("NEXUS file is missing TAXA/CHARACTERS fields")
    ntax = int(ntax_match.group(1))
    nchar = int(nchar_match.group(1))
    labels = [line.strip() for line in labels_match.group(1).splitlines() if line.strip()]
    label_set = set(labels)
    sequences: dict[str, str] = {}
    current: str | None = None
    for raw in matrix_match.group(1).splitlines():
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0] in label_set:
            current = parts[0]
            sequences[current] = "".join(parts[1:]).upper()
        else:
            if current is None:
                raise RuntimeError("NEXUS sequence continuation before first taxon")
            sequences[current] += "".join(parts).upper()
    if ntax != len(labels) or set(sequences) != label_set:
        raise RuntimeError("NEXUS taxon count/labels do not agree")
    return ntax, nchar, labels, sequences


def collection_number(label: str) -> str:
    match = re.search(r"(\d{3})", label)
    if not match:
        raise ValueError(f"cannot find Collection No. in taxon label {label!r}")
    return match.group(1)


def audit_nexus_species(
    filename: str,
    text: str,
    table_rows: list[dict[str, object]],
) -> dict[str, object]:
    ntax, nchar, labels, sequences = parse_nexus(text)
    table_counts = Counter({str(row["collection_no"]): int(row["n_samples"]) for row in table_rows})
    nexus_counts = Counter(collection_number(label) for label in labels)
    sequence_lengths = Counter(len(sequence) for sequence in sequences.values())
    char_counts = Counter("".join(sequences.values()))
    allowed = set(IUPAC_DIPLOID) | {"?", "-"}
    invalid = {key: value for key, value in char_counts.items() if key not in allowed}
    missing = [sequence.count("?") / nchar for sequence in sequences.values()]

    multiallelic = 0
    all_missing = 0
    for index in range(nchar):
        alleles: set[str] = set()
        for sequence in sequences.values():
            code = sequence[index]
            if code in IUPAC_DIPLOID:
                alleles.update(IUPAC_DIPLOID[code])
        if not alleles:
            all_missing += 1
        elif len(alleles) > 2:
            multiallelic += 1

    passed = (
        ntax == sum(table_counts.values())
        and table_counts == nexus_counts
        and sequence_lengths == Counter({nchar: ntax})
        and not invalid
        and multiallelic == 0
        and all_missing == 0
    )
    return {
        "filename": filename,
        "ntax": ntax,
        "nchar": nchar,
        "table_population_count": len(table_counts),
        "nexus_population_count": len(nexus_counts),
        "table_individual_count": sum(table_counts.values()),
        "population_sample_counts_exact_match": table_counts == nexus_counts,
        "sequence_lengths": {str(key): value for key, value in sequence_lengths.items()},
        "invalid_iupac_codes": invalid,
        "multiallelic_site_count": multiallelic,
        "all_missing_site_count": all_missing,
        "individual_missing_fraction": {
            "min": min(missing),
            "median": statistics.median(missing),
            "max": max(missing),
        },
        "passed": passed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare the Iranian Plateau TTF benchmark without using outcomes.")
    parser.add_argument("--supplement-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    s001 = args.supplement_dir / S001
    s002 = args.supplement_dir / S002
    if not s001.exists() or not s002.exists():
        raise RuntimeError("required Iranian Plateau supplement bundles are missing")

    table_bytes, table_rows = read_table_s1(s001)
    by_species: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in table_rows:
        by_species[str(row["species"])].append(row)
    if tuple(sorted(by_species)) != tuple(sorted(EXPECTED_SPECIES)):
        raise RuntimeError(f"unexpected Table S1 species: {sorted(by_species)}")
    collection_ids = [str(row["collection_no"]) for row in table_rows]
    if len(collection_ids) != len(set(collection_ids)):
        raise RuntimeError("Collection No. is not globally unique in Table S1")

    with ZipFile(s002) as archive:
        nexus_names = sorted(name for name in archive.namelist() if name.lower().endswith(".nex"))
        if len(nexus_names) != len(EXPECTED_SPECIES):
            raise RuntimeError(f"expected nine NEXUS files, found {nexus_names}")
        lookup = {normalized_name(Path(name).stem): name for name in nexus_names}
        species_audits = []
        for species in EXPECTED_SPECIES:
            key = normalized_name(species)
            if key not in lookup:
                raise RuntimeError(f"no NEXUS file matched {species}")
            name = lookup[key]
            text = archive.read(name).decode("utf-8")
            audit = audit_nexus_species(name, text, by_species[species])
            audit["species"] = species
            species_audits.append(audit)

    geometry = {
        "schema": "ttf_iran_plateau_geometry_v0.1",
        "source": {
            "pmcid": PMCID,
            "doi": DOI,
            "table_s1_sha256": digest_bytes(table_bytes),
            "supplement_s001_sha256": digest_bytes(s001.read_bytes()),
            "outcome_values_used": False,
            "boundary_labels_used": False,
        },
        "species": [
            {
                "species": species,
                "populations": [
                    {
                        "collection_no": row["collection_no"],
                        "manuscript_no": row["manuscript_no"],
                        "latitude": row["latitude"],
                        "longitude": row["longitude"],
                        "altitude_m": row["altitude_m"],
                        "n_samples": row["n_samples"],
                    }
                    for row in by_species[species]
                ],
            }
            for species in EXPECTED_SPECIES
        ],
    }
    audit = {
        "schema": "ttf_iran_plateau_input_audit_v0.1",
        "source": geometry["source"],
        "table": {
            "species_count": len(by_species),
            "population_count": len(table_rows),
            "individual_count": sum(int(row["n_samples"]) for row in table_rows),
            "collection_numbers_unique": len(collection_ids) == len(set(collection_ids)),
        },
        "nexus": species_audits,
        "passed": all(item["passed"] for item in species_audits),
    }
    if not audit["passed"]:
        raise RuntimeError("Iranian Plateau input audit failed")

    (args.output_dir / "geometry.json").write_text(json.dumps(geometry, indent=2, sort_keys=True) + "\n")
    (args.output_dir / "input_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "passed": audit["passed"],
        "species": audit["table"]["species_count"],
        "populations": audit["table"]["population_count"],
        "individuals": audit["table"]["individual_count"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
