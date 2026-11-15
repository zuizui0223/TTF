import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from primary_study_overlap_gate import FIELDS, audit_overlap, norm_doi, read_manifest


class MetadataOnlyOverlapGateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def manifest(self, file_name, rows, fields=FIELDS):
        path = self.root / file_name
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(fields)
            writer.writerows(rows)
        return read_manifest(path)

    def row(self, dataset, study, doi, author="Jones", year="2018", title="Transplant alpha", taxon="Plant species"):
        return (dataset, study, doi, author, year, title, taxon)

    def test_doi_alias_normalized(self):
        self.assertEqual(norm_doi("https://doi.org/10.1111/ABC"), "10.1111/abc")
        self.assertEqual(norm_doi("doi:10.1111/abc"), "10.1111/abc")

    def test_shared_primary_doi_overrides_distinct_dataset_dois(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A1", "10.1111/x")])
        b = self.manifest("b.csv", [self.row("10.5061/dryad.b", "B9", "https://doi.org/10.1111/X")])
        out = audit_overlap(a, b)
        self.assertEqual(out["confirmed_bibliographic_overlap"], 1)
        self.assertEqual(out["verified_independent_primary_studies"], 0)

    def test_title_year_duplicate_without_paper_doi(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A", "", title="Same Experiment")])
        b = self.manifest("b.csv", [self.row("10.5061/dryad.b", "B", "", title="Same--experiment!")])
        self.assertEqual(audit_overlap(a, b)["confirmed_bibliographic_overlap"], 1)

    def test_author_year_needs_manual_review(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A", "10.1111/a", title="Experiment one")])
        b = self.manifest("b.csv", [self.row("10.5061/dryad.b", "B", "10.1111/b", title="Experiment two")])
        self.assertEqual(audit_overlap(a, b)["possible_primary_study_overlap_needs_review"], 1)

    def test_different_authors_same_taxon_not_certified(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A", "10.1111/a")])
        b = self.manifest("b.csv", [self.row("10.5061/dryad.b", "B", "10.1111/b", "Smith", "2021", "Different study")])
        out = audit_overlap(a, b)
        self.assertEqual(out["unlinked_but_not_independence_certified"], 1)
        self.assertEqual(out["verified_independent_primary_studies"], 0)

    def test_repeated_taxa_do_not_inflate_study_count(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A", "10.1111/a"), self.row("10.5061/dryad.a", "A", "10.1111/a", taxon="Second species")])
        self.assertEqual(len(a["studies"]), 1)
        self.assertEqual(len(a["studies"]["A"]["taxa"]), 2)

    def test_fitness_column_rejected_before_row_read(self):
        path = self.root / "unsafe.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(FIELDS + ("fitness",))
            writer.writerow(self.row("10.5061/dryad.a", "A", "10.1111/x") + ("0.25",))
        with self.assertRaisesRegex(ValueError, "metadata-only header"):
            read_manifest(path)

    def test_missing_identifier_rejected(self):
        path = self.root / "missing.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(FIELDS)
            writer.writerow(self.row("10.5061/dryad.a", "A", "", author="", title=""))
        with self.assertRaisesRegex(ValueError, "Primary DOI"):
            read_manifest(path)

    def test_conflicting_study_metadata_rejected(self):
        path = self.root / "conflict.csv"
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(FIELDS)
            writer.writerow(self.row("10.5061/dryad.a", "A", "10.1111/a"))
            writer.writerow(self.row("10.5061/dryad.a", "A", "10.1111/b"))
        with self.assertRaisesRegex(ValueError, "Inconsistent bibliographic"):
            read_manifest(path)

    def test_same_dataset_disallowed(self):
        a = self.manifest("a.csv", [self.row("10.5061/dryad.a", "A", "10.1111/a")])
        with self.assertRaisesRegex(ValueError, "same source dataset"):
            audit_overlap(a, a)


if __name__ == "__main__":
    unittest.main()
