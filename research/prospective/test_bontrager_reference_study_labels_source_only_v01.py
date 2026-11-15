import csv
import io
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bontrager_reference_study_labels_source_only_v01 import (
    EXPECTED, binomial, canonical, source_study_inventory,
)


class OriginalStudyLabelTests(unittest.TestCase):
    def synthetic(self, data, extra_header=None):
        required = ["study", "functional_group", "source", "site", "fitness_type", "taxon",
                    "fitness", "lat_source", "long_source", "lat_site", "long_site"]
        extras = [f"extra_{i}" for i in range(47 - len(required))]
        headers = required + extras
        if extra_header:
            headers[7] = extra_header
        out = io.StringIO()
        w = csv.writer(out)
        w.writerow(headers)
        for item in data:
            d = dict(zip(required, item))
            w.writerow([d.get(key, "") for key in headers])
        return out.getvalue().encode("utf-8")

    def row(self, source, lon, group="annual", study="Doe, Ecology, 2004", taxon="Plant species"):
        return [study, group, source, "Garden", "surv", taxon,
                "NEVER_CONSUME_RESPONSES", "35", str(lon), "36", "140"]

    def test_two_original_sources_one_garden_yield_one_cell(self):
        r = source_study_inventory(self.synthetic([self.row("A", 139), self.row("B", 141)]))
        self.assertEqual(r["distinct_qualified_original_study_labels"], 1)
        self.assertEqual(r["distinct_qualified_plant_species"], 1)
        self.assertEqual(r["qualified_study_species_garden_fitness_type_cells"], 1)
        self.assertEqual(r["catalog"][0]["candidate_year_from_label_NOT_verified"], "2004")
        self.assertEqual(r["outcome_columns_selected"], [])
        self.assertIsNone(r["catalog"][0]["primary_paper_doi"])

    def test_bad_nonplant_records_never_increase_cohort(self):
        raw = self.synthetic([self.row("A", 139), self.row("B", 141),
                              self.row("C", 142, group="arthropod")])
        self.assertEqual(source_study_inventory(raw)["qualified_study_species_garden_fitness_type_cells"], 1)

    def test_duplicate_origin_coordinates_not_two_independent_sources(self):
        raw = self.synthetic([self.row("A", 139), self.row("B", 139)])
        self.assertEqual(source_study_inventory(raw)["qualified_study_species_garden_fitness_type_cells"], 0)

    def test_conflicted_origin_id_excluded(self):
        raw = self.synthetic([self.row("A", 139), self.row("A", 140), self.row("B", 141)])
        self.assertEqual(source_study_inventory(raw)["distinct_qualified_original_study_labels"], 0)

    def test_frozen_geometry_gate_fails_closed(self):
        raw = self.synthetic([self.row("A", 139), self.row("B", 141)])
        with self.assertRaisesRegex(ValueError, "not reproduced"):
            source_study_inventory(raw, EXPECTED)

    def test_invalid_binomial_does_not_count(self):
        raw = self.synthetic([self.row("A", 139, taxon="sp."),
                              self.row("B", 141, taxon="sp.")])
        self.assertEqual(source_study_inventory(raw)["distinct_qualified_plant_species"], 0)
        self.assertIsNone(binomial("sp."))

    def test_no_invented_paper_identity_from_study_label(self):
        raw = self.synthetic([self.row("A", 139, study="Doe, Ecology, 2004"),
                              self.row("B", 141, study="Doe, Ecology, 2004")])
        a = source_study_inventory(raw)["catalog"][0]
        self.assertIsNone(a["primary_paper_doi"])
        self.assertIsNone(a["primary_paper_title"])
        self.assertIs(a["independent_original_experiment_certified"], False)
        self.assertEqual(canonical("Doe, Ecology, 2004"), "doe ecology 2004")

    def test_unexpected_source_header_rejected(self):
        raw = self.synthetic([self.row("A", 139), self.row("B", 141)],
                             extra_header="invalid_source_latitude")
        with self.assertRaisesRegex(ValueError, "Unrecognized original"):
            source_study_inventory(raw)


if __name__ == "__main__":
    unittest.main()
