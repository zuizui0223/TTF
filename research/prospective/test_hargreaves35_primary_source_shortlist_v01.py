"""Source-only 10-paper DOI recovery shortlist integrity regression checks."""
import hashlib
import json
import re
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

class PrimarySourceShortlistTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index=json.loads((HERE/"hargreaves35_primary_article_recovery_shortlist_v01.json").read_text())
        cls.rbytes=(ROOT/"research/prospective/source_registry/hargreaves35_study_geometry_v01.json").read_bytes()
        cls.receipt=json.loads(cls.rbytes)

    def test_source_only_recovery_tied_to_immutable_35_studies(self):
        self.assertEqual(hashlib.sha256(self.rbytes).hexdigest(),
                         self.index["source_geometry_original_SHA256"])
        labels={s["original_study_label"] for s in self.receipt["studies"]}
        self.assertEqual(len(labels),35)
        chosen={r["original_Hargreaves_study_label"] for r in self.index["record_list"]}
        self.assertEqual(len(chosen),10)
        self.assertTrue(chosen.issubset(labels))

    def test_doi_metadata_not_computer_invented_geocoordinates(self):
        doi=set()
        for row in self.index["record_list"]:
            self.assertRegex(row["primary_article_doi"],r"^10\.\d{4,9}/\S+$")
            self.assertNotIn(row["primary_article_doi"],doi)
            doi.add(row["primary_article_doi"])
            self.assertTrue(row["first_author_year_journal_title_bibliographic_match"])
            self.assertFalse(row["primary_experiment_identity_certified"])
            self.assertFalse(row["historical_climate_contrast_qualified"])
            self.assertFalse(row["genetic_adaptation_plasticity_speciation_effect_opened"])
            self.assertIn("UNVERIFIED",row["original_source_geocoordinate_status"]) if "UNVERIFIED" in row["original_source_geocoordinate_status"] else self.assertEqual(row["original_source_geocoordinate_status"],"PUBLISHER_DISPLAY_COORDINATES_REQUIRE_VALIDATION")

    def test_opportunity_cell_totals_are_geometry_not_replication(self):
        counts={s["original_study_label"]:s["candidate_comparable_same_garden_cells"]
                for s in self.receipt["studies"]}
        rows=self.index["record_list"]
        self.assertEqual(self.index["summary"]["labeled_papers"],len(rows))
        self.assertEqual(self.index["summary"]["study_geometry_cells_sum_not_independent"],
                         sum(counts[r["original_Hargreaves_study_label"]] for r in rows))
        self.assertEqual(self.index["summary"]["study_geometry_cells_sum_not_independent"],106)
        self.assertEqual(self.index["summary"]["study_labels_with_original_source_lat_lon_verified"],0)
        self.assertEqual(self.index["summary"]["independent_experiments_certified"],0)

    def test_no_pr154_requalification(self):
        self.assertTrue(self.index["not_reused_to_rescue_pr154_Bontrager_44_of_50"])
        self.assertTrue(self.index["no_outcome_used_to_select_sources"])
        self.assertEqual(self.receipt["summary"]["verified_source_population_georeferenced_studies_in_Hargreaves_original"],0)

if __name__=="__main__":
    unittest.main()
