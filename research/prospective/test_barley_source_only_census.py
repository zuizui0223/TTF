"""Synthetic-only QA: source metadata gate never needs experimental ARR values."""
from __future__ import annotations
import csv
import io
import json
from pathlib import Path
import unittest

from barley_source_only_census import source_census, git_blob_sha1, CONTRACT


RULES = json.loads(CONTRACT.read_text())
FIELDS = ["study", "taxon", "phylum", "ecosystem", "source_population",
          "latitude", "longitude", "acclimation_time",
          "number_of_populations", "ARR", "thermal_limit_1", "thermal_limit_2"]


def sample_csv(arr=("0.10", "0.90", "0.12"), marine=False, discord=False):
    rows = [
        ["S1", "Testus alpha", "Arthropoda", "terrestrial", "P1", "30", "-99", "within-generation", "2", arr[0], "40", "42"],
        ["S1", "Testus alpha", "Arthropoda", "terrestrial", "P2", "31", "-98", "within-generation", "2", arr[1], "10", "12"],
        ["S2", "Testus beta", "Arthropoda", "marine" if marine else "terrestrial", "P1", "32", "-97", "within-generation", "1", arr[2], "55", "60"],
    ]
    if discord:
        rows.append(["S1", "Testus alpha", "Arthropoda", "terrestrial", "P1", "34", "-99", "within-generation", "2", "0.17", "1", "2"])
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(FIELDS)
    writer.writerows(rows)
    return out.getvalue().encode()


class SourceOnlyCensusTests(unittest.TestCase):
    def test_source_eligibility_does_not_change_with_any_phenotype_values(self):
        a = source_census(sample_csv(), RULES)
        b = source_census(sample_csv(("10000", "-50000", "NoPhenotypeAvailable")), RULES)
        self.assertEqual(a, b)
        self.assertEqual(a["original_study_species_groups_with_two_georeferenced_populations"], 1)
        self.assertEqual(a["eligible_terrestrial_species_two_populations"], 1)
        self.assertEqual(a["phenotype_ARR_values_evaluated"], False)
        self.assertEqual(a["status"], "PILOT_ONLY_NO_MACRO_HISTORY_VS_PLASTICITY_CLAIM")

    def test_aquatic_sources_are_not_terrestrial_climate_sites(self):
        out = source_census(sample_csv(marine=True), RULES)
        self.assertEqual(out["other_ecosystem_rows"].get("marine"), 1)
        self.assertEqual(out["terrestrial_rows"], 2)

    def test_inconsistent_origin_source_population_reported(self):
        out = source_census(sample_csv(discord=True), RULES)
        self.assertEqual(out["conflicting_origin_coordinates_for_same_pop_key"], 1)
        self.assertFalse(out["response_opening_authorized_by_this_census"])

    def test_source_git_blob_function_is_content_addressed(self):
        data = b"abc\n"
        self.assertEqual(git_blob_sha1(data), "8baef1b4abc478178b004d62031cf7fe6db6f903")
        self.assertNotEqual(git_blob_sha1(data), git_blob_sha1(b"abc\r\n"))


if __name__ == "__main__":
    unittest.main()
