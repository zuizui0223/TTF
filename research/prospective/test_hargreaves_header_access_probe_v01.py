"""Synthetic safety tests for the public Dryad source-only header probe."""
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hargreaves_header_access_probe_v01 import (
    MAX_HEADER_BYTES, URLS, parse_header_bytes, probe,
)


class SourceHeaderProbeTests(unittest.TestCase):
    def test_csv_header_only_ignores_fitness_row(self):
        result = parse_header_bytes(
            b'study,taxon,fitness\noriginal_study,plant,99.99\n'
        )
        self.assertEqual(result, ["study", "taxon", "fitness"])
        self.assertNotIn("99.99", result)

    def test_html_challenge_fails_closed(self):
        with self.assertRaises(ValueError):
            parse_header_bytes(b"<!DOCTYPE html>\n")

    def test_duplicate_column_fails(self):
        with self.assertRaises(ValueError):
            parse_header_bytes(b"study,study\n")

    def test_oversize_header_fails(self):
        with self.assertRaises(ValueError):
            parse_header_bytes(b"x" * (MAX_HEADER_BYTES + 1))

    def test_both_403_fail_access_without_claiming_overlap(self):
        from urllib.error import HTTPError
        def forbidden(request, timeout):
            raise HTTPError(request.full_url, 403, "Forbidden", None, None)
        with patch("hargreaves_header_access_probe_v01.urlopen", side_effect=forbidden):
            r = probe()
        self.assertEqual(r["status"], "SOURCE_ACCESS_GATED_NO_RECORD_LEVEL_OVERLAP_CLAIM")
        self.assertFalse(r["source_study_overlap_checked"])
        self.assertEqual(r["records_claimed_independent"], 0)
        self.assertEqual(len(r["attempts"]), len(URLS))


if __name__ == "__main__":
    unittest.main()
