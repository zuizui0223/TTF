"""Synthetic-only tests of Stage22 positional predictor source gate.

No original range files, source coordinates, genetic distances, or raw AUC.
"""
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from pyproj import CRS, Transformer
from shapely.geometry import LineString, box
from shapely.ops import transform
import shapely

import stage22_jittered_predictor_support as stage22


class Stage22Tests(unittest.TestCase):
    def test_geometry_hashes_are_not_optional(self):
        xy = np.zeros((76700, 2), float)
        cases = np.zeros(76700, dtype="u1")
        cases[:10013] = 1
        sp = np.zeros(76700, dtype="i4")
        records = [{"split_role": i % 2} for i in range(80)]
        with self.assertRaisesRegex(ValueError, "Stage20 coordinate identity STOP"):
            stage22.ensure_frozen_geometry(xy, cases, sp, records, {})

    def test_moved_cases_are_deterministic_background_unchanged(self):
        line = LineString([(0.0, 0.0), (0.5, 0.0)])
        polygon = box(-2.0, -2.0, 2.0, 2.0)
        center = polygon.representative_point()
        proj = CRS.from_proj4(
            f"+proj=aeqd +lat_0={center.y} +lon_0={center.x} "
            "+datum=WGS84 +units=m +no_defs"
        )
        to_xy = Transformer.from_crs(4326, proj, always_xy=True)
        to_ll = Transformer.from_crs(proj, 4326, always_xy=True)
        pline = transform(to_xy.transform, line)
        poly = transform(to_xy.transform, polygon)
        d = np.arange(5000.0, pline.length, 10000.0)
        pts = shapely.line_interpolate_point(pline, d)
        pts = pts[shapely.covers(poly, pts)]
        local = shapely.get_coordinates(pts)
        lon, lat = to_ll.transform(local[:, 0], local[:, 1])
        cases = np.array([1] * len(local) + [0], dtype="u1")
        xy = np.vstack([np.column_stack([lon, lat]), [[1.0, 0.0]]])
        sp = np.zeros(len(xy), dtype="i4")
        record = {
            "common_name": "Test",
            "scientific_name": "Testus syntheticus",
            "case_points": len(local),
            "split_role": 0,
        }
        contract = {
            "status": "LOCKED_PRE_EMPIRICAL_BREAK_AUC",
            "case_jitter": {
                "seed_namespace": "ttf-stage21-jensen50km-jitter|",
                "source_uncertainty_radius_km": 50,
                "max_attempts": 20,
                "nonzero_required_fraction": 0.9,
                "global_nonzero_median_km_at_least": 10,
            }
        }
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            for name in ("breaks.zip", "ranges.zip"):
                with zipfile.ZipFile(work / name, "w") as z:
                    z.writestr("placeholder.txt", "synthetic")
            shared = {"test": ("dummy", "dummy")}
            with patch.object(stage22, "members", return_value=shared), patch.object(
                stage22, "load_group", side_effect=[
                    (line, None), (polygon, None), (line, None), (polygon, None)
                ]
            ):
                moved, rows, stats = stage22.jitter_cases_from_verified_sources(
                    xy, cases, sp, [record], work, contract
                )
                again, again_rows, again_stats = stage22.jitter_cases_from_verified_sources(
                    xy, cases, sp, [record], work, contract
                )
        np.testing.assert_array_equal(moved, again)
        self.assertEqual(rows, again_rows)
        self.assertEqual(stats, again_stats)
        np.testing.assert_array_equal(moved[cases == 0], xy[cases == 0])
        self.assertGreater(stats["fraction_shifted"], 0.9)
        self.assertLessEqual(stats["shift_p95_km"], 50.0)
        self.assertEqual(len(rows), 1)

    def test_wrong_jitter_rule_stops_before_any_predictor(self):
        with self.assertRaisesRegex(ValueError, "50-km/20-draw"):
            stage22.jitter_cases_from_verified_sources(
                np.zeros((1, 2)), np.array([1], dtype="u1"),
                np.array([0]), [], Path("."), {
                    "status": "LOCKED_PRE_EMPIRICAL_BREAK_AUC",
                    "case_jitter": {
                        "source_uncertainty_radius_km": 49,
                        "max_attempts": 20,
                    },
                },
            )


if __name__ == "__main__":
    unittest.main()
