"""Source-only synthetic Dryad CSV header and ZIP central-directory tests."""
import io,unittest,zipfile
from dryad_headers_v01 import inspect_verified_source,classify_header

class HeaderTests(unittest.TestCase):
    def test_first_csv_line_only(self):
        raw=b"species,study,source_lat,source_lon,site_lat,site_lon,fitness\nA,S,1,2,3,4,9999999\n"
        result=inspect_verified_source(raw,"csv")
        self.assertEqual(result["headers"][0],"species")
        self.assertNotIn("9999999",str(result))
        self.assertFalse(result["row_2_or_later_fields_parsed"])
        names=classify_header(result["headers"])
        self.assertTrue(names["possible_origin_latlon_fields"])
        self.assertTrue(names["possible_garden_latlon_fields"])
    def test_archive_members_without_content(self):
        f=io.BytesIO()
        with zipfile.ZipFile(f,"w") as z:
            z.writestr("source/locations.csv","lat,raw_fitness\n42,7777777\n")
        x=inspect_verified_source(f.getvalue(),"zip")
        self.assertEqual(x["member_names"],["source/locations.csv"])
        self.assertNotIn("7777777",str(x))
        self.assertFalse(x["archive_member_contents_read"])
    def test_invalid_archive(self):
        with self.assertRaises(ValueError):inspect_verified_source(b"invalid","zip")
    def test_generic_lat_lon_does_not_prove_origin_or_garden(self):
        x=classify_header(["Study","Species","Latitude","Longitude","fitness"])
        self.assertFalse(x["possible_origin_latlon_fields"])
        self.assertFalse(x["possible_garden_latlon_fields"])

if __name__=="__main__":unittest.main()
