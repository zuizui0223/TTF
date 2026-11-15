"""Synthetic-only tests for exact-mirror transplant CSV headers, never record values."""
import io
import unittest
import zipfile
from bontrager_member_headers_v01 import csv_headers_from_archive,codebook_excerpts


class OriginalHeaderFirewall(unittest.TestCase):
    def test_first_row_only_in_selected_member(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:
            z.writestr('dryad/data/tp_locs.csv',
                       'study,site,lat,long,fitness\nS,A,1,2,999999\n')
        r=csv_headers_from_archive(b.getvalue(),['dryad/data/tp_locs.csv'])
        self.assertEqual(r[0]['headers'],['study','site','lat','long','fitness'])
        self.assertNotIn('999999',str(r))
        self.assertFalse(r[0]['source_rows_read'])
        self.assertFalse(r[0]['fitness_values_read'])

    def test_missing_selected_member_abstains(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:
            z.writestr('another.csv','study,fitness\nS,999999\n')
        with self.assertRaises(ValueError):
            csv_headers_from_archive(b.getvalue(),['dryad/data/tp_locs.csv'])

    def test_readme_only_short_codebook_lines(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:
            z.writestr('dryad/README.txt',
                       'Files describe origin location and transplant site.\nUnrelated sentence\n')
        x=codebook_excerpts(b.getvalue(),'dryad/README.txt',1000)
        self.assertEqual(x[0]['line'],1)
        self.assertNotIn('Unrelated sentence',str(x))


if __name__=="__main__":
    unittest.main()
