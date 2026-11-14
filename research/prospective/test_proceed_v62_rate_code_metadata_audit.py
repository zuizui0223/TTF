"""Only synthetic public-code classification tests. No PROCEED rate rows."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from proceed_v62_rate_code_metadata_audit import code_signals


class CodeOnlyTests(unittest.TestCase):
    def test_frozen_code_terms_and_no_execution(self):
        b=b"# ecological research code\nx <- abs(mean1-mean2)/generations\nlog10(x)\na <- df$gtime\n"
        s=code_signals(b,["generations","mean1","mean2","log10","gtime"],max_lines=4)
        self.assertEqual(s["line_count"],4)
        self.assertEqual(s["number_of_source_code_lines_matching_declared_terms"],3)
        self.assertTrue(s["code_executed"] is False)
        self.assertTrue(s["genetic_trait_rate_rows_read"] is False)
        self.assertTrue(all("line_number" in z for z in s["first_controlled_term_matches"]))

    def test_static_signal_selection_does_not_depend_on_any_response_rate(self):
        x=code_signals(b"response <- log10(abs(a))/gtime\n",
                       ["log10","gtime","darwin"])
        y=code_signals(b"response <- log10(abs(a))/gtime\n",
                       ["log10","gtime","darwin"])
        self.assertEqual(x,y)


if __name__=="__main__":
    unittest.main()
