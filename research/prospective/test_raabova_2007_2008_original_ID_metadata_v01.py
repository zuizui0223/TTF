import csv,copy,hashlib,io,json,sys,unittest
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from raabova_2007_2008_original_ID_metadata_v01 import (
    CONTRACT, classify_code, census,
)
class RaabovaStudyIDSourceOnlyTests(unittest.TestCase):
    def setUp(self):
        self.rule=copy.deepcopy(json.loads(CONTRACT.read_text(encoding="utf-8")))
        self.row2007=self.rule["validated_primary_papers"]["target_2007"]["study_label"]
        self.row2008=self.rule["validated_primary_papers"]["different_2008"]["study_label"]
    def source(self,rows):
        field=["study","taxon","functional_group","source","site","fitness"]
        header=field+[f"extra_{i}" for i in range(47-len(field))]
        s=io.StringIO(newline="")
        w=csv.writer(s,lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow([row.get(k,"") for k in header])
        raw=s.getvalue().encode("cp1252")
        self.rule["original_data"]["original_SHA256"]=hashlib.sha256(raw).hexdigest()
        return raw
    def r(self,study,source,site,taxon="Aster amellus",fitness="NEVER_READ"):
        return dict(study=study,source=source,site=site,taxon=taxon,functional_group="herbaceous_perennial",fitness=fitness)
    def test_exact_code_match_and_2008_distinct_paper(self):
        fields=self.rule["exact_id_rule"]["allowed_known_2007_population_codes"]
        self.assertEqual(classify_code("K1",fields),"K1")
        self.assertEqual(classify_code("Raabova2007_S2",fields),"S2")
        self.assertEqual(classify_code("Raabova2007-S3",fields),"S3")
        self.assertIsNone(classify_code("BlahKS1",fields))
        self.assertEqual(classify_code("S1_K2",fields),"AMBIGUOUS_MULTIPLE_CODES")
        self.assertNotEqual(self.rule["validated_primary_papers"]["different_2008"]["doi"],self.rule["validated_primary_papers"]["target_2007"]["doi"])
    def test_mutually_distinct_source_studies_not_conflated(self):
        raw=self.source([
            self.r(self.row2007,"p2007_S1","p2007_S1"),
            self.r(self.row2007,"p2007_S2","p2007_S2"),
            self.r(self.row2008,"p2008_1","p2008_1",taxon="Aster amellus 6x")
        ])
        # Synthetic original row-count guard is expected to stop (11088 fixture requirement).
        with self.assertRaisesRegex(ValueError,"row count drift"):
            census(raw,self.rule)
    def test_invalid_primary_reference_pair_is_hard_stop(self):
        raw=self.source([self.r(self.row2007,"S1","S1")])
        self.rule["validated_primary_papers"]["sudova_2010"]["table1_reference_50_is"]=self.rule["validated_primary_papers"]["target_2007"]["doi"]
        with self.assertRaisesRegex(ValueError,"not 2008"):
            census(raw,self.rule)
    def test_invalid_fitness_in_allowlist_fails(self):
        raw=self.source([self.r(self.row2007,"S1","S1")])
        self.rule["original_data"]["columns"].append("fitness")
        with self.assertRaisesRegex(ValueError,"allowlist"):
            census(raw,self.rule)
    def test_original_sha_change_fails(self):
        raw=self.source([self.r(self.row2007,"S1","S1")])
        with self.assertRaisesRegex(ValueError,"SHA mismatch"):
            census(raw+b"tamper",self.rule)
if __name__=="__main__":unittest.main()
