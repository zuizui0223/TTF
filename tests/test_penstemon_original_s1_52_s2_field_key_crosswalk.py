"""Synthetic-only Penstemon Table S1 verified grammar and S2 field identity tests.

No flower measurements, source PDF or genetic outcomes opened by these tests.
"""
from pathlib import Path
import importlib.util

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_s1_52_s2_field_key_crosswalk_20261118.py"
spec=importlib.util.spec_from_file_location("penstemon_s1_s2_visual",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_S2_source_field_unhyphenated_vouchers_are_distinct_from_genomes():
    text=("Table S2. Populations phenotyped for floral dimensions in field\n"
          "Population Species Locality Collector Herbarium deposited\n"
          "CAW015 P. virgatus Victor, CO CAW KANU\n"
          "CAW017 P. virgatus Lake George, CO CAW KANU\n"
          "CAW030 P. barbatus Gila NF, NM CAW KANU\n")
    result=m.parse_s2_field_identity(text)
    assert result["status"]=="SOURCE_S2_FIELD_IDENTITIES_PARSED"
    assert result["source_rows"]==3
    assert result["original_field_voucher_ids"]==["CAW015","CAW017","CAW030"]
    assert result["original_field_species_name_prefixes"]==[
        "P. virgatus","P. virgatus","P. barbatus"]
    assert result["no_source_numeric_effects_read"] is True


def test_original_source_voucher_sep_normalization_does_not_invent_species():
    assert m.normalize_voucher("CAW-015")=="CAW015"
    assert m.normalize_voucher("CAW015")=="CAW015"
    assert m.normalize_voucher("TTC026503")=="TTC026503"
    assert m.normalize_voucher("CAW-190")!="CAW015"


def test_partial_S2_line_must_report_incomplete_parser():
    text=("Table S2. Field plants\n"
          "CAW015 P. virgatus Victor, CO\n"
          "CAW020\n")
    result=m.parse_s2_field_identity(text)
    assert result["status"]=="SOURCE_S2_IDENTITY_PARSER_INCOMPLETE"
    assert result["unparsed_id_prefix_line_count"]==1
