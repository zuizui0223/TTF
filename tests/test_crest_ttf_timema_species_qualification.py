"""Synthetic, response-value-blind tests for Timema TTF-Q support geometry."""
from pathlib import Path
import importlib.util
import io
import zipfile

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"ttf_crest_timema_species_qualification_20261115.py"
spec=importlib.util.spec_from_file_location("ttf_timema_qual",FILE)
s=importlib.util.module_from_spec(spec)
spec.loader.exec_module(s)


def rows():
    out=[]
    n=0
    for prefix in ("ABC","DEF"):
        for p in prefix:
            n+=1;out.append({"pair_id":str(n),"species_1":p,"species_2":p,"population_pair":f"{p}0-{p}1"})
        for a,b in ((prefix[0],prefix[1]),(prefix[1],prefix[2])):
            n+=1;out.append({"pair_id":str(n),"species_1":a,"species_2":b,"population_pair":a+"0-"+b+"0"})
    return out


def test_components_two_independent_comparison_modules():
    a=s.connected_components(rows())
    assert len(a)==2
    assert all(x["species_count"]==3 for x in a)
    assert all(x["within_species_dyads"]==3 and x["between_species_dyads"]==2 for x in a)


def test_both_isolation_categories_have_disjoint_coverage():
    x=s.split_support(rows())
    assert x["status"]=="GEOMETRY_ONLY"
    assert x["source_geometry_qualified_partitions"]>=1
    b=x["best_both_types_supported_partition"]
    assert len(b["train_species"])==3 and len(b["test_species"])==3
    assert set(b["train_species"]).isdisjoint(b["test_species"])
    assert b["train"]["within_species"]>=3 and b["test"]["within_species"]>=3
    assert b["train"]["between_species"]>=2 and b["test"]["between_species"]>=2


def synthetic_xlsx():
    w=b'<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Figure3" sheetId="1" r:id="rId1"/></sheets></workbook>'
    rel=b'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>'
    sh=b'<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c r="K1" t="inlineStr"><is><t>Habitat Isolation</t></is></c></row><row r="2"><c r="K2"><v>0.123456789</v></c><c r="L2"><v>-0.192837465</v></c></row><row r="3"><c r="K3"><v>0.9</v></c></row></sheetData></worksheet>'
    b=io.BytesIO()
    with zipfile.ZipFile(b,"w") as z:
        z.writestr("xl/workbook.xml",w)
        z.writestr("xl/_rels/workbook.xml.rels",rel)
        z.writestr("xl/worksheets/sheet1.xml",sh)
    return b.getvalue()


def test_outcome_presence_does_not_decode_effect_sizes():
    present=s.source_response_cell_presence(synthetic_xlsx())
    assert present["habitat_isolation"]==2
    assert present["sexual_isolation"]==1
    assert present["both_cells_structurally_present"]==1
    assert present["numeric_literal_storage_cells_not_validated_trials"] == {
        "habitat_isolation":2,"sexual_isolation":1
    }
    assert "0.123456789" not in repr(present)
    assert "-0.192837465" not in repr(present)
