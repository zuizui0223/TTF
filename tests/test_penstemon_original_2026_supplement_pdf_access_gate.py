"""PDF-source-only preflight with in-memory synthetic PDF, no plant outcomes."""
import importlib.util
from pathlib import Path
from io import BytesIO
import pytest
from pypdf import PdfWriter

FILE=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_2026_supplement_pdf_access_gate_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_pdf_gate",FILE)
p=importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def synthetic_pdf():
    f=BytesIO()
    writer=PdfWriter()
    writer.add_blank_page(width=612,height=792)
    writer.write(f)
    return f.getvalue()


def test_pdf_signature_and_metadata_safe():
    output=p.inspect_pdf_bytes(synthetic_pdf())
    assert output["page_count"]==1
    assert output["status"]=="ORIGINAL_SUPPLEMENT_PDF_CONTENT_VERIFIED_INDEX_ONLY"
    assert output["no_original_specimen_rows_or_phenotype_values_exported"] is True
    assert output["source_pdf_not_equivalent_to_34_phenotype_specimen_manifest"] is True
    assert output["table_index_mentions_page_number_1_based"]=={
        "Table S1":[],"Table S2":[],"Table S3":[]
    }


def test_html_response_cannot_be_called_original_pdf():
    with pytest.raises(ValueError,match="PDF signature"):
        p.inspect_pdf_bytes(b"<html>Table S1</html>")


def test_truncated_or_missing_source_is_not_biological_result():
    assert "PMC_supplement_generic" in p.SOURCES
    assert "bioRxiv_original_supplement" in p.SOURCES
    assert p.MAX_PDF_BYTES<=20_000_000
