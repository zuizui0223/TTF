"""Original PubMed Central supplementary PDF source integrity only."""
import importlib.util
from pathlib import Path
from io import BytesIO
from pypdf import PdfWriter
import pytest

P=Path(__file__).resolve().parents[1]/"research"/"prospective"/"penstemon_original_europepmc_table_s1_index_20261115.py"
spec=importlib.util.spec_from_file_location("penstemon_s1",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def make_pdf():
    w=PdfWriter()
    w.add_blank_page(width=612,height=792)
    buff=BytesIO()
    w.write(buff)
    return buff.getvalue()


def test_non_source_fake_pdf_rejected_on_original_size():
    with pytest.raises(ValueError,match="not original"):
        m.index_pdf(make_pdf())


def test_source_archive_sha_and_size_pin_are_immutable():
    assert m.EARLIER_ZIP_SHA=="42abf929e1428e9f643e76ed80057574b66c4bcd2c47b77ff50ef370c4a01278"
    assert m.PDF_MEMBER=="media-1.pdf"
    assert m.PDF_DECLARED_BYTES==2808017
