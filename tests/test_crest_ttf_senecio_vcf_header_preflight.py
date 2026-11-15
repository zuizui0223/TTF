"""Offline synthetic-only VCF header and source hash gate checks."""
import importlib.util
import io
import tarfile
from pathlib import Path

FILE = Path(__file__).resolve().parents[1] / "research" / "prospective" / "ttf_crest_senecio_vcf_header_preflight_20261115.py"
spec = importlib.util.spec_from_file_location("senecio_vcf_header", FILE)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


def fake_archive(samples=("D00-001", "H00_003", "D14-01", "H15-02"), duplicate=False):
    all_samples = tuple(samples) + ((samples[0],) if duplicate else ())
    header = "##fileformat=VCFv4.2\n##source=synthetic\n" + \
        "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",*all_samples]) + "\n"
    # Invalid-body sentry: a header-only validator must not decode any records.
    body = b"not_a_VCF_variant_record\n" * 12
    raw = header.encode() + body
    f = io.BytesIO()
    with tarfile.open(fileobj=f, mode="w:gz") as tf:
        info = tarfile.TarInfo("synthetic.vcf")
        info.size = len(raw)
        tf.addfile(info, io.BytesIO(raw))
    return f.getvalue()


def with_source_sha(data: bytes):
    old = v.GIT_BLOB_SHA
    v.GIT_BLOB_SHA = v.git_blob_sha(data)
    try: return v.audit_archive(data)
    finally: v.GIT_BLOB_SHA = old


def test_header_only_counts_and_unverified_flags():
    result = with_source_sha(fake_archive())
    assert result["vcf_sample_count"] == 4
    assert result["vcf_meta_header_lines"] == 2
    assert result["population_sample_counts"] == {"D00": 1, "D14": 1, "H00": 1, "H15": 1}
    assert result["genotypes_read"] is False
    assert result["snp_variant_values_read"] is False
    assert result["phenotype_genotype_same_individuals"] == "NOT_VERIFIED"


def test_duplicate_sample_fails_closed():
    try:
        with_source_sha(fake_archive(duplicate=True))
    except ValueError as e:
        assert "duplicate VCF sample" in str(e)
    else: raise AssertionError("duplicate sample accepted")


def test_changed_blob_fails_closed():
    try: v.audit_archive(fake_archive())
    except ValueError as e:
        assert "Git blob SHA mismatch" in str(e)
    else: raise AssertionError("changed original accepted")


def test_non_header_line_prevents_variant_read():
    try:
        v.read_vcf_header_only(io.BytesIO(b"##fileformat=VCFv4.2\n1\t2\t3\n"))
    except ValueError as e:
        assert "unexpected non-header" in str(e)
    else: raise AssertionError("variant before header accepted")


def test_appledouble_resource_fork_does_not_duplicate_vcf():
    data = fake_archive()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as old:
        original = old.extractfile("synthetic.vcf").read()
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as tf:
        for name, payload in (
            ("._synthetic.vcf", b"AppleDouble resource fork"),
            ("synthetic.vcf", original),
        ):
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tf.addfile(info, io.BytesIO(payload))
    result = with_source_sha(stream.getvalue())
    assert result["vcf_sample_count"] == 4
    assert result["archive_member_name"] == "synthetic.vcf"
