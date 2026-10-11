"""All mask source ZIP tests use synthetic fixture data, never phylogatR sequences."""
import hashlib
import stat
import zipfile

import pytest

from scripts.verify_and_extract_historical_host_mask_source import ARCHIVE_SHA, ARCHIVE_BYTES, extract_verified_archive


def make_zip(tmp_path, members):
    archive=tmp_path/"small.zip"
    with zipfile.ZipFile(archive,"w",compression=zipfile.ZIP_DEFLATED) as z:
        for name,value in members.items():
            z.writestr(name,value)
    raw=archive.read_bytes()
    return archive,hashlib.sha256(raw).hexdigest(),len(raw)


def test_exact_original_identity_and_root_required(tmp_path):
    a,sha,n=make_zip(tmp_path,{"root/genes.txt":"gene\n", "root/cite.txt":"citation\n",
                               "root/Insecta/Foo/Foo-COI.afa":">h1\nACGT\n"})
    root=extract_verified_archive(a,tmp_path/"destination",expected_sha=sha,expected_size=n)
    assert root.name=="root"
    assert (root/"genes.txt").is_file()


def test_archive_byte_mutation_fails_before_any_extraction(tmp_path):
    a,sha,n=make_zip(tmp_path,{"root/genes.txt":"gene\n","root/cite.txt":"citation\n"})
    a.write_bytes(a.read_bytes()+b"x")
    with pytest.raises(RuntimeError,match="identity mismatch"):
        extract_verified_archive(a,tmp_path/"dst",expected_sha=sha,expected_size=n)
    assert not (tmp_path/"dst").exists()


def test_zip_path_traversal_fails_closed(tmp_path):
    a,sha,n=make_zip(tmp_path,{"../outside.txt":"bad",
                               "root/genes.txt":"ok", "root/cite.txt":"ok"})
    with pytest.raises(RuntimeError,match="traversal"):
        extract_verified_archive(a,tmp_path/"dst",expected_sha=sha,expected_size=n)
    assert not (tmp_path/"outside.txt").exists()


def test_symlink_member_fails_closed(tmp_path):
    path=tmp_path/"sym.zip"
    info=zipfile.ZipInfo("root/source-link")
    info.create_system=3
    info.external_attr=(stat.S_IFLNK|0o777)<<16
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("root/genes.txt","gene")
        z.writestr("root/cite.txt","cite")
        z.writestr(info,"../elsewhere")
    data=path.read_bytes()
    with pytest.raises(RuntimeError,match="symlink"):
        extract_verified_archive(path,tmp_path/"dst",
                                 expected_sha=hashlib.sha256(data).hexdigest(),
                                 expected_size=len(data))


def test_multiple_phylogatr_roots_fail_closed(tmp_path):
    a,sha,n=make_zip(tmp_path,{"a/genes.txt":"a","a/cite.txt":"a",
                               "b/genes.txt":"b","b/cite.txt":"b"})
    with pytest.raises(RuntimeError,match="one phylogatR root"):
        extract_verified_archive(a,tmp_path/"dst",expected_sha=sha,expected_size=n)


def test_canonical_exact_source_identity_is_full_sha256():
    assert len(ARCHIVE_SHA)==64
    assert ARCHIVE_SHA=="5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5"
    assert ARCHIVE_BYTES==274_988_692
