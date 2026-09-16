from ttf.direct_insdc_metadata import iter_genbank_metadata


def test_metadata_parser_ignores_origin_identity():
    prefix = '''LOCUS       X 10 bp DNA
ACCESSION   ABC123
  ORGANISM  Testus animalis
            Eukaryota; Metazoa; Animalia; Arthropoda.
FEATURES             Location/Qualifiers
     source          1..10
                     /lat_lon="35.0 N 135.0 E"
     gene            1..10
                     /gene="COX1"
ORIGIN
'''
    a = list(iter_genbank_metadata((prefix + '        1 aaaaaaaaaa\n//\n').splitlines(True)))
    b = list(iter_genbank_metadata((prefix + '        1 cccccccccc\n//\n').splitlines(True)))
    assert a == b
    assert a[0].accession == 'ABC123'
    assert a[0].organism == 'Testus animalis'
    assert a[0].coi_family_annotated is True


def test_non_coi_is_retained_as_metadata_but_not_eligible_annotation():
    text = '''LOCUS       X 10 bp DNA
ACCESSION   XYZ9
  ORGANISM  Testus animalis
            Eukaryota; Metazoa; Animalia.
FEATURES             Location/Qualifiers
     source          1..10
                     /lat_lon="-12.5 130.1"
     gene            1..10
                     /gene="CYTB"
ORIGIN
        1 acgtacgtac
//
'''
    rec = list(iter_genbank_metadata(text.splitlines(True)))[0]
    assert rec.coi_family_annotated is False
    assert 'Animalia' in rec.taxonomy
