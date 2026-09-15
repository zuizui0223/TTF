from scripts.census_direct_insdc_metadata import parse_latlon


def test_latlon_parser_common_genbank_forms():
    assert parse_latlon('35.0 N 135.0 E') == (35.0, 135.0)
    assert parse_latlon('12.5 S 130.1 E') == (-12.5, 130.1)
    assert parse_latlon('-12.5 130.1') == (-12.5, 130.1)
    assert parse_latlon('95 N 10 E') is None
