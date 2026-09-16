import importlib.util
from pathlib import Path

_SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'census_direct_insdc_metadata.py'
_SPEC = importlib.util.spec_from_file_location('census_direct_insdc_metadata', _SCRIPT)
_MOD = importlib.util.module_from_spec(_SPEC)
assert _SPEC and _SPEC.loader
_SPEC.loader.exec_module(_MOD)
parse_latlon = _MOD.parse_latlon


def test_latlon_parser_common_genbank_forms():
    assert parse_latlon('35.0 N 135.0 E') == (35.0, 135.0)
    assert parse_latlon('12.5 S 130.1 E') == (-12.5, 130.1)
    assert parse_latlon('-12.5 130.1') == (-12.5, 130.1)
    assert parse_latlon('95 N 10 E') is None
