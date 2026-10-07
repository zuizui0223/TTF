from scripts.census_historical_host_neotoma import api_dataset_type, norm_name, parse_records


def test_api_dataset_type_normalizes_frozen_label():
    assert api_dataset_type("plant macrofossils") == "plant macrofossil"
    assert api_dataset_type("pollen") == "pollen"


def test_parse_occurrence_payload_extracts_taxon_and_site():
    payload = {
        "data": [
            {
                "ages": {"age": 22000, "ageolder": 22500, "ageyounger": 21500},
                "sample": {"taxonname": "Quercus alba"},
                "site": {"siteid": 42, "datasettype": "plant macrofossil"},
                "taxon": {"taxonname": "Quercus alba"},
            },
            {
                "sample": {"taxonname": ""},
                "site": {"siteid": 99},
                "taxon": {"taxonname": ""},
            },
        ]
    }
    got = parse_records(payload)
    assert got == [
        {
            "taxonname": "Quercus alba",
            "siteid": 42,
            "datasettype": "plant macrofossil",
            "age": 22000,
            "ageolder": 22500,
            "ageyounger": 21500,
        }
    ]


def test_norm_name_collapses_whitespace():
    assert norm_name("  Pinus   banksiana ") == "Pinus banksiana"
