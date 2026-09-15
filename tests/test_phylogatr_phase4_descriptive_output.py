from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ttf.genetic_empirical_score import DescriptiveTotalTransfer
from test_genetic_empirical_score import _geometries, _design


def _runner():
    path = Path('scripts/run_phylogatr_phase4_empirical_test.py')
    spec = importlib.util.spec_from_file_location('phase4_descriptive_runner', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _argv(root, output):
    flags = ['geometry', 'phase1-manifest', 'phase2-manifest', 'phase3-rule',
             'phase3-authorization', 'references', 'qualification', 'self-rule',
             'self-references', 'self-qualification', 'phase4-rule',
             'phase4-authorization', 'opening-state']
    return ['runner', '--root', str(root), '--output', str(output)] + [
        value for flag in flags for value in ('--' + flag, 'synthetic-fixture.json')]


@pytest.mark.parametrize('primary_positive,self_positive,self_qualified,expected', [
    (True, False, False, 'TRANSFERABLE_PLACE_COMPONENT'),
    (False, True, True, 'LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN'),
    (False, True, False, 'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING'),
    (False, False, True, 'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING'),
])
@pytest.mark.parametrize('total_value', [None, -0.9, 0.9])
def test_runner_decision_ignores_total_descriptor(
    tmp_path, monkeypatch, primary_positive, self_positive, self_qualified, expected, total_value
):
    """Integration of result assembly; authorization itself has separate real-chain tests.

    All source/extraction seams below use constructed arrays, never empirical data.
    Keep the real primary/self numerical adapters and profiled-private inference.
    """
    module = _runner()
    geometries = _geometries()
    names, design = _design(geometries)
    rng = np.random.default_rng(418)
    distances = {name: rng.uniform(0, 1, len(g.edge_nodes)) for name, g in geometries.items()}
    (tmp_path / 'dummy.fasta').write_text('synthetic stub, no sequence')
    (tmp_path / 'dummy.csv').write_text('synthetic stub')
    context = SimpleNamespace(
        train_species=design.train_species, eval_species=design.eval_species,
        geometries=geometries, frozen_latlon={name: name for name in names},
        phase1={'provenance': {'cite_sha256': 'fixture', 'genes_sha256': 'fixture'},
                'selected_panels': {name: {'aligned_fasta_relative': 'dummy.fasta',
                    'occurrence_relative': 'dummy.csv', 'occurrence_sha256': 'fixture',
                    'aligned_header_sha256': 'fixture'} for name in names}},
        phase2={'species_ledger': {name: {'survives': True, 'canonical_mask_sha256': 'fixture'} for name in names}},
        phase3_rule={'core_method': {'bandwidth_km': 500, 'prior_strength': 0.25,
                    'segment_points': 5, 'strength_neighbours': 4,
                    'profiled_private_selected_configurations': 2},
                    'geometry_contract': {'minimum_endpoint_disjoint_ibd_training_edges': 5},
                    'qualification': {'profile_strength_draws': 20}},
        phase3_authorization={'execution': {'edge_chunk_size': 32, 'train_chunk_size': 4096},
                              'geometry_fingerprint_sha256': 'synthetic'},
        references={'references': {label: {'training_strength': [0.0] * 120,
                    'statistic': [-1.0 if primary_positive else 1.0] * 120} for label in ['A', 'B']}},
        self_rule={'geometry': {'bandwidth_km': 500, 'prior_strength': 0.25,
                               'prior_mean': 0.5, 'segment_points': 5}, 'inference': {'alpha': 0.05}},
        self_references={'statistics': [-1.0 if self_positive else 1.0] * 120},
        self_qualification={'passed': self_qualified, 'status': 'PASS' if self_qualified else 'SELF_DETECTABILITY_NOT_QUALIFIED'},
        phase4_rule={'sequence_contract': {'minimum_comparable_fraction_of_frozen_alignment_length': 0.5},
                     'primary_estimand': {'alpha': 0.05}},
    )
    monkeypatch.setattr(module, 'load_phylogatr_phase4_context', lambda *a, **kw: context)
    monkeypatch.setattr(module, 'sha256_path', lambda path: 'fixture')
    monkeypatch.setattr(module, 'fasta_headers_only', lambda path: ())
    monkeypatch.setattr(module, 'fasta_header_sha256', lambda headers: 'fixture')
    monkeypatch.setattr(module, 'read_canonical_mask_alignment', lambda path: None)
    monkeypatch.setattr(module, 'canonical_mask_sha256', lambda mask: 'fixture')
    monkeypatch.setattr(module, 'extract_species_frozen_edge_distances', lambda f, o, name, g, **kw:
                        SimpleNamespace(genetic_distance=distances[name], valid_pair_counts=np.ones(len(distances[name]), dtype=int)))
    monkeypatch.setattr(module, 'score_total_genetic_distance_mapping', lambda *a, **kw:
                        DescriptiveTotalTransfer(total_value, {name: total_value for name in design.eval_species}))
    output = tmp_path / 'result.json'
    monkeypatch.setattr(sys, 'argv', _argv(tmp_path, output))
    assert module.main() == 0
    payload = json.loads(output.read_text())
    assert payload['decision'] == expected
    assert payload['primary_place_beyond_ibd']['positive'] is primary_positive
    secondary = payload['secondary_total_genetic_transfer']
    assert secondary['statistic'] == total_value
    assert secondary['used_for_primary_decision'] is False
    assert 'p_value' not in secondary and 'positive' not in secondary
    assert set(secondary['heldout_species_scores']) == set(design.eval_species)
    assert payload['serialized_edge_genetic_distance_vectors'] is False
    assert payload['serialized_sequence_identity'] is False
    assert 'NaN' not in output.read_text()


def test_runner_does_not_compute_descriptor_before_authorization(tmp_path, monkeypatch):
    module = _runner()
    def denied(*args, **kwargs):
        raise RuntimeError('Phase-3 Gate-D did not PASS')
    def forbidden(*args, **kwargs):
        pytest.fail('attempted source opening or descriptive scoring without authorization')
    monkeypatch.setattr(module, 'load_phylogatr_phase4_context', denied)
    monkeypatch.setattr(module, 'extract_species_frozen_edge_distances', forbidden)
    monkeypatch.setattr(module, 'score_total_genetic_distance_mapping', forbidden)
    output = tmp_path / 'result.json'
    monkeypatch.setattr(sys, 'argv', _argv(tmp_path, output))
    with pytest.raises(RuntimeError, match='did not PASS'):
        module.main()
    assert not output.exists()
