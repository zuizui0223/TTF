from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

import pytest

from ttf.genetic_manuscript_export import (
    PREFIX, OPEN_FLAGS, build_closed_export, build_empirical_export, write_export,
)


def _write(path, payload):
    path.write_text(json.dumps(payload, sort_keys=True) + '\n')
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _chain(tmp_path, primary_p=0.05, self_p=0.01, qualified=True, unavailable=False):
    """Entirely fabricated reporting fixtures; these are not qualification evidence."""
    paths = {key: tmp_path / f'{key}.json' for key in (
        'result', 'authorization', 'qualification', 'self_qualification', 'phase4_rule', 'self_rule')}
    closed = {key: False for key in OPEN_FLAGS}
    fingerprint = 'a' * 64
    train = [f'train_{i:03d}' for i in range(75)]
    evaluation = [f'eval_{i:03d}' for i in range(75)]
    self_status = 'PASS' if qualified else 'SELF_DETECTABILITY_NOT_QUALIFIED'
    hashes = {}
    hashes['qualification'] = _write(paths['qualification'], {
        'schema': PREFIX + 'phase3_qualification_v0.1', 'status': 'PASS', 'passed': True,
        'phase4_identity_opening_eligible': True, 'geometry_fingerprint_sha256': fingerprint, **closed,
        'type1_gate': {'pass': True, 'wilson95_upper_ceiling': 0.1, 'max_observed_wilson95_upper': 0.08},
        'power_gate': {'pass': True, 'wilson95_lower_floor': 0.8, 'observed_wilson95_lower': 0.9},
    })
    hashes['self_qualification'] = _write(paths['self_qualification'], {
        'schema': PREFIX + 'phase3_self_qualification_v0.1', 'status': self_status,
        'passed': qualified, 'geometry_fingerprint_sha256': fingerprint, **closed,
    })
    for key, source in [
        ('phase4_rule', 'docs/supporting/genetic_phylogatr_phase4_response_rule_v0.1.json'),
        ('self_rule', 'docs/supporting/genetic_phylogatr_phase3_self_detectability_rule_v0.1.json')]:
        paths[key].write_bytes(Path(source).read_bytes())
        hashes[key] = hashlib.sha256(paths[key].read_bytes()).hexdigest()
    auth = {
        'schema': PREFIX + 'phase4_identity_opening_authorization_v0.1',
        'status': 'AUTHORIZE_EXACT_FRESH_NUCLEOTIDE_IDENTITY_OPENING',
        'geometry_fingerprint_sha256': fingerprint,
        'species': {'train_species': train, 'eval_species': evaluation},
        'phase3_gate_d': {'status': 'PASS', 'passed': True, 'type1_gate_pass': True, 'power_gate_pass': True},
        'fresh_self_detectability': {'status': self_status, 'passed': qualified},
        'outcome_firewall_before_execution': {**closed, 'decker_empirical_genetic_outcomes_opened': False},
        'phase3_qualification_sha256': hashes['qualification'],
        'phase3_self_qualification_sha256': hashes['self_qualification'],
        'phase4_rule_sha256': hashes['phase4_rule'], 'phase3_self_rule_sha256': hashes['self_rule'],
    }
    auth_hash = _write(paths['authorization'], auth)
    decision = ('TRANSFERABLE_PLACE_COMPONENT' if primary_p <= 0.05 else
                'LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN' if qualified and self_p <= 0.05 else
                'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING')
    primary = dict(statistic=0.2, profiled_private_p_value=primary_p, alpha=0.05,
                   positive=primary_p <= 0.05, selected_configurations=['A0', 'A1'],
                   component_p_values={'A0': primary_p, 'A1': primary_p / 2},
                   heldout_species_scores={name: 0.2 for name in evaluation})
    self_result = dict(statistic=0.4, p_value=self_p, alpha=0.05, positive=self_p <= 0.05,
                       synthetic_method_qualified=qualified, synthetic_qualification_status=self_status,
                       species_scores={name: 0.4 for name in evaluation})
    total_scores = {name: 0.9 for name in evaluation}
    if unavailable:
        total_scores[evaluation[0]] = None
    total = dict(name='total_genetic_distance_transfer',
                 status='NOT_EVALUABLE_DESCRIPTIVE' if unavailable else 'DESCRIPTIVE_ONLY',
                 statistic=None if unavailable else 0.9, heldout_species_scores=total_scores,
                 inferentially_qualified=False, used_for_primary_decision=False,
                 biological_ibd_adjustment=False, same_primary_geometry_split_kernel=True,
                 training_geometry_control='within_species_edge_length_rank_orthogonalization',
                 undefined_species=[evaluation[0]] if unavailable else [],
                 frozen_evaluation_species_count=len(evaluation),
                 finite_evaluation_species_count=len(evaluation) - int(unavailable))
    result = {
        'schema': PREFIX + 'phase4_empirical_result_v0.1',
        'status': 'EMPIRICAL_RESULT_OPENED_UNDER_FROZEN_PHASE4_AUTHORIZATION',
        'phase4_authorization_sha256': auth_hash, 'geometry_fingerprint_sha256': fingerprint,
        **{key: True for key in OPEN_FLAGS}, 'decker_empirical_genetic_outcomes_opened': False,
        'serialized_sequence_identity': False, 'serialized_edge_genetic_distance_vectors': False,
        'primary_place_beyond_ibd': primary, 'within_species_self_diagnostic': self_result,
        'secondary_total_genetic_transfer': total, 'decision': decision,
        'source_integrity': {'cite_sha256_verified': True, 'genes_sha256_verified': True,
            'all_survivor_headers_occurrences_masks_verified': True,
            'species': {name: {'header_sha_verified': True, 'occurrence_sha_verified': True,
                'phase2_mask_sha_verified': True} for name in train + evaluation}},
        'interpretation': 'UNTRUSTED PROSE: a universal causal barrier has been proven',
    }
    _write(paths['result'], result)
    return paths


@pytest.mark.parametrize('p,sp,qualified,decision', [
    (0.05, 0.01, True, 'TRANSFERABLE_PLACE_COMPONENT'),
    (0.050001, 0.05, True, 'LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN'),
    (0.2, 0.01, False, 'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING'),
    (0.2, 0.2, True, 'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING'),
])
@pytest.mark.parametrize('unavailable', [False, True])
def test_export_all_branches_and_full_species_table(tmp_path, p, sp, qualified, decision, unavailable):
    paths = _chain(tmp_path, p, sp, qualified, unavailable)
    report = build_empirical_export(**paths)
    assert report.manifest['decision'] == decision
    rows = list(csv.DictReader(io.StringIO(report.species_csv)))
    assert len(rows) == 75 and rows[0]['species'] == 'eval_000'
    assert (rows[0]['total_transfer_descriptive'] == '') is unavailable
    assert 'UNTRUSTED PROSE' not in report.markdown
    assert 'no significance test' in report.markdown
    if decision.startswith('LINEAGE'):
        assert 'consistent with' in report.markdown and 'does not establish zero transfer' in report.markdown
    before = {key: path.read_bytes() for key, path in paths.items()}
    output = write_export(report, tmp_path / 'report')
    assert {p.name for p in output.iterdir()} == {'results.md', 'species_scores.csv', 'export_manifest.json'}
    manifest = json.loads((output / 'export_manifest.json').read_text())
    for name, digest in manifest['output_sha256'].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    assert before == {key: path.read_bytes() for key, path in paths.items()}
    with pytest.raises(ValueError, match='already exists'):
        write_export(report, output)


@pytest.mark.parametrize('mutation', ['unopened', 'hash', 'geometry', 'species', 'mean', 'decision',
                                     'p_max', 'alpha', 'self_support', 'secondary_p', 'secondary_mean', 'secondary_count', 'nan', 'bool_score', 'source_integrity', 'secondary_control'])
def test_export_rejects_inconsistent_result(tmp_path, mutation):
    paths = _chain(tmp_path)
    r = json.loads(paths['result'].read_text())
    p, s, t = r['primary_place_beyond_ibd'], r['within_species_self_diagnostic'], r['secondary_total_genetic_transfer']
    if mutation == 'unopened': r[OPEN_FLAGS[0]] = False
    elif mutation == 'hash': r['phase4_authorization_sha256'] = 'b' * 64
    elif mutation == 'geometry': r['geometry_fingerprint_sha256'] = 'b' * 64
    elif mutation == 'species': del p['heldout_species_scores']['eval_000']
    elif mutation == 'mean': p['statistic'] = 0.3
    elif mutation == 'decision': r['decision'] = 'NOT_EVALUABLE_FOR_LINEAGE_CONDITIONING'
    elif mutation == 'p_max': p['component_p_values']['A1'] = 0.1
    elif mutation == 'alpha': p['alpha'] = 0.1
    elif mutation == 'self_support': s['synthetic_method_qualified'] = False
    elif mutation == 'secondary_p': t['p_value'] = 0.01
    elif mutation == 'secondary_mean': t['heldout_species_scores']['eval_000'] = None
    elif mutation == 'secondary_count': t['finite_evaluation_species_count'] = 74
    elif mutation == 'nan': p['statistic'] = float('nan')
    elif mutation == 'source_integrity': r['source_integrity']['species']['eval_000']['phase2_mask_sha_verified'] = False
    elif mutation == 'secondary_control': t['biological_ibd_adjustment'] = True
    else: p['heldout_species_scores']['eval_000'] = True
    _write(paths['result'], r)
    with pytest.raises(ValueError): build_empirical_export(**paths)
    assert not (tmp_path / 'report').exists()


def test_export_rejects_changed_receipt_even_if_geometry_matches(tmp_path):
    paths = _chain(tmp_path)
    q = json.loads(paths['qualification'].read_text())
    q['status'] = 'NOT_EVALUABLE'
    _write(paths['qualification'], q)
    with pytest.raises(ValueError, match='hash mismatch'): build_empirical_export(**paths)


def test_export_rejects_failed_gate_even_with_matching_hash_chain(tmp_path):
    paths = _chain(tmp_path)
    q = json.loads(paths['qualification'].read_text())
    q['type1_gate']['max_observed_wilson95_upper'] = 0.101
    digest = _write(paths['qualification'], q)
    a = json.loads(paths['authorization'].read_text())
    a['phase3_qualification_sha256'] = digest
    digest = _write(paths['authorization'], a)
    r = json.loads(paths['result'].read_text())
    r['phase4_authorization_sha256'] = digest
    _write(paths['result'], r)
    with pytest.raises(ValueError, match='Type-I gate inconsistent'): build_empirical_export(**paths)


def test_canonical_closed_receipt_produces_only_status(tmp_path):
    report = build_closed_export(Path('benchmarks/frozen/genetic_empirical_opening_state_v0.3.json'))
    output = write_export(report, tmp_path / 'closed')
    assert {p.name for p in output.iterdir()} == {'status.md', 'export_manifest.json'}
    assert report.manifest['empirical_results_generated'] is False
    assert report.manifest['decision'] == 'AWAITING_FRESH_ARCHIVE'


@pytest.mark.parametrize('schema,status', [
    ('confirmatory_phase1_geometry_v0.1', 'NOT_EVALUABLE_PHASE1_PANEL_TOO_SMALL'),
    ('confirmatory_phase2_mask_v0.1', 'NOT_EVALUABLE_PHASE2_PANEL_TOO_SMALL'),
    ('phase3_qualification_v0.1', 'NOT_EVALUABLE'),
])
def test_failed_phases_do_not_become_biological_nulls(tmp_path, schema, status):
    path = tmp_path / 'closed.json'
    payload = {'schema': PREFIX + schema, 'status': status, **{key: False for key in OPEN_FLAGS},
               'passed': False, 'phase4_identity_opening_eligible': False}
    _write(path, payload)
    report = build_closed_export(path)
    assert report.species_csv == ''
    assert 'does not establish absence' in report.markdown
    payload[OPEN_FLAGS[0]] = True
    _write(path, payload)
    with pytest.raises(ValueError): build_closed_export(path)



def test_negative_self_statistic_is_reported_null_relative(tmp_path):
    paths = _chain(tmp_path, primary_p=0.2, self_p=0.01, qualified=True)
    payload = json.loads(paths['result'].read_text())
    payload['within_species_self_diagnostic']['statistic'] = -0.25
    for name in payload['within_species_self_diagnostic']['species_scores']:
        payload['within_species_self_diagnostic']['species_scores'][name] = -0.25
    _write(paths['result'], payload)
    report = build_empirical_export(**paths)
    assert 'structural-null reference' in report.markdown
    assert 'raw numerical sign is not interpreted against zero' in report.markdown
    assert report.manifest['decision'] == 'LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN'

def test_duplicate_json_keys_rejected(tmp_path):
    path = tmp_path / 'bad.json'
    path.write_text('{"status":"PASS","status":"NOT_EVALUABLE"}')
    with pytest.raises(ValueError, match='duplicate JSON'): build_closed_export(path)


def test_cli_closed_mode_and_incomplete_empirical_mode(tmp_path):
    cmd = [sys.executable, 'scripts/export_genetic_ttf_manuscript.py']
    done = subprocess.run(cmd + ['--closed-receipt', 'benchmarks/frozen/genetic_empirical_opening_state_v0.3.json',
                                 '--output-dir', str(tmp_path / 'closed')], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    paths = _chain(tmp_path)
    failed = subprocess.run(cmd + ['--result', str(paths['result']), '--output-dir', str(tmp_path / 'empirical')], capture_output=True, text=True)
    assert failed.returncode != 0 and not (tmp_path / 'empirical').exists()
    args = [value for key, path in paths.items() for value in ('--' + key.replace('_', '-'), str(path))]
    done = subprocess.run(cmd + args + ['--output-dir', str(tmp_path / 'empirical')], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert (tmp_path / 'empirical' / 'species_scores.csv').is_file()
