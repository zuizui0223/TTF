"""Deterministic design anatomy and idealized NEW paired-block planning.

Inputs are covariates and group IDs only. No phenotype outcome interface,
no simulation, no coefficient fitting, and no authorization to reopen v1.
All reported spans are fixed in advance, not selected for a positive result.
"""
from __future__ import annotations
from collections import Counter
import numpy as np
from scipy.stats import nct, t

ALPHA = 1 / 60
RANK_RTOL = 1e-10


def _labels(values, n):
    a = np.asarray(values)
    if a.shape != (n,) or any(v is None or not str(v).strip() for v in a):
        raise ValueError('Nonempty aligned source group labels required')
    return a.astype(str)


def _incidence(labels):
    return (labels[:, None] == np.unique(labels)[None, :]).astype(float)


def span_summary(matrix, history):
    """Source predictor residual energy; it is not biological variance explained."""
    A = np.asarray(matrix, float)
    h = np.asarray(history, float)
    if (A.ndim != 2 or len(A) < 3 or h.shape != (len(A),) or
            not np.isfinite(A).all() or not np.isfinite(h).all()):
        raise ValueError('Finite aligned source design and history required')
    norms = np.linalg.norm(A, axis=0)
    use = norms > 0
    if np.any(use):
        U, s, _ = np.linalg.svd(A[:, use] / norms[use], full_matrices=False)
        Q = U[:, s > RANK_RTOL * s[0]]
    else:
        Q = np.zeros((len(h), 0))
    residual = h - Q @ (Q.T @ h)
    energy = float(residual @ residual)
    identified = energy > RANK_RTOL**2 * max(1., float(h @ h))
    return {'rank': Q.shape[1], 'history_energy': energy,
            'history_identifiable': bool(identified),
            'residual_df_after_history': len(h) - Q.shape[1] - int(identified)}


def _ratio(numerator, denominator):
    return numerator / denominator if denominator > 1e-20 else None


def analyse_geometry(X, history, species, studies):
    """All ten declared nuisance spans and same-study/same-species cell support."""
    X = np.asarray(X, float)
    h = np.asarray(history, float)
    if (X.ndim != 2 or X.shape[1] != 21 or len(X) < 3 or
            h.shape != (len(X),) or not np.isfinite(X).all() or not np.isfinite(h).all()):
        raise ValueError('Exactly 21 finite source columns and aligned history required')
    n = len(h)
    sp, st = _labels(species, n), _labels(studies, n)
    cells = list(zip(sp.tolist(), st.tolist()))
    ordered = sorted(set(cells))
    cell_id = np.array([ordered.index(k) for k in cells])
    S, T, C = _incidence(sp), _incidence(st), _incidence(cell_id)
    one = np.ones((n, 1))
    G = np.column_stack([S, T])
    matrices = {
        'intercept': one,
        'environment11': np.column_stack([one, X[:, :11]]),
        'baseline21': np.column_stack([one, X]),
        'species': S, 'study': T, 'species_study': G,
        'groups_environment11': np.column_stack([G, one, X[:, :11]]),
        'groups_baseline21': np.column_stack([G, one, X]),
        'study_species_cells': C,
        'cells_baseline21': np.column_stack([C, one, X]),
    }
    spans = {name: span_summary(A, h) for name, A in matrices.items()}
    counts = Counter(cells)
    repeated = [key for key, size in counts.items() if size >= 2]
    varying = []
    for k, key in enumerate(ordered):
        hc = h[cell_id == k]
        if len(hc) > 1 and float(((hc - hc.mean())**2).sum()) > RANK_RTOL**2 * max(1., float(hc @ hc)):
            varying.append(key)
    g, e, f = (spans[k] for k in ('species_study', 'groups_environment11', 'groups_baseline21'))
    for a, b in ((g, e), (e, f)):
        if b['rank'] < a['rank'] or b['history_energy'] > a['history_energy'] + 1e-8 * max(1., a['history_energy']):
            raise ValueError('Nested source span numerical inconsistency')
    result = {
        'n': n, 'species': len(np.unique(sp)), 'studies': len(np.unique(st)),
        'spans': spans,
        'cell_support': {
            'occupied_cells': len(counts),
            'singleton_cells': sum(size == 1 for size in counts.values()),
            'repeated_cells': len(repeated),
            'rows_in_repeated_cells': sum(counts[k] for k in repeated),
            'species_with_repeated_cell': len({s for s, _ in repeated}),
            'studies_with_repeated_cell': len({r for _, r in repeated}),
            'cells_with_history_variation': len(varying),
            'species_with_within_cell_history_variation': len({s for s, _ in varying}),
            'studies_with_within_cell_history_variation': len({r for _, r in varying}),
            'rows_per_cell_histogram': {str(size): number for size, number in sorted(Counter(counts.values()).items())},
        },
        'rank_budget': {
            'within_group_dimensions_before_baseline': n - g['rank'],
            'environment11_additional_rank_after_groups': e['rank'] - g['rank'],
            'method10_additional_rank_after_environment_and_groups': f['rank'] - e['rank'],
            'history_additional_rank': int(f['history_identifiable']),
            'residual_df_after_history': f['residual_df_after_history'],
        },
        'history_energy_ratios': {
            'group_only_over_intercept': _ratio(g['history_energy'], spans['intercept']['history_energy']),
            'environment_given_groups_over_groups': _ratio(e['history_energy'], g['history_energy']),
            'method_given_groups_environment_over_groups_environment': _ratio(f['history_energy'], e['history_energy']),
            'groups_given_baseline21_over_baseline21': _ratio(f['history_energy'], spans['baseline21']['history_energy']),
            'all_controls_given_groups_over_groups': _ratio(f['history_energy'], g['history_energy']),
        },
        'decomposition_order_is_descriptive_not_causal': True,
        'all_spans_reported_not_selected_as_winning_models': True,
        'biological_effect_tested': False,
        'empirical_opening_authorized': False,
    }
    return result


def paired_power(blocks, delta, gamma=.5, sigma=.5):
    """Ideal two-observation blocks, equal nuisance within block, iid normal error.

Each block contains h=(-delta/2,+delta/2), a separate intercept, and a common
history slope. Energy=blocks*delta**2/2; unknown-variance t df=blocks-1.
"""
    if (not all(np.isfinite(v) for v in (blocks, delta, gamma, sigma)) or
            blocks < 2 or int(blocks) != blocks or delta <= 0 or sigma <= 0):
        raise ValueError('Integer blocks >=2, positive history contrast/noise required')
    df = int(blocks) - 1
    cutoff = t.isf(ALPHA / 2, df)
    nc = gamma * np.sqrt(blocks * delta**2 / 2) / sigma
    return float(nct.cdf(-cutoff, df, nc) + nct.sf(cutoff, df, nc))


def required_blocks(delta, target=.8):
    if not np.isfinite(target) or not ALPHA < target < 1:
        raise ValueError('Power target must exceed alpha and be less than one')
    for blocks in range(2, 5001):
        power = paired_power(blocks, delta)
        if power >= target:
            return {'delta_history_in_original_training_SD': delta,
                    'blocks': blocks, 'rate_comparison_records': 2 * blocks,
                    'energy': blocks * delta**2 / 2, 'residual_df': blocks - 1,
                    'power': power, 'gamma': .5, 'sigma': .5,
                    'new_independent_perfectly_matched_blocks_assumed': True,
                    'not_a_guaranteed_empirical_sample_size': True}
    raise ValueError('Required ideal blocks exceed bounded planning range')
