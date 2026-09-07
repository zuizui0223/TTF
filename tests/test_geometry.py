import numpy as np

from ttf import (
    SpeciesGeometry,
    SpeciesSample,
    geometry_fingerprint,
    run_geometry_calibration,
    simulate_fixed_geometry_boundary_world,
)
from ttf.core import edge_turnover
from ttf.inference import (
    heldout_species_bootstrap_from_prepared,
    heldout_species_bootstrap_test,
)
from ttf.nulls import edges_on_fixed_graphs, fixed_graphs
from ttf.transfer import prepare_transfer


def make_geometries(n_species: int = 12, n_records: int = 18):
    rng = np.random.default_rng(20260907)
    geometries = []
    for i in range(n_species):
        center = np.array([0.25 * np.cos(i), 0.25 * np.sin(i)])
        coords = center + rng.normal(0.0, 0.7, size=(n_records + i % 3, 2))
        blocks = np.asarray([f"site_{j // 3}" for j in range(len(coords))], dtype=object)
        geometries.append(
            SpeciesGeometry(
                species=f"sp_{i:02d}",
                coordinates=coords,
                blocks=blocks,
            )
        )
    return tuple(geometries)


def test_fixed_geometry_world_preserves_sampling_frame_exactly():
    geometries = make_geometries()
    assert all(not hasattr(item, "trait") for item in geometries)
    world = simulate_fixed_geometry_boundary_world(
        geometries,
        shared_fraction=0.5,
        amplitude=2.0,
        seed=17,
    )
    source = {item.species: item for item in geometries}
    for sample in world.samples:
        original = source[sample.species]
        assert np.array_equal(sample.coordinates, original.coordinates)
        assert np.array_equal(sample.blocks, original.blocks)
        assert sample.trait.shape == (len(original.coordinates),)
        assert np.isfinite(sample.trait).all()


def test_full_shared_world_uses_one_transition_hyperplane():
    geometries = make_geometries()
    world = simulate_fixed_geometry_boundary_world(
        geometries,
        shared_fraction=1.0,
        amplitude=2.0,
        seed=31,
    )
    assert len(world.shared_species) == len(geometries)
    for species in world.shared_species:
        assert np.array_equal(world.boundary_normal[species], world.shared_normal)
        assert world.boundary_offset[species] == world.shared_offset


def test_geometry_fingerprint_is_order_invariant_and_coordinate_sensitive():
    geometries = make_geometries()
    a = geometry_fingerprint(geometries)
    b = geometry_fingerprint(tuple(reversed(geometries)))
    assert a == b

    changed = list(geometries)
    coords = changed[0].coordinates.copy()
    coords[0, 0] += 1e-6
    changed[0] = SpeciesGeometry(
        species=changed[0].species,
        coordinates=coords,
        blocks=changed[0].blocks,
    )
    assert geometry_fingerprint(tuple(changed)) != a


def test_prepared_geometry_bootstrap_matches_direct_pipeline():
    geometries = make_geometries(n_species=12, n_records=16)
    world = simulate_fixed_geometry_boundary_world(
        geometries,
        shared_fraction=0.5,
        amplitude=2.0,
        noise_sd=0.5,
        transition_width=0.25,
        seed=101,
    )
    train = tuple(f"sp_{i:02d}" for i in range(6))
    evaluation = tuple(f"sp_{i:02d}" for i in range(6, 12))
    direct = heldout_species_bootstrap_test(
        world.samples,
        train_species=train,
        eval_species=evaluation,
        k=3,
        bandwidth=0.5,
        n_bootstrap=99,
        seed=17,
    )

    geometry_map = {item.species: item for item in geometries}
    used = train + evaluation
    template = [
        SpeciesSample(
            species=name,
            coordinates=geometry_map[name].coordinates,
            trait=np.arange(len(geometry_map[name].coordinates), dtype=float),
            blocks=geometry_map[name].blocks,
        )
        for name in used
    ]
    graphs = fixed_graphs(template, k=3)
    template_edges = edges_on_fixed_graphs(template, graphs)
    prepared = prepare_transfer(
        [template_edges[name] for name in train],
        [template_edges[name] for name in evaluation],
        bandwidth=0.5,
    )
    sample_map = {sample.species: sample for sample in world.samples}
    turnover = {
        name: edge_turnover(sample_map[name], graphs[name])
        for name in used
    }
    fast = heldout_species_bootstrap_from_prepared(
        prepared,
        train_turnover={name: turnover[name] for name in train},
        eval_turnover={name: turnover[name] for name in evaluation},
        n_bootstrap=99,
        seed=17,
    )
    assert np.isclose(fast.observed.statistic, direct.observed.statistic)
    assert np.allclose(fast.species_scores, direct.species_scores)
    assert np.allclose(fast.null_statistics, direct.null_statistics)
    assert np.allclose(fast.null_studentized, direct.null_studentized)
    assert fast.p_value == direct.p_value


def test_geometry_calibration_smoke_uses_heldout_species_inference():
    geometries = make_geometries(n_species=12, n_records=16)
    cells = run_geometry_calibration(
        geometries,
        shared_fractions=(0.0, 1.0),
        amplitudes=(2.0,),
        n_replicates=2,
        n_bootstrap=99,
        eval_fraction=0.5,
        k=3,
        bandwidth=0.5,
        noise_sd=0.5,
        transition_width=0.25,
        seed=91,
    )
    assert len(cells) == 2
    assert {cell.shared_fraction for cell in cells} == {0.0, 1.0}
    assert all(cell.n_replicates == 2 for cell in cells)
    assert all(0.0 <= cell.rejection_rate <= 1.0 for cell in cells)
    assert all(np.isfinite(cell.mean_statistic) for cell in cells)
    assert all(np.isfinite(cell.mean_null_statistic) for cell in cells)
