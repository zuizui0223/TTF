import numpy as np

from ttf import (
    SpeciesGeometry,
    geometry_fingerprint,
    run_geometry_calibration,
    simulate_fixed_geometry_boundary_world,
)


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
