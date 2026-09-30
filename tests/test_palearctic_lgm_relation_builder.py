import numpy as np

from ttf.lgm_climate_envelope import fit_common_climate_space


def test_lgm_design_common_space_retains_all_four_frozen_climate_axes():
    rng=np.random.default_rng(8)
    x=rng.normal(size=(100,4))
    space=fit_common_climate_space(x,axes=4)
    assert space.axes==4
    assert space.eigvec.shape==(4,4)
    assert len(space.eigval)==4
