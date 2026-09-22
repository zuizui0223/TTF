import numpy as np
from ttf.relational_environment import nearest_coverage


def test_directed_coverage_is_target_fraction():
    source=np.asarray([[0.,0.,0.],[100.,0.,0.]])
    target=np.asarray([[1.,0.,0.],[99.,0.,0.],[1000.,0.,0.],[1200.,0.,0.]])
    assert np.isclose(nearest_coverage(target,source,radius=10),0.5)
    # Reversing direction changes the denominator and therefore the estimand.
    assert np.isclose(nearest_coverage(source,target,radius=10),1.0)
