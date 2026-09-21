from ttf.lepidoptera_host_resource_qualification import REFERENCE_TAG as V01_REFERENCE_TAG
from ttf.lepidoptera_host_resource_qualification_v02 import (
    EVALUATION_TAG,
    REFERENCE_TAG,
    frozen_v02_seed,
)


def test_v02_seed_namespace_is_disjoint_from_v01():
    assert REFERENCE_TAG=="lepidoptera-host-resource-v02-reference"
    assert EVALUATION_TAG=="lepidoptera-host-resource-v02-evaluation"
    assert REFERENCE_TAG!=V01_REFERENCE_TAG
    a=frozen_v02_seed(20260921,REFERENCE_TAG,"private",0)
    b=frozen_v02_seed(20260921,EVALUATION_TAG,"private",0)
    assert a!=b
