import numpy as np

from ttf.climate_occurrence_filter import climate_validity_mask


def test_climate_validity_removes_nonfinite_and_nodata_cells():
    x = np.asarray([
        [1.0, 2.0, 3.0],
        [np.nan, 2.0, 3.0],
        [1.0, -9999.0, 3.0],
        [1.0, 2.0, np.inf],
        [4.0, 5.0, 6.0],
    ])
    keep = climate_validity_mask(
        x,
        nodata_values=[None, -9999.0, None],
    )
    assert keep.tolist() == [True, False, False, False, True]
