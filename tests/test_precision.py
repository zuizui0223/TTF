from ttf.calibration import CalibrationCell
from ttf.precision import qualify_calibration_precision, wilson_interval


def cell(shared: float, amplitude: float, n: int, rejected: int) -> CalibrationCell:
    rate = rejected / n
    return CalibrationCell(shared, amplitude, n, 0.05, rate, 0.0, 0.0, 0.5)


def test_wilson_interval_contains_estimate():
    interval = wilson_interval(35, 500)
    assert interval.low < interval.estimate < interval.high
    assert interval.estimate == 0.07
    assert interval.high < 0.10


def test_precision_gate_requires_zero_shared_upper_bound_and_power_lower_bound():
    cells = (
        cell(0.0, 0.5, 500, 24),
        cell(0.0, 1.0, 500, 28),
        cell(0.0, 2.0, 500, 25),
        cell(0.0, 3.0, 500, 30),
        cell(1.0, 2.0, 500, 475),
    )
    report = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    )
    assert report.passed
    assert report.type1_pass
    assert report.power_pass
    assert report.max_zero_shared_upper95 < 0.10
    assert report.full_shared_moderate_lower95 > 0.80


def test_precision_gate_fails_even_when_point_type1_is_below_ceiling_if_ci_is_not():
    cells = (
        cell(0.0, 1.0, 100, 7),
        cell(1.0, 2.0, 100, 96),
    )
    report = qualify_calibration_precision(
        cells,
        moderate_amplitude=2.0,
        type1_upper_ceiling=0.10,
        power_lower_floor=0.80,
    )
    assert report.max_zero_shared_estimate == 0.07
    assert report.max_zero_shared_upper95 > 0.10
    assert not report.type1_pass
    assert not report.passed
