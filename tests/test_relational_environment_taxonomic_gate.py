import importlib.util
from pathlib import Path


SCRIPT = Path("scripts/build_relational_environment_relation.py")


def load():
    spec = importlib.util.spec_from_file_location("envrelation", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def guardrail():
    return {
        "largest_single_order_fraction_max": 0.50,
        "order_fraction_threshold": 0.05,
        "minimum_orders_at_or_above_fraction_threshold": 4,
    }


def test_taxonomic_breadth_gate_passes_four_balanced_orders():
    m = load()
    names = [f"sp{i}" for i in range(20)]
    orders = ["A"] * 8 + ["B"] * 4 + ["C"] * 4 + ["D"] * 4
    metadata = {name: {"order": order} for name, order in zip(names, orders)}
    out = m.taxonomic_breadth_summary(names, metadata, guardrail())
    assert out["pass"] is True
    assert out["largest_order_fraction"] == 0.4
    assert out["orders_at_or_above_fraction_threshold_count"] == 4


def test_taxonomic_breadth_gate_fails_single_order_majority():
    m = load()
    names = [f"sp{i}" for i in range(20)]
    orders = ["A"] * 11 + ["B"] * 3 + ["C"] * 3 + ["D"] * 3
    metadata = {name: {"order": order} for name, order in zip(names, orders)}
    out = m.taxonomic_breadth_summary(names, metadata, guardrail())
    assert out["pass"] is False
    assert out["largest_order_fraction"] == 0.55


def test_taxonomic_breadth_gate_fails_fewer_than_four_five_percent_orders():
    m = load()
    names = [f"sp{i}" for i in range(100)]
    orders = ["A"] * 48 + ["B"] * 47 + ["C"] * 5
    metadata = {name: {"order": order} for name, order in zip(names, orders)}
    out = m.taxonomic_breadth_summary(names, metadata, guardrail())
    assert out["pass"] is False
    assert out["largest_order_fraction"] == 0.48
    assert out["orders_at_or_above_fraction_threshold_count"] == 3
