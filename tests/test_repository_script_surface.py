from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path("scripts/audit_repository_script_surface.py")


def load_module():
    spec = importlib.util.spec_from_file_location("surface_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_repository_script_surface_has_no_missing_references():
    module = load_module()
    result = module.build_surface(Path(".").resolve())
    assert result["schema"] == "ttf_repository_script_surface_audit_v0.1"
    assert result["missing_script_references"] == []
    assert result["script_count"] >= result["active_referenced_count"]
    active = {row["path"] for row in result["active_referenced"]}
    assert "scripts/run_phylogatr_phase4_empirical_test.py" in active
    assert "scripts/run_relational_environment_qualification.py" in active
    assert "scripts/run_v11_fresh_density_scaled_statistics.py" in active


def test_repository_script_surface_is_deterministic():
    module = load_module()
    first = module.build_surface(Path(".").resolve())
    second = module.build_surface(Path(".").resolve())
    assert first == second
