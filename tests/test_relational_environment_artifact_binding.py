import importlib.util
from pathlib import Path


SCRIPT = Path("scripts/freeze_relational_environment_artifact_binding.py")


def load_module():
    spec = importlib.util.spec_from_file_location("env_binding", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_binding_schema_and_artifact_name_are_fixed():
    module = load_module()
    assert module.SCHEMA == "ttf_relational_environment_relation_artifact_binding_v0.3"
    assert module.ARTIFACT == "relational-environment-relation-v0.3"


def test_qualification_workflow_verifies_bound_relation_hashes():
    workflow = Path(".github/workflows/relational-environment-qualification-v03.yml").read_text()
    assert "relational_environment_relation_artifact_binding_v0.3.json" in workflow
    assert "files_sha256" in workflow
    assert "relational_environment_design_v0.3.json" in workflow
    assert "relational_environment_design_v0.3.npz" in workflow
    assert "occurrence_ledger_v0.2.json" in workflow
