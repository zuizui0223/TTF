from pathlib import Path

def test_projected_phase4_wrapper_reconstructs_exact_phase1_source_view() -> None:
    source = Path("scripts/run_phylogatr_projected_phase4_empirical_test.py").read_text()

    for required in (
        "project_genes_file_in_place",
        "raw_archive_sha256",
        "phase1_manifest_sha256",
        "projection_rule_sha256",
        "archive_root_relative",
        "_require_same_projection",
        "run_phylogatr_phase4_empirical_test.py",
        "--root",
    ):
        assert required in source

    assert "sequence_identity_opened" in source
    assert "pairwise_genetic_distances_opened" in source
    assert "temporary_projection_view_deleted_after_phase4" in source


def test_phase4_authorizer_freezes_projected_archive_execution_surface() -> None:
    source = Path("scripts/authorize_phylogatr_phase4_identity_opening.py").read_text()
    for required in (
        "src/ttf/phylogatr_source_projection.py",
        "docs/supporting/genetic_phylogatr_source_projection_rule_v0.1.json",
        "scripts/run_phylogatr_projected_phase4_empirical_test.py",
    ):
        assert required in source
