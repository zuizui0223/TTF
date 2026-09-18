from pathlib import Path


MANUSCRIPT = Path("manuscript/genetic_ttf_flagship_v0.2.md")
GOAL = Path("docs/GENETIC_TTF_DEVELOPMENT_GOAL.md")


def test_genetic_empirical_manuscript_current_state_is_not_predatified() -> None:
    manuscript = MANUSCRIPT.read_text()
    goal = GOAL.read_text()

    assert "LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN" in manuscript
    assert "T = 0.0325655" in manuscript
    assert "upper-tail `p = 0.008991`" in manuscript
    assert "genetic_phase4_v0.1/figure1_qualification_margins.svg" in manuscript
    assert "genetic_phase4_v0.1/figure2_post_ibd_species_scores.svg" in manuscript

    assert "10.1016/j.ympev.2009.09.016" in manuscript
    assert "a future Phase-4 authorization binds" not in goal
    assert "the completed empirical opening" in goal
    assert "the one-shot empirical result are now frozen" in goal
