from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_ttf_q_theory_contains_no_accidental_ascii_control_characters():
    text=(ROOT/"docs/TTF_Q_THEORY_V01.md").read_text(encoding="utf-8")
    forbidden={
        "\x08":"backspace from an unescaped \\b sequence",
        "\x0b":"vertical tab from an unescaped \\v sequence",
        "\x0c":"form feed from an unescaped \\f sequence",
    }
    for char,reason in forbidden.items():
        assert char not in text,reason

    for token in (
        r"\mathbb{R}",
        r"\frac",
        r"\beta",
        r"\tau",
        r"\mathcal{E}",
        r"\operatorname{col}",
        r"\begin{cases}",
    ):
        assert token in text
