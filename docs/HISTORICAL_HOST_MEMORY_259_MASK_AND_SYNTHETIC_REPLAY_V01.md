# Historical-host memory — confirmatory 259-species mask and post-mask qualification

**Status as of 2026-10-11:** `PASS_TO_SEPARATE_CONFIRMATORY_IDENTITY_AUTHORIZATION_REVIEW`. This is a **source-checked pre-genetic qualification result**, not a positive historical-host genetic effect.

## What was independently replayed

The full original phylogatR archive resides **privately** in Library at
`/TTF/source/phylogatr_results_exact_20260918.zip` (274,988,692 bytes,
SHA256 `5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5`).
It was not uploaded to public GitHub.

On the exact original **321 confirmatory species**, independently re-running
the Boolean A/C/G/T validity mask against the frozen locality and edge geometry
produced **259 eligible species**, **61 failed because at least one unchanged
edge lacked jointly canonical sites**, and **1 malformed/duplicate aligned
FASTA**. There were **521 invalid original edges** among mask-checkable
species; this does not include edges in the malformed alignment, and no
survivor edge was removed or rewired. The retained species names matched
the pre-existing survivor roles CSV **259 out of 259**.

The exact post-mask files held privately in Library:
`/TTF/historical_host_memory/confirmatory_mask_259_20261011/`
- `survivor_roles_v0.1.csv` SHA256 `fbeee8baa0529e162c0331a65717ec313e625e057e33ecafc48e99fa2e3c2535`
- `survivor_edges_v0.1.csv` SHA256 `91b59ca469bcf80bb2e5e63d609eabd6dfaacc1643bde82b1455c62cfeafeabf`
- `independent_mask_ledger_v0.1.json` SHA256 `b8b80e3edeae07d0a61606ec800d4a0c61296ebb06b379df3a1140ca3ecfa3d5`
- `independent_survivor_synthetic_v0.1.json` SHA256 `ba5bb37514feae9a2031fd2d3d5af00bbcc8f51308d97a8314722eb7b1dc4482`
- both private replay scripts are saved alongside these receipts. The original ZIP is retained separately.

**Byte-level identity:** the two retained CSVs are *exactly*, byte-for-byte,
the previously frozen 71,520-edge predictor CSV and 641-species assignment CSV
respectively, filtered by the surviving 259 species. There are
**20,362 surviving original edges**; the 62 failed species account for
21,230 other original confirmatory edges. The original confirmatory
panel had 41,592 edges.

## The two surviving gates

The **unchanged** predictor information criterion evaluates as follows:
259 species (required >=200), median independent host-history information
`0.240586632194` (required >=0.10), 244/259 species >=0.05
(unique fraction **0.9420849420849421**, required >=0.70), and 259/259
with condition number <=30 (required >=0.90). **PASS.**

The pre-frozen independent survivor seed namespace
`historical-host-memory-survivor-v0.1` was then used on the exact
259-species / 20,362-edge external predictor design, with the unchanged
1,999 null-reference / 500 evaluation null / 500 positive worlds and
endpoint-safe leave-two-locality-out IBD correction.

- Null false-positive rate: **30/500 = 6.0%**, Wilson 95% upper **0.084361**,
  required <=0.10 — PASS.
- Positive power: **485/500 = 97.0%**, Wilson 95% lower **0.951096**,
  required >=0.80 — PASS.

These are **simulation** power/error properties at the previously declared
positive world effect, **not** an observed insect genetic effect.

## Source provenance and response firewall

An aggregate public receipt is
`benchmarks/frozen/historical_host_memory_private_confirmatory_mask_259_survivor_replay_v0.1.json`.
Its public CI validates canonical hashes, arithmetic, subset counts, Wilson
limits and the response firewall. It **does not** reprocess the original
private ZIP on GitHub. An independent local Boolean-mask replay *did*
reprocess the 321 confirmatory originals, with exact 259-species agreement.

The original historical host-climate model is a **proxy for climate
opportunity relative to present recorded host taxa**, not verified plant
occupancy or host interaction at 21 ka. Cultivated-host flags remain
descriptive and cannot be post-hoc used to revise the frozen panel.

**Still permanently closed:** development-panel empirical responses and
confirmatory nucleotide identities, p-distances, empirical post-IBD turnover,
or the coefficient relating historical host opportunity to genetic structure.
The next authorized *review* is a separate exact-source, one-shot identity
opening authorization. Do not infer ecological support from the synthetic
qualification and do not unblind through the reproducibility receipt.
