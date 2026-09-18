# Frozen fresh genetic TTF Phase-4 result

Status: **EMPIRICAL_PHASE4_COMPLETE**

Source of truth: `benchmarks/frozen/genetic_phylogatr_phase4_empirical_handoff_v0.1.json`.

- exact source archive: 274,988,692 bytes; SHA-256 `5a0fd9ac25893c749d14186fbcce4a46b99163c9d810b36e40eebce7bece61a5`;
- Phase 1: 250 species, byte-identical manifest `daeedabbeb99fa304568f0f2b2e1c46e0dd6221fbac02b0491faf2435ef7b7b3`;
- Phase 2: 211 survivors, inherited 103 training / 108 evaluation split;
- Phase 3 cross-species Gate-D: PASS (max private-null Wilson95 upper 0.0680778; shared-A2 Wilson95 lower 0.971387);
- Phase 3 self-detectability: PASS (null Wilson95 upper 0.0889524; private-A2 Wilson95 lower 0.992376);
- primary `place_beyond_ibd`: `T = 0.0325655`, profiled-private `p = 0.739261`, not significant;
- empirical within-species self diagnostic: `S = -0.226098`, upper-tail `p = 0.00899101`, significant relative to the frozen structural-null reference;
- secondary total genetic transfer: `0.0155016`, descriptive only;
- frozen decision: `LINEAGE_CONDITIONED_SPATIAL_STRUCTURE_WITHIN_TESTED_DOMAIN`.

The self statistic is not tested against zero, so its negative raw sign is not a negative biological effect label. The supported interpretation is null-relative within-species detectability together with non-detectable cross-species transfer. This does **not** prove universal zero transfer or identify a particular historical mechanism.

No sequence identities or edge-level genetic-distance vectors are serialized here. Post-result retuning and result-selection reruns are forbidden by the frozen one-shot contract.
