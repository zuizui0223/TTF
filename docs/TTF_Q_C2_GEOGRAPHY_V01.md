# TTF-Q C2 — geography-specific information loss

C2 is frozen **before the bound C1 workflow result is opened**. It is nested inside C1 rather than chosen in response to C1.

C1 asks how much historical climate-displacement similarity remains after current climate, lineage, occurrence-count imbalance, and exact source/target identity. C2 adds exactly one control:

> the directed fraction of target occurrence coordinates lying within 500 km of any source occurrence coordinate.

The 500-km radius is inherited from the already frozen B/C spatial-opportunity scale. C2 does not tune a distance threshold, does not discard low-overlap dyads, and does not use genetic edge geometry.

Let `U1` be C1's total unique historical-relation variance fraction and `U2` the corresponding fraction after adding external geographic co-opportunity. Because the models are nested,

```
geography-retained fraction = U2 / U1
geography-attributable fraction = 1 - U2 / U1
```

The first quantity is the fraction of history-specific ecological information that remains after geography. The second is the portion of C1's apparent history-specific information that can be reproduced by present spatial co-opportunity.

Both development and confirmatory panels are reported. Neither panel may be selected by its result.

This is still entirely response blind: nucleotide identity, pairwise genetic distance, `T_st`, and `beta_hist` remain unopened, and relational v0.3 remains `CLOSED_NO_EVALUABLE_TEST`.

Implementation:

- rule: `docs/supporting/ttf_q_c2_geographic_opportunity_v0.1.json`
- geometry core: `src/ttf/external_geographic_opportunity.py`
- design builder: `scripts/attach_ttf_q_external_opportunity.py`

The resulting design is compatible with the common `run_ttf_q_characterization.py --kind historical` runner, so the same C2 geometry can later receive the predeclared detectability surface without redefining the ecological relation.
