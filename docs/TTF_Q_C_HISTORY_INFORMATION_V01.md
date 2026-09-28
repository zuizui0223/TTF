# TTF-Q C1 — history-specific ecological information

This is the first empirical target inside the fresh TTF-Q methods program.

It deliberately stops **before** geographic-opportunity reconstruction and
before any genetic response. It can therefore be evaluated from the
response-blind historical-relation design alone.

For each source-target dyad, C1 asks how much historical climate-displacement
similarity remains after conditioning on:

- current-climate similarity;
- class, order and family identity;
- occurrence-count imbalance;
- exact source identity;
- exact target identity.

The principal continuous quantity is
`total_unique_variance_fraction`: the fraction of standardized
historical-relation variance remaining after those adjustments.

C1 also reports a directly interpretable ecological contrast among dyads in the
upper quartile of present-day niche similarity:

- **high-current / low-history**: present similarity reached through dissimilar
  late-Quaternary displacement histories;
- **high-current / high-history**: present similarity accompanied by similar
  late-Quaternary displacement histories.

This is not a genetic-transfer test and does not alter
`CLOSED_NO_EVALUABLE_TEST`. Geographic opportunity is added only in the later
full TTF-Q qualification stage.

Runner:

```bash
python scripts/run_ttf_q_historical_information.py \
  --historical-design relational_historical_design_v0.1.npz \
  --panel development \
  --output ttf_q_c1_history_information.json
```
