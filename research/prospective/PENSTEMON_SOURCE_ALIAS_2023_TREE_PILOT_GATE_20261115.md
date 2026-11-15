# Penstemon original voucher identity gate: 2026 NCBI versus 2023 author tree (2026-11-15)

**Predefined source-only diagnostic, before opening optional NCBI submitter-original sample identifiers.** The full 52 original SAMN accessions and their 21 source taxon labels have been authenticated via NCBI official ESearch→ESummary, Actions #43710257923. That is genome-source support only, not 34 phenotyped accessions.

**First-step fixed pilot:** the **lexicographically first five** of the 52 ENA-source SAMN IDs (SAMN60992931..SAMN60992935), without selecting samples by plant species, trait, or evolutionary origin. Fetch ESummary after ESearch as before, and inspect only the public `sourcesample` field (if a simple short source-key string), as well as the *structural type* of the `identifiers` field. Do not decode or export `sampledata` or other unrequested BioSample properties. Match a `sourcesample` to the author 2023 public nuclear CDS tip names **by exact string identity only**. If matching requires guessed normalization/alias expansion or if `sourcesample` is not the original voucher ID, record NOT_EVALUABLE; do not claim same individuals.

**Independent 2023 source:** `lcwheeler/penstemon_novaseq`, last modified 2023-06-02, CDS tree 45 original tip names, blob `ca1d2072e7df72649a8ec8699e31d8397fc46cae`. The historical genome source is not automatically the final 2026 author dataset. A species label match alone is never an individual/voucher match.

**Long-term biological target:** 52 genome vouchers → 34 directly measured floral vouchers → independent ancestry. Even successful genotype-voucher keys would NOT establish measured floral phenotypes or TTF-G causality. The 2026 study code/phenotype source is described as forthcoming after journal publication, not confirmed available.

No biological sequence, flower trait, fitness, reproductive-isolation outcome or source phylogenetic topology is opened under this metadata-only pilot. Keep all findings in Draft PR #159; frozen TTF science stays intact.
