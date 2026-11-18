# Penstemon 2026: prespecified full 52-BioSample taxon provenance gate (2026-11-15)

**Status: SOURCE-ONLY, prospective before opening the remaining 47 taxon metadata summaries.** First-five fixed pilot has already succeeded at [Actions #43710076449](https://github.com/zuizui0223/TTF/actions/runs/43710076449), verifying source-accession/NCBI numeric UID conversion, and documenting 3 *P. virgatus*, 2 *P. barbatus* from the first five original ENA BioSample accessions. That pilot result is opened and is not a new independent test.

## Original source and fixed 52-sample protocol

- BioProject `PRJNA1479764`, original 2026 *Penstemon* *Habroanthus* genome sequencing (Stone et al. 2026 bioRxiv DOI `10.64898/2026.07.21.739878`). Independent ENA `read_run` metadata establishes exactly 52 unique run IDs and 52 unique `SAMN...` BioSample accessions, distinct from the 34 floral-phenotyped accessions.
- Freeze an exact source-only algorithm: use the full lexicographically sorted 52 valid ENA BioSample IDs; resolve them with one NCBI ESearch query of **`SAMN...[Accession] OR...`**; require 52 unique numeric Entrez UIDs; fetch ESummary by numeric UIDs in one request; require 52 unique returned BioSample accessions with exact identity equality to the original ENA set. No exclusions, taxonomic synonym merging, accession dropping, or source-specific threshold tuning.
- Read only NCBI `accession`, `organism`, `title`-presence and field names, **not** raw XML/sampledata attributes, genotype/sequence, flower chemistry, geography, or mating measurements. Strict report is frequency counts of **source taxon strings** with original NCBI accession and source mapping reproducible in audit code, plus whether the following source species pairs have at least one original BioSample on each side:
  - `barbatus` × `virgatus`
  - `cardinalis` × `glaber`
  - `eatonii` × `laevis`
  - `labrosus` × `putus` (alternative `deaveri` separately)
- These contrasts are the published 2026 **four candidate independent evolutionary origins**, not four new independent trait responses. Source taxon coverage alone does not demonstrate four measured flowering traits or independent inferred causes, and `34` floral measures must not be equated to the `52` sequenced samples without a sample × phenotype manifest.
- Report all 52 source accession mappings as a metadata-only receipt, or stop `NOT_EVALUABLE` on any absent, conflicting or invalid ID; no residual partial subset is silently accepted.

## Confirmatory claim boundary

Even if all four source groups are present, TTF-P/TTF-G/TTF-RI remains untested without a strictly joined public original phenotype / neutral phylogeny / other outcome resource. The older 2023 45-tip CDS/window trees and the 2026 52-genome BioProject may refer to related but nonidentical specimens; do not map `CAW...` to `SAMN...` by species alone. The wide 239-species Dryad nuclear tree remains publicly metadata-visible but original tree bytes give HTTP 401.

This study does **not** alter the frozen TTF main science. The only authorized operation under this gate is source taxonomy and original genome-sample provenance; preserve all original source labels.
