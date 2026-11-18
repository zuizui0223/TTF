# Penstemon four-origin holdout: source-only eligibility rule (2026-11-15)

**Registered before opening any original floral phenotype values.** The 52 BioSample × taxon source receipt is already authenticated and is the sole data source for this gate. The 34 phenotyped accessions are **not mapped to the 52 genome accessions**. No trait, fitness, genomic variant, barrier, or ancestry response is fitted.

## Four biologically defined blocks (paper-defined, not selected from phenotype outcomes)

Stone et al. (2026), bioRxiv DOI 10.64898/2026.07.21.739878, inferred four independently derived hummingbird-syndrome origins within *Penstemon* sect. *Habroanthus*. Define four indivisible source units containing derived species and their reported sister comparator(s):

- B1: *P. barbatus* + *P. virgatus* (both species included in one block despite reported east/west genealogical intermingling).
- B2: *P. cardinalis* + *P. glaber*.
- B3: *P. eatonii* + *P. laevis*.
- B4: *P. labrosus* + both *P. deaveri* **and** *P. putus* (both bee comparators from the same sister grouping, not two independent shifts).

Use exactly four leave-one-origin-block-out folds. Exclude all source taxa outside these four blocks from **both** arms; do not silently treat unrelated species as bee comparators. No row-level or random individual splits. Unknown accession×phenotype membership means only **source-genome coverage**, not TTF-P eligibility, can be determined here.

## Full genome source authority

[52-authenticated NCBI BioSample taxon receipt](penstemon_habroanthus_52_ncbi_biosample_taxon_receipt_20261115.json), immutable Git blob SHA `acdf2997e4d338fcb879d0323dc4c6a54f4fe413`, verified original ENA PRJNA1479764 IDs via NCBI source accession lookup and successful run #43710257923.

The manuscript reports 34 phenotyped accessions from among 52 sequenced; the source provides **no authentic specimen-level phenotype identity table in this repository yet**. For any block with n genome specimens, an upper bound on phenotyped specimens is min(n,34), but its lower bound is max(0, n−18). These are **sharp combinatorial bounds**, not estimates. Sum across the four blocks lower bounded by max(0,34−source-samples-outside-blocks). Do not interpret lower bounds as measured flower observations.

## Admission

A source-only positive gate requires:
1. 52 original BioSample rows, no duplicate BioSample accessions and 21 literal NCBI taxon labels (20 ingroup + palmeri outgroup);
2. every reported hummingbird origin and every bee sister comparator in its *own* block represented by at least one genome accession;
3. groups disjoint at species string and accession level; sum of all four blocks does not exceed 52;
4. all four test folds retain both hummingbird and bee **genome source taxa** in the held-out block **and** training complement.

**Stop** before TTF-P, TTF-G or CREST predictions if exact 34 phenotype accessions, trait endpoint, and independently supported genealogy are not source-matched. Four independent origins means a maximum of **four evolutionary evaluation units**, irrespective of 52 genomes or multiple sampled flowers.

The first frozen TTF-Q main manuscript, existing P0/P1 Senecio results and PR #156 are untouched.
