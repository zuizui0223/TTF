# One-shot retry: official Penstemon Table S1 visual-only source diagnosis (2026-11-18)

The prior source-only render jobs are documented as failing for distinct non-biological causes: missing `pypdf` dependency on the first attempt, then timeout on the official Europe PMC supplemental ZIP on the second. The reader now retries that same original endpoint at most three times while validating PDF SHA256 `3682720dad4f0621bebef3be7014aa667067c5f0dce9160ca8f4c4a758e7f943`.

**No phenotype value or variant data are opened.** The only outputs are verified page-11 Table S1 and page-12 pre-S2 crop PNGs, to identify whether the original PDF actually contains a missed 52nd voucher or a layout-induced read error. This is not permission to change the strict 52/34 source gate based on assumptions.

If a third render still fails, stop retrying the same network path and log the remote-access issue; do not relax the source identity or count criterion.
