# External reference acquisition: interrupted-transport recovery

This is a transport-only correction to PR #145, not a revised ecological design or a new confirmatory analysis.

The initial workflow limited the whole job to 25 minutes while keeping every completed species in memory until all 641 futures returned. Thus a job-level timeout could discard previously retrieved observation references. Runtime progress is not a source-backed saved result.

The correction persists each species atomically before processing the next species. Checkpoints bind the exact original geometry, the unchanged acquisition contract, the 641-species identity, the original internal quadrature implementation, and the exact query/filter/thinning/retry function source. Both successful and failed species attempts are retained; normal resume requests only previously unattempted species. Unexpected exceptions cannot invalidate already saved species. Corrupt or differently bound checkpoints fail closed.

A soft 180-second submission budget pauses new requests and lets at most four current species finish their existing bounded retries. The workflow retains partial snapshots, with all 641 species represented and unattempted species explicitly labelled NOT_ATTEMPTED_IN_CURRENT_SNAPSHOT (not zero observations). The 45-minute job ceiling leaves room for in-flight bounded requests and artifact delivery. Further runs resume, not replace, these checkpoints. There is no automatic retry of a completed REQUEST_ERROR attempt; any such retry needs a separately recorded technical decision.

No species membership, host association, taxonomy condition, dataset, date range, uncertainty cutoff, three-window rule, 10-km thinning, 50-km primary separation, support threshold, or maximum of four workers changes. A local source comparison verified byte-identical function bodies for window_offsets, unit_vectors, nearest_km, accepted_key, clean_record, thin_records, get_json and acquire_one. The original 26 tests and 10 new interruption/resume tests pass (36 total). No claim is made here about new empirical counts or CI completion before an actual saved snapshot is read.

Results remain exploratory ecological data support only. Host observation matching, real climate recalculation, and genetic validation are not performed by this workflow. No old NOT_EVALUABLE_EMPIRICAL_RESPONSE decision is changed.
