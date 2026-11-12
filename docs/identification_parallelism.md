# QC identification scheduling

Both `qc-illumina` and `qc-nanopore` schedule SpeciesFinder/ANI after
read preparation and assembly. Up to two independent sample jobs run concurrently.
The actual count is capped by sample count, the run's CPU budget, and available
RAM. Each job gets an equal integer share of the available threads. Smaller
hosts and unrecognized database layouts retain serial processing.

The conservative memory scheduling estimate is at least 32 GiB per job, or
twice the SpeciesFinder component inventory size plus 4 GiB, whichever is larger.
This is an estimate, not a measured peak or hard memory limit. The existing
CPU/cgroup-aware memory helper reserves host headroom; explicit Illumina memory
budgets remain respected. Resource allocations are printed in the QC log and
saved in `qc_out/identification_resources.json`.

Jobs own separate sample workspaces and read the common reference indexes without
modifying them. A failed nomination affects only that sample. Report rows remain
in sample-sheet order even when jobs finish out of order. Batch reconciliation,
CheckM, publication and optional broader confirmation/trees retain their existing
ordering and scheduling. Identification thresholds and biological interpretation
are unchanged. Legacy KmerFinder and Lambda confirmation stay serial.

Optional progress events count Lambda as completed when its confirmation finishes
(whether QC passes or fails), or failed when the confirmation has a technical
error. Illumina samples not identified because read QC failed are explicitly
skipped. Processing completion is not a QC verdict.

Focused tests cover resource limits, overlapping jobs, failure isolation, and the
actual control/progress branches of both CLI scripts. No production speed-up
factor is claimed without a measured benchmark.
