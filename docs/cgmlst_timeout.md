# cgMLST runtime safeguard

`bactipipe relate` gives each genome's cgMLST invocation 900 seconds by default.
Use `--cgmlst-timeout SECONDS` to select a different positive limit. This bounds
the entire invocation, including conda, Prodigal and KMA descendants; cancellation
of the parent cohort also stops these descendants.

A timeout leaves that sample's cgST unavailable (`NA`) and does not block other
samples, ANI, SKA or report creation. The time limit is an operational safeguard,
not a diagnostic threshold. The original KMA parameters and successful typing
results are unchanged.

Evidence is retained in `cgMLST/<sample>/cgmlst.stdout.txt` and
`cgMLST/<sample>/cgmlst_status.json`; the latter records execution status, elapsed
seconds, timeout and any failure reason. Cohort `analysis_result.json` carries
these failures as sample-scoped warnings and records the configured limit.
`--resume` reuses successful cgST values, but retries unavailable values.
