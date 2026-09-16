# SpeciesFinder assembly-fragment limit review

Date: 2026-09-12. Status: source implementation; deployment and full QC reruns pending.

## Decision

Use a maximum **20,000 bp** per non-overlapping assembly fragment for
SpeciesFinder nomination on both Illumina and Nanopore. This is a locally
benchmarked engineering setting, not a diagnostic cutoff or universal optimum.
No reference database, ANI threshold, aligned-fraction gate, assembly or
published run was changed by this comparison.

The upstream KMA whole-contig alignment was sensitive to circular chromosome
origin/strand representation. Direct skani comparisons confirmed all six
P2533317 assemblies against CP040464.1 at 99.18% ANI, including those with no
intact-query SpeciesFinder nominations. Lossless fragmentation resolves the
observed nomination failure without changing the sequence used for ANI.
The exact upstream seed/alignment heuristic was not instrumented; this is not
evidence for a documented hard maximum FASTA record length.

## Comparison

Twenty cases: seven non-control samples from each of `X_250821_8_SEQ10_LMO` and
`X_260910_8_SEQ10_LMO`, four Illumina controls from `X_260903_32_SEQ04_LMO`, plus
two half-rotated Campylobacter variants. Thus there are 18 distinct samples, not
20 independent biological specimens. Eighteen cases have an expected target;
the two P2631472 samples are unresolved controls. Expectations are review context,
not independently validated diagnostic truth.

All 100 jobs completed using SpeciesFinder 4.0.4/KMA 1.6.13 and the immutable
`BactiPipe-SpeciesFinder-v4.0.4-20260112-3a72a1e6a9c7` database. Each comparison
retained the existing cap of 20 taxa and three exact templates per taxon; skani
always received the original intact assembly.

| Maximum query length | Target nominated | Confident ANI calls / 20 cases | Additional nominated taxa | Total nomination time |
| --- | ---: | ---: | --- | ---: |
| Intact contigs | 9/18 | 6 | None | 445.56 s |
| 5,000 bp | 18/18 | 17 | Mannheimia in unresolved Illumina P2631472 | 438.00 s |
| 10,000 bp | 18/18 | 17 | None | 410.84 s |
| **20,000 bp** | **18/18** | **17** | **None** | **413.11 s** |
| 50,000 bp | 18/18 | 17 | Acinetobacter in Salmonella P2633261-5 | 387.14 s |

Times sum individual nomination durations from a two-worker comparison, not wall
time; minor differences are not a replicated performance experiment. Additional
nominations did not become confident identifications and are not proven
contamination. The 20-kb setting is the largest tested size with no additional
taxon in this panel, recovers all target nominations and preserves every previous
confident species call. Both circular-origin variants are consistent. Neither
species nor reference caps are exhausted across taxa at that setting; within-taxon
template ranking can change, so a changed best ANI need not imply changed sequence.

Key findings at 20 kb:

- All six original Campylobacter samples confirm C. hyointestinalis at 99.18% ANI.
- All six Salmonella samples confirm S. enterica at 99.88–99.93% ANI. Intact
  nominations had missed three completely and selected insufficient reference
  evidence for the other three in this assembly-only comparison.
- Illumina E. cecorum P2627591 and P2630704 remain confirmed at 95.89% and 95.30%.
  The latter changes from 95.54% with intact nomination because the nominated
  reference set changes, not because its assembly was fragmented for ANI.
- Illumina PA16802CY260022 retains M. haemolytica at 99.65% ANI.
- P2631472 remains unresolved on both platforms; no unsupported ID is created.
- Nanopore P2529124 now nominates E. cecorum but remains unresolved: the highest
  tested ANI is 93.83%, with 53.79% reference and 51.84% query aligned fraction.
  None of its three references passes ANI >=95% and both aligned fractions >=65%.

These are nomination/ANI tests, not complete QC or contamination benchmarks:
no FASTQs, assemblies, CheckM, competitive read mapping, virulence or typing were
regenerated. Existing published run results can therefore differ in other ways.

## Implementation and audit

`bactipipe/speciesfinder_query.py` generates numeric fragment IDs and retains
source headers, contig numbers, lengths and zero-based half-open coordinates.
Short contigs retain their sequence; headers are standardized, not treated as
biological evidence. Ambiguous nucleotide bases and final tails are retained.
Malformed or empty FASTA fails explicitly. A new output directory is required;
the helper cannot overwrite the original assembly or an earlier audit.

`species_workflow.py` applies this only to the SpeciesFinder database path, never
to the legacy KmerFinder or reviewed-reference workflow. Nomination gets the
derived query; skani, placement and read confirmation get the original assembly
and existing inputs. The derived query is deleted after successful nomination.
Analysis provenance schema 4 binds the method, limit and hashes; the database
identity and BactiPipe's pre-stable application version remain unchanged.

There is no identification-result cache in this execution path: each invocation
stages a new result and replaces the old directory only after success. Selective
sample reanalysis likewise stages a replacement; it cannot reuse an old intact
nomination merely because an assembly is unchanged. No automatic reanalysis of
historical results is introduced.

## Evidence and next step

Private comparison root:
`/home/ngslab/temp_data/speciesfinder-fragment-limits-20260912/`.
Retain `comparison.tsv`, `summary.json`, `manifest.json`, source receipts,
per-case `result.json`, compact tool outputs, logs and status. Large derived
query FASTAs and extracted candidate references were removed after each case's
successful checkpoint. No FASTQs were downloaded for this work.

Focused regression tests cover boundaries/tails, ambiguity, duplicate headers,
losslessness, determinism, malformed inputs, both-platform routing, original
assembly ANI, cleanup and failed-attempt preservation. Full-workflow source-staged
Docker smoke checks use one previously missed Nanopore Campylobacter and one
confirmed Illumina E. cecorum, without changing production results.

After deployment is authorized, rerun `X_250821_8_SEQ10_LMO` and
`X_260910_8_SEQ10_LMO` and review their complete QC/secondary-evidence outputs.
Do not activate the candidate database or weaken identification gates on the
strength of this engineering comparison alone.
