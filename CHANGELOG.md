# Changelog

All notable BactiPipe application changes will be documented in this file.
Releases follow [Semantic Versioning](VERSIONING.md), and Git tags use the
`vMAJOR.MINOR.PATCH` format.

## Unreleased

- Correct detect organism-list parsing and show a plain, copyable list of readable
  organism names instead of tool/database details; accept legacy group aliases.
  Unsupported AMRFinder organism names now select ABRicate/NCBI without `-O`;
  Enterococcus uses the appropriate VirulenceFinder subsets, with VFDB fallback
  for unsupported virulence organisms. Preserve actual database provenance and
  propagate virulence-tool failures instead of silently reporting no hits.
- Simplify ordinary CLI help to five analysis commands; expose placement/rendering
  through identify options and SpeciesFinder preparation/validation through
  update-db speciesfinder. Retain historical command aliases for compatibility.

- Convert optional primary-QC workbook R/D requests to downstream cohorts with
  organism-aware defaults and reference selectors. Multiple references share one
  typing pass, a single SKA cohort calculation, and consolidated SNP/ANI column pairs
  in TSV/PDF reports; accept repeatable CLI reference flags and retain per-reference
  evidence in the relate-v2 contract.
  Sort downstream reports by supplied display names with worksheet-order fallback;
  accept a reference Specimen description and the QC NB Index header alias.

- Emit opt-in timestamped QC/assembly milestones and bounded resource/workload
  metadata for execution-duration auditing; forward parallel Illumina worker
  events live without changing ordinary CLI output or analytical results.

- Allow named GTDB MAG/environmental representatives to establish species
  confirmation; retain resolved taxonomy, ANI/AF and competitor safeguards
  and preserve source category in the reference audit.

- Remove redundant nearest-reference text from tree leaf annotations; retain
  the blue highlight and legend, accession, GTDB label, and outgroup annotation.

- Adopt ANI >=95% and max(query AF, reference AF) >=50% for new species decisions,
  matching GTDB-Tk's AF convention on both platforms. Record the aggregation mode,
  retain both measured AFs, preserve historical both/65% evidence, and update ANI
  plots/CLI controls. Reference eligibility and competing-species safeguards remain.

- Label Illumina quality/depth log messages by sample and report measured yield
  depth after preparation instead of incorrectly claiming assembly is pending.

- Fix SPAdes inheriting the Illumina QC worker's one-thread OpenMP cap despite
  its allocated --threads value. Override OpenMP only for assembler children,
  retain Python/BLAS limits, and verify effective threads in real-tool smoke tests.

- Link marker-tree and ANI figures through a bounded, exact-accession reference
  panel with shared IDs, explicit ANI-only references, and missing-value labels.
  New compact trees retain available top ANI nominees and a defensible outgroup.
  Retain low-alignment evidence for illustration without relaxing species gates;
  support paired SVG, PDF and PNG presentation.

- Emit bounded, structured execution milestones and sample completion events for
  Nanopore/Illumina QC, assembly, identification and final QC. The web integration
  keeps per-attempt progress distinct from QC Pass/Fail and from final species
  investigations. These optional events do not change analytical decisions.
- Use at most six concurrent Flye assemblies with ten threads by default, subject
  to CPU/memory bounds and explicit overrides; report completion-order progress.
- Add broader GTDB whole-assembly ANI confirmation before optional bac120 placement;
  confident ANI bypasses tree construction. Add compact retained-evidence tree and
  ANI visualizations with `render-placement`, without inventing outgroup/bootstrap
  support. Preserve original assemblies, reference identities and both ANI fractions.

- Bound concurrent Illumina read QC/assembly for both SKESA and SPAdes by sample,
  CPU and estimated RAM budgets (up to three jobs, twenty threads each by default).
  Keep identification serial and CheckM/reporting batch-wide. Add per-sample worker
  audits, retain failed assembler logs, and include resource settings in assembly
  checkpoint identity. New CLI controls: --sample-jobs, --sample-threads, --memory-gb.

- Use lossless, non-overlapping 20-kb assembly fragments for SpeciesFinder
  nomination on both platforms to avoid origin-sensitive long-contig misses.
  Keep original assemblies for ANI, QC and downstream analysis; retain fragment
  checksums/coordinates, remove the temporary query after successful nomination,
  and advance analysis provenance to schema 4. Identification thresholds and
  candidate database activation are unchanged.

- List Minimum depth immediately before mean quality in PDF QC thresholds:
  one Nanopore value or the Illumina 30/40/50X policy plus E. coli adjustment.
  Verify shared Lambda control reporting for future Illumina runs; do not add
  synthetic control rows to historical runs without Lambda inputs.

- Conditionally filter Illumina reads with fastp when combined raw base-weighted
  mean Phred quality is below minimum Q + 1. Use selected reads consistently for
  QC, assembly and read evidence; retain full-precision raw/selected metrics.
  Assess yield depth against actual assembly length, with 30/40/50X quality tiers
  and +10X only for confidently identified E. coli. Remove the catalog-size
  pre-assembly depth gate and bind assembly/CheckM reuse to their inputs.

- Compact QC PDF headers into aligned Run Metadata, Tool Versions and QC Thresholds
  columns with matching embedded heading fonts. Display retained Lambda identity
  and genome coverage together, with independent failure coloring and unavailable
  historical coverage explicitly marked N/A; provisional thresholds remain 99%/90%.

- Use assembly FASTA exclusively for SpeciesFinder candidate nomination on both
  Illumina and Nanopore; retain ANI validation and read confirmation of nominated
  taxa without searching reads for additional unassembled organisms.

- Add `qc-nanopore --input-run-name` so controlled benchmark reanalysis can
  reuse an existing physical read directory while publishing under a new run
  name.
- Integrate SpeciesFinder 4.0.4 and KMA 1.6.13 as the broad candidate-nomination
  layer for the inactive species-identification workflow. Exact nominated KMA
  templates are extracted for skani ANI/aligned-fraction confirmation; no
  SpeciesFinder score or coverage value is treated as a diagnostic call.
- Add resumable official database retrieval, immutable SpeciesFinder database
  preparation/validation, bounded candidate limits, fail-closed upstream-output
  checks, a row-level source/taxonomy audit, exclusion of MAG/SAG/single-cell
  templates from diagnostic ANI, and explicit preservation of clinically
  important *Streptococcus equi* subspecies.
- Keep the legacy KmerFinder production path and previous custom Kraken
  candidate path available until both-platform scientific validation approves
  a replacement.
- Add an inactive evidence-based species-identification workflow with skani
  ANI/aligned fraction, conservative platform-specific mixed-signal assessment,
  competitive mapped-read/breadth confirmation, and
  optional GTDB bac120 placement for unresolved non-mixed samples.
- Add taxid-controlled reporting groups for ambiguous complexes while keeping
  familiar clinical names such as *Escherichia coli* in user-facing results.
- Compare historical Salmonella serovar labels at species resolution without
  representing ANI as serovar confirmation, while preserving clinically
  important distinctions such as *Streptococcus equi* subsp. *equi* versus
  subsp. *zooepidemicus*.
- Prevent conserved reads shared by clinically retained subspecies under one
  NCBI species parent from creating a false mixed-culture result; the
  subspecies remain distinct and are resolved by whole-assembly ANI and
  bidirectional aligned fraction.
- Reconcile weak secondary evidence across the complete run, retaining
  donor-linked possible carryover without hiding substantial or donor-free
  secondary signals. Because this is an isolate workflow, only strong residual
  evidence after cohort review is labeled suspected culture contamination;
  weaker evidence remains review-only and is never called a mixed infection.
- Require independent CheckM corroboration before competitive short-read
  mapping can label an Illumina isolate as suspected culture contamination.
  Preserve a confident primary ANI identification separately when isolate QC
  supports contamination, and keep low-ANI taxa unresolved rather than
  weakening the global species threshold.
- Defer QC-run bac120 placement until cohort reconciliation and repeat that
  decision after selective sample updates, preventing carryover from
  incorrectly suppressing an eligible unresolved-sample tree.
- Add an explicit `--gtdbtk-in-memory` mode for controlled placement runs on
  high-memory hosts while retaining disk-backed scratch as the default.
- Keep KmerFinder out of the new workflow while retaining its parser and the
  existing production path for legacy comparison and rollback.
- Add a single reviewed-manifest database contract and immutable,
  non-activating `plan-identify-refs` planner and dual skani/Kraken2
  `build-identify-db` candidate builder.
- Add `plan-identify-inventory` to derive an auditable, bounded proposal of
  controlled veterinary species and same-genus GTDB representative neighbors.
- Add checkpointed `fetch-identify-reports` retrieval with per-batch checksums
  and explicit missing-accession audits.
- Add current-accession reconciliation with canonical fallback promotion and
  fail-closed exclusion of suppressed, unnamed, or taxonomically conflicting
  assemblies, using pinned NCBI merged-taxid resolution and deliberate direct-
  taxid support for controlled subspecies.
- Add exact-taxid controlled-gap recovery and checkpointed per-assembly NCBI
  FASTA retrieval with multi-contig preservation, length verification,
  checksums, resumable batches, and reviewed-manifest generation.
- Add explicit full-manifest and immutable candidate-database validation commands.
- Bind the reviewed ANI and candidate-only Kraken2 indexes into one versioned
  database directory; omitting it preserves the production KmerFinder path.
- Isolate per-sample QC, assembly, and organism-identification failures so other batch samples can finish.
- Add repeatable `--only-sample` selectors, transactional `--update-existing`
  publication, and `sample_failures.tsv` for targeted QC reanalysis.
- Keep Nanopore sample-folder/barcode resolution aligned with the web launcher.
- Record deterministic application source and formal database identities in downstream result provenance.
- Include the formal BactiPipe database identity in QC and downstream reports.
- Resolve controlled organism names once per sample and remove the silent 5 Mb
  fallback used by Nanopore QC and Flye.
- Run unknown-organism Nanopore samples without size-dependent Filtlong or Flye
  arguments; keep Illumina assembly size-independent.
- Recalculate coverage and its verdict after a previously unknown organism is
  identified, and record exact per-sample size decisions in
  `sample_execution_plan.tsv`.
- Bound Nanopore QC and assembly workers by sample count, CPU allocation, and
  available-memory limits.

## 0.1.0 - 2026-08-31

Initial formally versioned development release.

- Provides Illumina and Nanopore QC and assembly workflows.
- Provides cohort-level `relate` strain-relatedness analysis.
- Provides `detect` virulence and antimicrobial-resistance analysis.
- Records the authoritative BactiPipe application version in generated reports
  and downstream provenance.
- Exposes the application version through `bactipipe -v` and
  `bactipipe --version`.
