# BactiPipe species identification

September 15 GTDB confirmation revision: MAG/environmental origin no longer
disqualifies an installed GTDB representative. Resolved species names/taxids,
clinically meaningful subspecies, ANI >=95%, either aligned fraction >=50%, and
competing-species checks remain required. Preserve source categories as evidence.
This supersedes older GTDB MAG-exclusion notes, not viral database policies or
the historical custom database builder. Changes await the next deployment;
published reports are not automatically rewritten.

## Status

The broad SpeciesFinder/KMA candidate-nomination layer, exact-template
extraction, ANI/aligned-fraction decision model, mixed-signal confirmation,
bac120 placement wrapper, and opt-in QC integration are implemented but are
**not active in production**. The official bacterial SpeciesFinder database is
retrieved separately, checksum verified, and prepared as an immutable inactive
version. The initial ANI/AF and mixed-signal thresholds remain validation
candidates. Existing KmerFinder behavior stays in production until a prepared
SpeciesFinder version passes the focused cases, complete regression cohort, and
independent Illumina/Nanopore diagnostic review.

The 2026-09-07 GTDB R232 inventory/current-NCBI reconciliation and controlled-
gap recovery and FASTA retrieval produced 4,218 validated references across
3,453 species. It represents 1,578 of 1,582 controlled catalog entries. Four
species remain explicitly unresolved because NCBI has no assembly passing the
quality and taxonomy gates. Eight conflicting assemblies were excluded rather
than silently relabeled. Candidate1 built and validated successfully but exposed
a systematic false-mixture call for pure *Escherichia coli* controls. Candidate2
corrected that decision defect and improved the repeated 11-case benchmark from
six expected, three manual-review, and two unexpected outcomes to eight expected
and three manual-review outcomes, with no unexpected known-positive calls. It
remains inactive pending broader diagnostic validation. Cohort-aware
donor/recipient reconciliation is now implemented and is being benchmarked
against the weak *Salmonella enterica* signal in C04 PrimeLyse and a strong
Salmonella sample from the same multiplexed run.

## Diagnostic result model

BactiPipe separates three questions that must not be collapsed into one score:

1. Does the assembly have a confident whole-genome species match?
2. Do the reads support a single isolate or more than one bacterial species?
3. If no species call is defensible, where does the sample fall among known
   bacterial species?

The primary species evidence is skani ANI plus query/reference aligned fractions.
On 2026-09-14 the user approved ANI >=95% and **max(query AF, reference AF) >=50%**
for new decisions on both Illumina and Nanopore. Either fraction can meet 50%; they
are not averaged and both need not pass. No rounding is applied before comparison.
This replaces the provisional requirement that both fractions reach 65%.

The AF rule matches [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/faq.html) and its
installed `external/skani.py` implementation (`max(af_r, af_q)`). [GTDB documents](https://gtdb.ecogenomic.org/faq)
the change from 65% to 50% at RS207. GTDB uses species-specific ANI radii (usually
95%, sometimes higher); BactiPipe retains the user-selected uniform 95% floor.
These are tool-aligned operating defaults, not a claim of completed diagnostic
validation or exact equivalence to the complete GTDB classifier.

Existing competitor rules (0.5 ANI points and 2 AF points), reference eligibility,
MAG restrictions, isolate-QC rules and clinical subspecies handling remain in place.
AF-based ranking/margins now consistently use the configured aggregation. New
outputs store `aligned_fraction_mode=either` with both measured AF values. Historical
thresholds lacking this field retain `both` and their saved cutoff (65% fallback)
when replayed; historical ANI plots retain the smaller-AF axis. New plots use the
larger-AF axis with a 50% line. The standalone identification command also supports
`--aligned-fraction-mode both` for deliberate legacy comparisons.

No existing reports were rewritten or runs requeued for this change. The planned
post-deployment comparison will prioritize all available non-Salmonella/non-E. coli
runs, with independent laboratory identification as review evidence. Thresholds may
be revisited after that review.

Source verification on September 14: 81 analytical tests and 50 focused web tests
passed, covering inclusive 95/50 boundaries, either AF direction on both platforms,
legacy replay, retained MAG/competitor safeguards, plot aggregation and shared-panel
consistency. Logs are `/home/ngslab/status/species-af50-20260914-tests.log` and
`/home/ngslab/status/species-af50-20260914-web-tests.log`. No build/deployment or
scientific reruns were started by this change; source parity and compilation passed.

As of 2026-09-12, SpeciesFinder receives only the assembly FASTA for both
Illumina and Nanopore isolate analysis, even when filtered FASTQs are available.
There is no FASTQ fallback or extra read-based candidate-discovery screen.
ANI validation is unchanged. Competitive mapping still uses filtered reads to
assess support for assembly-nominated taxa, not to discover additional unassembled
organisms. The source-staged fragmentation refinement uses deterministic,
non-overlapping fragments of at most **20,000 bp** for nomination only, preserving
all bases, short contigs and final tails. Original assemblies remain authoritative
for ANI, placement, QC metrics and downstream typing. Provenance schema 4 records
`speciesfinder_input_kind=fragmented_assembly_fasta` with the method, limit,
source/query/sequence checksums and counts. `speciesfinder_query/query.json` and
`fragments.tsv` retain zero-based half-open source coordinates, while the derived
`query.fasta` is removed after successful nomination. Each identification executes
in a fresh staging directory; prior intact-nomination results cannot satisfy a
new invocation. Failed execution preserves the previous published result.
Existing results are not rewritten automatically. See
[fragment-limit comparison](speciesfinder_fragment_limit_review.md).

SpeciesFinder/KMA proposes up to 20 taxa and up to three exact database
templates per taxon by default. BactiPipe extracts only those templates from
the same KMA index and runs skani against them. A SpeciesFinder score, query
identity, query coverage, or template coverage never establishes the final
species: only BactiPipe's ANI and bidirectional aligned-fraction gates can do
that. In KMA's alignment output, `Template_Coverage` is the aligned-path length
divided by template length, while `Query_Coverage` is the reciprocal ratio.
Insertions or length differences can therefore make either value exceed 100.
Both are retained under their upstream names as nomination evidence and must
not be relabeled as conventional genome breadth, assembly coverage, or ANI.

KmerFinder is not part of this new candidate path. Its parser remains available
only to compare already-produced legacy results during validation. Its query
and template coverage are preserved under those names and must never be
labeled percent identity.

Possible secondary organisms in an isolate are evaluated from filtered reads. SpeciesFinder/KMA
only nominates taxa; its score is never a positive identification or proof of a
contaminated culture. KMA `Score` is retained explicitly as a nomination score and is never
stored or displayed as a read count. Reads are then mapped competitively to one
extracted reference for every bounded candidate taxon. A secondary signal must
pass competitively mapped read-count,
mapped abundance, and distributed reference-breadth gates. The independently
assigned competitive read count supplies the support gate that the prior
Kraken path obtained from distinct minimizers. When the previous custom Kraken
candidate is selected, its stricter distinct-minimizer evidence is still
required. Missing required evidence produces an indeterminate review state; it
never proves or dismisses contamination.

Some named species cannot be resolved reliably by routine whole-genome ANI or
competitive read mapping. These are handled only through a controlled,
taxid-based reporting-group table. The first internal group contains
*Escherichia coli* and the four traditional *Shigella* species. Member
references, taxids, and evidence remain distinct in the audit; only ANI
competition and mixed-signal decisions are collapsed within the group. The
group label is supporting provenance, not the user-facing identification. A
known expected member retains its clinical name and uses its best passing
same-species reference for the reported ANI. Evidence from an unrelated
organism, such as *Salmonella enterica*, is still evaluated separately.

Illumina and Nanopore use separate provisional secondary-signal profiles. The
Illumina profile starts with stricter read, distinct-minimizer, and breadth
requirements because short exact matches are particularly vulnerable to
non-specific classification among closely related bacteria. These settings
must be estimated from pure isolates, defined mixtures, and background data;
they are not controlled diagnostic defaults yet.

After cohort cross-talk review, competitive mapping alone is still not enough
to label an Illumina isolate as suspected culture contamination. A strong
residual Illumina signal must also be corroborated by independent assembly QC;
the current implementation uses CheckM contamination above the configured run
limit. If CheckM is unavailable or remains within its limit, the secondary
signal stays `mixed_signal_review`. Nanopore retains its platform-specific
long-read evidence rule while both platforms preserve the same cohort review.

Within-sample mapping can establish that a secondary signal has distributed
genome support, but it cannot determine whether that signal came from the
sample. BactiPipe Isolate therefore never interprets a secondary as a mixed
infection. Final culture-contamination interpretation requires a run-wide
comparison first.
A weak secondary is eligible for `possible_cross_talk` only when the same taxon
is a much stronger primary signal in another multiplexed sample and the weak
signal is small relative to both that donor and its recipient's primary signal.
The evidence remains visible. Signals without a credible cohort donor remain
under secondary-signal review unless they also meet the stricter culture-
contamination evidence tier. That tier currently requires the secondary to be
at least 10% of competitively assigned reads, at least 25% as abundant as the
primary, and distributed across at least 50% of its reference, while the primary
must cover at least 70% of its reference. A strong signal is reported only as
`suspected_culture_contamination`, never as a definitive mixed infection. No negative control is required
for this donor/recipient comparison. The initial validation policy requires at
least 100 competitively mapped reads in the primary donor, a recipient-to-donor
ratio <=10%, and a recipient-secondary-to-recipient-primary ratio <=20%. These
are provisional software-validation gates, not approved diagnostic cutoffs.

The result states are:

- `confirmed_species`
- `no_confident_identification`
- `ambiguous_species`
- `mixed_signal_review`
- `suspected_culture_contamination`

`confirmed_species` supplies ANI to the legacy `% Match` field. A sample with
independently corroborated `suspected_culture_contamination` also retains a
confident primary ANI identification when one exists; contamination remains a
separate isolate-QC failure rather than erasing the primary organism. Other
states retain their closest-reference measurements in `closest_species.tsv`
while the main result reports `N/A` match. A known expected organism that
cannot be confirmed fails taxonomic confirmation; an Unknown expectation
remains `N/A` unless mixed-isolate evidence makes the isolate QC invalid.

A species-specific low ANI is not repaired by lowering the global 95% ANI
gate. For diverse taxa such as *Enterococcus cecorum*, the candidate database
must instead be expanded with quality-controlled representatives across known
clades. Until a reference passes the global ANI and bidirectional aligned-
fraction gates, the report remains unresolved and shows the closest supported
species for review.

The QC commands enable this path only when `--species-database-root` is
supplied. A SpeciesFinder database root binds the official bacterial KMA index,
its taxonomy TSV, exact template-name lookup, archive checksum, tool versions,
and component inventory into one immutable identity. An inactive database
additionally requires `--allow-candidate-species-database`, which is intended
only for controlled benchmarking. Omitting these options preserves current
KmerFinder production behavior until activation is approved. The older custom
reviewed-manifest/skani/Kraken2 candidate contract remains accepted for
comparison and rollback, but it is no longer the planned broad nomination
strategy.

The first immutable Candidate 4 database was prepared on 2026-09-09 as
`BactiPipe-SpeciesFinder-v4.0.4-20260112-3a72a1e6a9c7`. Its verified source
archive SHA-256 is
`3a72a1e6a9c77de3417b0d96f03cec3c91db2d2cdac04d2aaef6b51c0b3023a3`.
The database contains 36,230 KMA templates. The prepared source audit records
27,108 with a name and numeric species taxid, 553 excluded as explicit
MAG/SAG/metagenome or single-cell sources, 8,123 lacking a usable species name,
and 446 having a name but no numeric species taxid. Runtime validation further
keeps 6,335 placeholder labels such as `Genus sp.`, family-level `bacterium` or
`archaeon`, uncultured taxa, and unnamed symbionts out of the callable candidate
set, leaving 20,773 currently reportable templates. Every excluded row remains
visible in the taxonomy audit or raw SpeciesFinder results. Missing taxonomy is
failed closed until reviewed accession-level NCBI enrichment is available; it
is never replaced with a synthetic taxid or inferred from filename patterns.

## SpeciesFinder database contract

The source database is the official SpeciesFinder database archive. Retrieval
uses `scripts/fetch_speciesfinder_database.sh`, which resumes a partial archive,
verifies the expected byte count, writes a SHA-256 checksum and archive
inventory, and extracts only the bacterial KMA database and taxonomy file. It
does not activate anything.

Run the long retrieval in a detached session and inspect its durable status:

```bash
tmux new-session -d -s speciesfinder_db_fetch \
  'scripts/fetch_speciesfinder_database.sh /path/to/staging > /path/to/fetch.log 2>&1'
tail -n 5 /path/to/staging/status.tsv
```

Success is the `ready_for_preparation` status, which also records the exact
database-prefix and taxonomy-file paths. Prepare those paths as an inactive
candidate:

```bash
bactipipe prepare-speciesfinder-db \
  --database-prefix /path/to/staging/extracted/bacteria/bacteria_14_AT_110126 \
  --taxonomy-file /path/to/staging/extracted/bacteria/bacteria_14_AT_110126.tax \
  --archive-sha256-file /path/to/staging/downloads/speciesfinder_db.tar.gz.sha256 \
  --output-root /path/to/speciesfinder_db

bactipipe validate-identify-db \
  --database-root /path/to/speciesfinder_db/versions/VERSION \
  --allow-candidate
```

`bactipipe prepare-speciesfinder-db` publishes those verified components under
an immutable `versions/<version>/` directory. It records the archive checksum,
SpeciesFinder and KMA versions, hashes and sizes for every component, and builds
a SQLite accession-to-KMA-template lookup from the exact `.name` order. The
build also writes a row-level taxonomy/source audit. Explicit MAG, SAG,
metagenome-assembled, and single-cell templates cannot become callable ANI
references. Templates without a usable species name and numeric species taxid
remain in the upstream database and raw output but are excluded from BactiPipe
candidate nomination. The default status is candidate/inactive. `bactipipe validate-identify-db
--allow-candidate` validates this contract; production validation refuses it
unless a separately reviewed activation has set both production status and the
active flag.

For each sample, SpeciesFinder runs on filtered reads when they are available,
or on the assembly for the standalone fallback. BactiPipe parses the extended
TSV directly and verifies that SpeciesFinder did not mask an internal KMA
failure. Candidate order follows KMA query count when present, then score, and
is bounded before extraction. The current real KMA output has no read-count
column, so only the competitive mapping stage supplies mapped-read counts and
percentages. `kma seq2fasta` retrieves the exact nominated
templates; each extracted FASTA is checksummed and written into a sample-local
manifest used by skani and competitive mapping. This avoids building another
large custom Kraken database and prevents ANI from silently using references
that SpeciesFinder did not nominate.

Raw SpeciesFinder output may include placeholder taxa such as `Enterobacter
sp.` or `Enterobacteriaceae bacterium`. They remain visible as upstream context
but cannot consume a callable candidate slot, become an ANI identification, or
participate in mixed-isolate decisions. Only a reportable species label with a
numeric taxid can enter those steps.

The official taxonomy generally reports species rank. Explicit, taxid-based
rules retain *Streptococcus equi* subsp. *equi* (148942) and subsp.
*zooepidemicus* (40041) from template description/taxonomy evidence. These
remain different ANI/reporting units but share their NCBI species parent only
for mixed-signal assessment, preventing conserved reads from manufacturing a
false mixture.

## Previous reviewed-reference candidate contract

Database builds consume one reviewed TSV manifest. Every callable accession must
be an isolate-derived assembly with a current checksum. Each callable species
has exactly one canonical reference for tree display; additional diversity
references may improve ANI sensitivity but are collapsed to one result per
species. The explicit `include_in_kraken` field selects the reviewed subset
used for raw-read candidate nomination. Every callable species must have at
least one such reference, and a Kraken reference must itself be callable.
MAG/SAG references cannot be used for a positive diagnostic call.

The intended source hierarchy is:

1. a suitable type-material assembly;
2. an NCBI RefSeq reference genome;
3. an NCBI RefSeq representative genome; or
4. a reviewed isolate-derived GenBank fallback.

The broad backbone must include named non-veterinary species and near
neighbors, so an unusual isolate is not forced toward the closest organism in
the laboratory's controlled vocabulary. The veterinary layer identifies
important species, synonyms, deliberate subspecies, difficult complexes, and
species needing additional diversity representatives. GTDB supplies an
independent species-cluster/taxonomy check, while displayed names remain
reviewed NCBI/veterinary names.

The first inventory is generated from pinned GTDB R232 bacterial metadata with
`bactipipe plan-identify-inventory`. It scans the complete metadata stream but
keeps only two bounded groups: eligible assemblies belonging to the controlled
veterinary species catalog, and eligible GTDB species representatives from the
same NCBI or GTDB genera. Unnamed GTDB placeholder clusters are excluded from
this diagnostic-call database unless they resolve directly to a controlled
species. Common organisms therefore cannot contribute an unbounded number of
assemblies. The command writes an audit TSV, a deduplicated NCBI assembly-
accession list, an explicit missing-controlled-species TSV, and a JSON count
summary. It does not download genomes, create a reviewed manifest, or activate
a database.

Controlled membership is matched by NCBI taxid, never by name alone. Obsolete
taxids are resolved through the pinned NCBI `merged.dmp`. The normal unit is a
species-ranked taxid; direct `ncbi_taxid` matches preserve deliberate catalog
exceptions such as veterinary subspecies. Official taxid merges may make two
legacy catalog rows aliases of one current species, and that relationship is
retained in the audit rather than creating duplicate callable species. This
prevents an unnamed cluster under a named parent from being promoted merely
because its submitted text resembles a controlled name.

`bactipipe fetch-identify-reports` takes the proposed assembly-accession list
and retrieves current NCBI assembly reports in independently resumable batches.
Each completed batch is checksummed, missing or withdrawn accessions are
reported explicitly, current assembly-version substitutions are recorded in
`accession_resolution.tsv`, and a combined JSONL report is published
atomically. This metadata-only stage does not retrieve genomic FASTA files.

`bactipipe reconcile-identify-inventory` joins those current reports back to
the proposal. It promotes an available diversity reference when a canonical
accession disappeared and excludes suppressed assemblies, unnamed placeholder
taxa, failed NCBI taxonomy checks, and current taxid conflicts. Scientific-name
changes with the same resolved taxid are accepted as nomenclature updates. The
focused review and missing-species audits retain every rejected decision.

Controlled gaps are queried independently against current NCBI assemblies.
Only an exact controlled taxid with isolate provenance, an acceptable NCBI
taxonomy check, completeness >=90%, contamination <=5%, no more than 500
contigs, N50 >=20 kb, and total length >=100 kb can supplement the inventory.
The current pass recovered *Bacillus fonticola*, *Serratia odorifera*,
*Staphylococcus schleiferi*, and *Vibrio alfacsensis*. *Actinomyces pacaensis*,
*Clostridium carboxidivorans*, *Clostridium manihotivorum*, and *Pasteurella
stomatis* remain absent and auditable; the database build does not weaken its
gates to force their inclusion.

`bactipipe download-identify-refs` downloads the accepted assemblies into one
FASTA per assembly. Downloads are independently checkpointed in bounded
batches, failed batches split to isolate problematic accessions, every FASTA
is counted and checksummed, and its total length must equal current NCBI
metadata before a reviewed manifest can be written. Bacterial multi-contig
assemblies are preserved. Only canonical references enter the compact custom
Kraken2 candidate library; all reviewed references remain available to skani.

`bactipipe build-identify-db` copies or links exact FASTAs, verifies checksums,
builds a skani sketch and a custom Kraken2 index, records the pinned NCBI
taxonomy checksums and tool/database identity, and publishes an immutable
candidate under `versions/`. Kraken FASTA headers receive explicit manifest
taxids. References are imported in bounded temporary batches so the build does
not create another full-size concatenated FASTA; each successful batch is
deleted immediately, and Kraken build intermediates are cleaned after the
final index is verified. An optional hash-size ceiling exists only for
controlled sensitivity experiments. The command never activates its candidate.

## Phylogenetic placement

Unresolved, non-mixed assemblies may be placed with GTDB-Tk's bac120 marker
set. The wrapper runs placement with `--place_species`, finds the tree
containing the sample, keeps only the nearest named species, removes duplicate
strain representatives, and writes Newick, SVG, and neighbor TSV artifacts.
NCBI species names are preferred for display. A mixed assembly suppresses the
tree because concatenated markers from multiple organisms can give a
misleading placement. Full QC runs defer placement until run-wide
donor/recipient reconciliation is complete, so a weak secondary later
classified as possible carryover cannot incorrectly suppress an otherwise
eligible unresolved-sample tree. Selective sample updates repeat that decision
against the merged full cohort. Too few recovered markers or named neighbors
yields an explicit unavailable result.

GTDB-Tk uses disk-backed pplacer scratch by default. Controlled validation on
a host with sufficient free RAM may opt into in-memory placement with
`--gtdbtk-in-memory`; this is an explicit resource decision, not an automatic
fallback after a disk-space failure.

## Output contract

Each sample receives:

- `species_identification.json`
- `species_identification.tsv`
- `closest_species.tsv`
- the raw `skani.tsv`
- `speciesfinder_candidates.tsv` and the raw SpeciesFinder/KMA output when the
  SpeciesFinder database is selected
- the exact extracted candidate-reference manifest and competitive-mapping evidence
- optional Kraken2 evidence only for the previous custom candidate contract
- `cross_sample_signals.tsv` when run-wide donor/recipient evidence is available
- optional `phylogenetic_placement/closest_species_tree.nwk`
- optional `phylogenetic_placement/closest_species_tree.svg`
- optional `phylogenetic_placement/closest_species_tree_neighbors.tsv`

Taxonomic confirmation is evaluated at the resolution claimed by this
workflow. Historical *Salmonella enterica* serovar labels, including the older
`Salmonella typhimurium` form, are compared with *Salmonella enterica* at the
species level; a Pass does not constitute serovar confirmation. Clinically
important controlled subspecies remain separate. *Streptococcus equi* subsp.
*equi* and subsp. *zooepidemicus* therefore require distinct reported names and
evidence; neither is reduced to generic *S. equi* for confirmation.

Close ANI candidates are considered ambiguous unless the best reference also
has the configured minimum bidirectional aligned-fraction advantage (two
percentage points by default). This retains ambiguity when both genome identity
and breadth are effectively tied, while permitting a materially broader best
match to resolve closely related taxa. The initial complete-cohort validation
uses this specifically to retest *S. equi* subsp. *equi* versus subsp.
*zooepidemicus*; expected organism names never force the result.

Kraken candidate extraction follows the reviewed callable taxids, including a
clinically retained subspecies even when NCBI assigns it a rank below species.
In that mode direct taxon counts prevent a parent species row from
double-counting reads assigned to a retained child taxon.

Competitive read mapping is not used to distinguish clinically retained
subspecies that share one NCBI species parent. Their read signals are grouped
only for mixed-culture assessment because conserved reads otherwise map across
both references and falsely resemble a mixture. The units remain separate in
the manifest and ANI decision, where whole-assembly identity and bidirectional
aligned fraction determine the reported subspecies. Unrelated species remain
separate mixed-signal candidates.

The candidate3c regression completed on 2026-09-09. Its focused two-platform
gate confirmed all eight routine *S. equi* subspecies datasets exactly and did
not force the poor-quality challenge. All 52 complete-cohort cases executed;
comparison with candidate2 showed only the intended subspecies, historical
Salmonella-label, and related display-name changes. This is a software
regression pass, not diagnostic activation approval. Independently
characterized positives, negatives, and controlled mixtures on Illumina and
Nanopore remain required before thresholds are locked.

## Validation and activation gate

Before production activation:

1. **Completed 2026-09-09:** retrieve and checksum the official SpeciesFinder
   bacterial database; prepare and validate the immutable inactive version.
2. **Completed for the first assembly case:** run a format/extraction smoke
   test against the real KMA database and confirm every nominated accession
   resolves to the exact extracted template. The first case nominated
   *Kosakonia oryzae* and its exact `CP065358.1` template; skani measured
   98.61% ANI with 92.51% reference and 89.28% query aligned fraction.
   The subsequent read-backed smoke confirmed the same species at 98.92% ANI,
   93.74% reference and 88.51% query aligned fraction, with 272,886
   competitively mapped reads and no qualifying secondary signal.
3. Run the focused Candidate 4 cases, then the complete 52-case regression
   cohort, without changing the production KmerFinder path.
4. Benchmark independently characterized positive and negative veterinary isolates on Illumina and
   Nanopore.
5. Include low-coverage, contaminated, mixed-culture, and unresolved samples,
   including P2631472.
6. Review difficult species complexes and define explicit taxid-based
   group-level or specialist-test rules. Candidate2 validated the first rule
   for the *E. coli/Shigella* complex on the focused two-platform regression
   cohort; broader validation remains required.
7. Measure SpeciesFinder candidate recall and competitive-mapping false-signal
   behavior by platform, then approve separate Illumina and Nanopore confirmation gates.
8. Confirm that existing confident reports retain their column contract and
   that low-coverage KmerFinder hits are no longer positive identifications.
9. Validate the implemented run-wide donor/recipient carryover classification
   on complete multiplexed cohorts. It must retain possible carryover evidence,
   must not hide substantial secondary signals, and must not dismiss a signal
   when no credible primary donor exists. The first 12-case cohort gate passed
   on 2026-09-08: C04 PrimeLyse remained *Escherichia coli* and its weak
   *Salmonella enterica* evidence was retained as donor-linked possible
   cross-talk rather than a true mixed-isolate call.
10. Activate only by a separate reviewed maintenance step; retire the
   KmerFinder runtime database only after rollback no longer depends on it.

The first full-R232 placement challenge completed on 2026-09-08. P2631472
recovered 119 unique bac120 markers (98.48% of the alignment) and placed within
*Bibersteinia*, nearest to *Bibersteinia trehalosi*. GTDB-Tk classified the
assembly as a species-level novelty by RED, while direct skani comparison to
the retained *B. trehalosi* reference did not meet the search floor. The
pipeline therefore correctly keeps `no_confident_identification` and exposes
the genus-level placement/tree only as supporting evidence.

The next inactive validation stage uses 52 samples from complete Illumina and
Nanopore cohorts so that expected isolates, unresolved cases, a QC challenge,
and donor-linked secondary signals are evaluated with their real cohort
context. The organism entries in historical sample sheets are expectations,
not independently verified truth, so discordances require scientific review.
