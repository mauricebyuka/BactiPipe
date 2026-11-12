# Conventional MLST scheme catalogue

Existing `relate` organism tokens keep their configured MLST, cgMLST and serotyping
choices. For an organism outside that registry, `--organism "Bacillus subtilis"`
looks up the exact species in the installed MLST catalogue. Case, whitespace and
underscores are normalized; spelling, species, subspecies, and genus boundaries
are not guessed. This lookup enables conventional MLST only, never a new cgMLST
scheme or serotyper. It does not perform species identification.

## Build an expanded installation

Use `bactipipe update-db mlst --all-schemes --seed ACTIVE_MLST_DB --output NEW_VERSION
--work BUILD_WORKSPACE`. Output and workspace must be separate paths; an existing
output is refused. The active directory is never changed by this command.
Run a large install detached with a durable log. Reuse the same workspace to
resume verified per-URL downloads and per-scheme receipts. Do not change a seed
or preferences during a build; use a new workspace for a new snapshot.

Discovery combines PubMLST's classic-scheme XML list (also containing Pasteur
schemes), current PubMLST API listings, and BactiPipe's established scheme URLs.
Only conventional schemes with ST profiles and 2–20 loci are eligible. cgMLST,
wgMLST and rMLST are excluded. Some provider schemes may be unreachable or need
authorization; these remain explicitly unavailable in the inventory. Discovery
errors are recorded, so `READY_FOR_REVIEW` does not mean universal coverage.
The legacy XML list is not maintained for newly added schemes; API discovery
supplements it. Other providers/new scheme families need explicit catalogue
support, not an assumption that every online typing scheme was discovered.

Downloads use public HTTPS endpoints only. PubMLST may restrict unauthenticated
data to submissions on or before 2024-12-31. Provider access notices are retained;
these snapshots must not be described as unrestricted/current complete datasets.
Credentials are not collected or inferred by this installer.

An optional `--preferences FILE.tsv` accepts `species` and `source` columns, where
source is the authoritative BIGSdb scheme URL. This sets an explicit preferred
scheme when several apply. Existing configured scheme choices remain unchanged.
No arbitrary first-match ranking is used. Ambiguity, an unavailable preferred
scheme, or absent catalogue produces no MLST call and an explicit selection
status; another scheme is not silently substituted.

## Outputs and activation

- `config` and per-scheme FASTA/profile files: CGE `mlst.py` compatible layout.
- `mlst_catalog.json`: exact species, candidates, preferences, source URLs,
  installed/unavailable state, access notices and file checksums.
- `species_to_mlst.tsv`: readable species-to-preferred-scheme table.
- `bactipipe_mlst_sources.json`: source provenance for installed schemes.
- Workspace `inventory.json`, `progress.json`, receipts, download cache and
  `status.json`: resumability and audit trail. Terminal states are
  `READY_FOR_REVIEW` or `FAILED`.

With `--seed`, established configured schemes are copied without changing their
allele/profile bytes. New and reused schemes are checked for profile/locus
consistency and required FASTA loci. A build missing an established scheme is
not released for review. Review the table, unresolved choices, access notices,
and real-tool validation before arranging a separate idle-time activation at
`BACTIPIPE_DB_DIR/mlst_db`. Do not replace active files while analyses use them.
The configured-only updater refuses to overwrite an expanded catalogue.

`bactipipe relate --list-mlst-schemes --db-dir DB_ROOT` lists catalogue names and
preferred schemes. In a web override manifest use `organism=Other` and the full
species in `other_organism`; no new worksheet column is needed.

`relate` records its selection and catalogue checksum in run configuration,
execution logs and `analysis_result.json`. It verifies selected automatic scheme
file checksums before typing. Resume cannot reuse STs after a changed automatic
selection/catalogue. Existing configured-only historical runs retain their
previous behavior. Scheme installation does not rerun any cohort or change its
saved results; reanalysis must be requested explicitly.
