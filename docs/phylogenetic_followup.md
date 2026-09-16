# Separate, iteration-aware phylogenetic placement

September 15: ANI-first GTDB confirmation permits named MAG/environmental
representatives under the same ANI/AF and competing-species gates as isolate
references. Unknown species labels remain non-reportable. Source categories
remain in the audit. The implementation fingerprint invalidates cached old
confirmation decisions; a newly confirmed sample bypasses marker/tree work.

`bactipipe place-assembly` runs an explicit supporting analysis **after** QC.
It neither regenerates QC reports nor changes the SpeciesFinder/ANI species call.
Use original assemblies, not the fragmented SpeciesFinder nomination FASTA.

```bash
bactipipe place-assembly \
  --assembly results/assemblies/genomes/Sample01.fasta \
  --cache-root placement/Sample01 \
  --reference-root /databases/gtdb/r232/data \
  --metadata /databases/gtdb/r232/species_metadata/bac120_metadata_r232.tsv.gz \
  --threads 8
```

This explicit CLI command is synchronous; an external scheduler can execute it
independently of QC. The web application uses the same engine on its own queue.
Keep a separate persistent cache for each sample. GTDB-Tk runs in the `gtdbtk`
Conda environment. The reference installation must be immutable and include its
verified publisher package `MD5SUM.txt` beside `data/`.

An identical assembly reuses verified results. If the assembly changed, GTDB
extracts and aligns markers again. An identical masked bacterial marker alignment
can reuse the compact tree while refreshing marker counts. Changed markers,
references, software or processing settings require new placement. Cached files
are checksum-verified before reuse. Failed/corrupt work is never treated as success.

The supporting `classify` stage uses actual predicted proteins with `--genes`,
so it does not reuse or recompute whole-genome ANI. BactiPipe's preceding ANI
confirmation remains independent. Tree distance is not percent identity, and a
nearby leaf does not confirm that species. Mixed/contaminated assemblies require
review before calling this low-level command; the web scheduler applies eligibility
after cohort reconciliation. This command alone does not infer that eligibility.

Output includes SVG/Newick trees, a nearest-species TSV, marker/provenance JSON
and cache receipts. Only compact results are retained; large temporary GTDB work
is removed. The command reports the output directory and whether reuse occurred.
The older explicitly opted-in QC GTDB flags remain backward compatible and
synchronous; omit them when using this separate follow-up command.

The new marker-only stage remains pending real-tool/scientific activation review.
It has not been enabled in the production web application.
