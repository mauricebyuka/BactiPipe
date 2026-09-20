# BactiPipe  

**BactiPipe** was developed and is actively used by the **NGS section of the Penn State Animal Diagnostic Laboratory** to bring **consistency** to the analysis of **BACTERIAL** whole-genome sequencing (WGS) data. Multiple tools and parameters are available for genomic data analysis but this flexibility often leads to variability in results from run to run. In diagnostic laboratories, such variability may compromize one of the most critical aspects of quality control in a diagnostic setting: **reproducibility and standardization**. BactiPipe addresses this challenge by wrapping together the core tools used in our laboratory with carefully chosen and fixed parameters, ensuring consistent analysis across all sequencing runs.  

BactiPipe provides **end-to-end quality control and genome assembly workflows** for both **Illumina** and **Oxford Nanopore** sequencing data. It automates read QC, genome assembly, organism identification, and assembly quality checks using tools such as **fastp, filtlong, Spades, Unicycler, Flye, skani, SpeciesFinder/KMA,** and **CheckM**. KmerFinder remains the legacy production identifier while the new workflow is validated. The pipeline supports data from **local directories, network-mounted storage, and Amazon S3**, making it suitable for both on-premises and cloud-based environments. Outputs include detailed logs, quality metrics plots, assembled genomes, and comprehensive QC summary reports.

For a detailed description of the pipeline, workflows, and usage examples, please see the [**BactiPipe Wiki Page**](https://github.com/ADL-genomics/BactiPipe/wiki).  

The user-facing CLI has five analysis commands:

- `bactipipe qc-illumina` and `bactipipe qc-nanopore` perform read QC, assembly,
  identification, and assembly QC.
- `bactipipe relate` performs organism-appropriate typing plus optional ANI,
  SKA, and cgMLST analysis on a related cohort of assembled genomes.
- `bactipipe detect` reports AMR determinants and virulence genes from assembled
  genomes.
- `bactipipe identify` runs the evidence-based species-identification
  workflow against an explicitly selected database version.

Separate placement and figure redrawing are options of `identify`:
`bactipipe identify --place --help` and `bactipipe identify --render --help`.
These work with local assemblies and retained evidence; no web application is required.

Prepare a verified SpeciesFinder installation with
`bactipipe update-db speciesfinder --help`. The existing required index, taxonomy,
archive checksum and output-root options are unchanged. This does not download
another large database or activate a candidate automatically. Validate with
`bactipipe update-db speciesfinder --verify --database-root PATH` (add
`--allow-candidate` for an inactive installation). `update-db all` continues to
update only the legacy tool-managed databases, not immutable SpeciesFinder versions.

Long historical database-planning/build commands and `place-assembly`,
`render-placement`, `prepare-speciesfinder-db`, `validate-identify-db` remain
accepted compatibility aliases, but are no longer advertised in ordinary help.
The older custom-reference inventory builders are maintenance/experimental tools,
not steps required to run QC, relate, detect or identify on a workstation.

`bactipipe detect --organism-list` prints one simple organism name per line, with
no tool/database table. For example, use `--organism Salmonella` or
`--organism "Enterococcus faecium"`. Species distinctions are retained when
needed, and older underscore names remain accepted. Tool selection is internal.
Only an exact supported AMRFinderPlus organism (ignoring case and
space/underscore differences, plus reviewed legacy group aliases) is passed to its `-O` option. Other organisms use
ABRicate `--db ncbi` without an organism option. Enterococcus uses VirulenceFinder
`virulence_ent`; faecium/lactis additionally use `virulence_entfm_entls`. Other
VirulenceFinder-supported groups retain their organism-specific databases;
unsupported groups use ABRicate `--db vfdb`. The existing Escherichia combined
VirulenceFinder/VFDB/ecoli_vf policy is preserved. ResFinder/CARD remain opt-in.
AMR-only and virulence-only modes remain independent. ABRicate/NCBI is a gene
screen, not equivalent to AMRFinderPlus organism-specific mutation screening;
outputs record the actual tools/databases. A failed tool is not a negative result.

The new identification workflow is not yet the production QC default. The
current candidate uses only the assembly FASTA for SpeciesFinder/KMA to nominate a bounded set of taxa and
exact database templates, then uses whole-genome ANI/aligned fraction for the
actual identification. There is no FASTQ nomination fallback or additional read-based discovery screen.
For nomination only, contigs are divided into non-overlapping pieces of at most
20,000 bp on both platforms, retaining short contigs and every final tail.
ANI, assembly metrics and downstream typing still use the original assembly.
The derived FASTA is removed after nomination; checksums and fragment coordinates
remain in `species_identification/<sample>/speciesfinder_query/`.
Read-level support for assembly-nominated taxa is confirmed by competitive mapping,
and run-wide donor/recipient comparison distinguishes weak possible multiplex
carryover from secondary signals requiring review or suspected culture contamination.
The evidence remains available. The workflow can place unresolved non-mixed
assemblies with GTDB bac120 markers. The older custom Kraken candidate remains
available only for comparison and rollback. See
[`docs/species_identification.md`](docs/species_identification.md).

Run `bactipipe <command> --help` for the current input contract and options.

For either QC command, repeat `--only-sample SAMPLE_ID` to select samples and
add `--update-existing` to safely replace only those samples in an existing
run output. BactiPipe stages the reanalysis separately and publishes it with a
directory swap; failed reanalysis leaves the existing output unchanged.

For controlled Nanopore reanalysis, `--input-run-name SOURCE_RUN` reuses that
physical input directory while `--run_name` continues to name the independent
output. This is primarily used by the web application's Admin-only benchmark
workflow.

## Primary-sheet downstream requests and report ordering

The web integration can derive downstream cohorts from the primary QC workbook's
`R/D*` column: `D` requests detection, and `R1`–`R4` request relatedness using
reference-library names in F5:F8 (labels in E5:E8). A blank request means QC only.
Multiple reference names, separated by commas, semicolons or newlines, share one
cohort with separate SNP-distance/ANI column pairs. Combined requests such as `R1,D` are supported. An explicit
downstream manifest takes precedence over these automatic requests.

The shared adapter is `bactipipe.primary_downstream.primary_downstream_tsv`.
Standalone QC commands do not automatically schedule downstream work; `relate`
and `detect` remain separate CLI commands consuming assemblies and command sheets.
Both commands sort report samples by supplied display name, followed by unnamed
samples in sheet order; `--sample-order sheet` preserves the entire input order.
Reference rows remain first. `relate --reference-specimen TEXT` sets the reference
row's Specimen label; the web application supplies its reference-library Display name.
Repeat `--reference`, `--reference-id` and `--reference-specimen` in matching order
for multiple references. Typing runs once per genome and SKA once per cohort;
ANI is calculated against each reference. All references remain at the top of the
consolidated table. PDFs display two reference pairs per panel, continuing onto
another page when needed; the TSV keeps every pair side by side. Multi-reference
`analysis_result.json` uses `bactipipe-relate-v2` with explicit reference IDs and
evidence paths. Single-reference columns and invocations remain compatible.

## Controlled organism sizes

The QC sample sheet accepts the organism names in
`bactipipe/data/pathogenic_bacteria.txt` plus `Lambda` and `Unknown`. Matching
ignores capitalization and extra whitespace, but BactiPipe does not guess at
unsupported taxonomic names. Unrecognized/blank entries use Unknown mode and
retain an input warning; they do not halt assembly.

Each controlled organism is resolved once and the exact taxid, reference
accession, genome size, catalog version, and catalog checksum used by the run
are recorded in `sample_execution_plan.tsv`. Illumina assemblers do not receive
a genome-size argument. Final reported depth on both platforms is QC-retained
read bases divided by actual assembly length, regardless of identification.
Unknown Nanopore assemblies use the bounded adaptive assembly workflow.

### Illumina filtering and depth

Raw mean quality is the base-weighted arithmetic mean of Phred scores across
both mates combined. Filtering runs only below `--minqual + 1` (default: below
Q29); exactly Q29 skips filtering. When needed, fastp applies a per-read mean
Q gate of `--minqual` (rounded upward if fractional), its normal adapter/quality
filters, and no minimum length. Empty output or malformed/unbalanced pairs fail
the sample. All subsequent QC, assembly and read evidence use this selected
dataset, which must meet the configured sample mean Q (default Q28).

Depth is selected read bases / total assembled FASTA length, not mapped depth
or a catalog-size estimate. No estimated-depth gate prevents assembly. Required
depth is 30X at mean Q≥30, 40X at 29≤Q<30, and 50X below Q29; confidently
identified E. coli requires another 10X. Expected and closest-reference names
do not add this requirement. Cohort review and selective retries update the
requirement when the identification changes. Decisions use unrounded values.

`qc_out/<sample>/illumina_read_selection.json` retains raw/selected counts,
full-precision quality, filtering command/version and input signatures.
`depth_metrics.json` records assembly-based depth and its threshold. Filtered
FASTQs are temporary; compact fastp reports/logs remain. Raw FASTQs are untouched.

### Illumina sample parallelism

SKESA and SPAdes now run read QC/assembly for up to three samples concurrently.
`--threads` is the total run CPU budget, capped by affinity/container quota;
`--sample-threads` defaults to 20 per sample and `--sample-jobs` defaults to 3.
Available CPU and RAM can reduce concurrency automatically. An optional
`--memory-gb` further limits the preparation budget. Use `--sample-jobs 1
--sample-threads 60 --threads 60` for a serial 60-thread comparison.

Identification remains serial after preparation; CheckM, cohort checks and final
reporting remain batch-wide. Failed samples do not stop successful samples.
`illumina_resources.json` and sample worker JSON/log files record allocations
and outcomes. Assembly cache identity includes thread/memory settings.
SPAdes receives a memory limit; SKESA's memory flag only applies in its
sorted-counter mode, so its normal mode still requires peak-memory monitoring.
Unicycler stays serial; Nanopore/downstream behavior is unchanged.

The installed application version is available through either alias:

```bash
bactipipe -v
bactipipe --version
```

BactiPipe follows the repository's [versioning policy](VERSIONING.md), with
release notes maintained in the [changelog](CHANGELOG.md).

For the `relate`/`detect` helper environment, create the compatible base and
then install the exact legacy CGE/SeqSero wrappers without their obsolete
dependency metadata:

```bash
mamba env create -f envs/bactipipe-kma.yml
mamba env create -f envs/genepid.yml
bash envs/install_genepid_tools.sh genepid
mamba env create -f envs/viramr.yml
mamba env create -f envs/kleborate.yml
bash envs/install_kleborate_compat.sh kleborate
```

Install or refresh the databases used by downstream typing and detection with:

```bash
bactipipe update-db mlst serotypefinder virulencefinder cgmlstfinder amrfinder abricate
```

The former CGE `mlst_db` Bitbucket repository has been removed. BactiPipe therefore rebuilds the
same directory format expected by CGE `mlst.py` from the official PubMLST and Institut Pasteur BIGSdb
REST APIs. Source endpoints, dates, loci, and server access notices are saved in
`mlst_db/bactipipe_mlst_sources.json`. If authenticated BIGSdb access has been issued to you, set the
complete HTTP Authorization header in the host-specific `BACTIPIPE_PUBMLST_AUTHORIZATION` or
`BACTIPIPE_PASTEUR_AUTHORIZATION` variable so credentials are sent only to the issuing service.
