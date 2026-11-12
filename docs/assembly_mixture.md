# Assembly-first secondary-signal review

## Purpose and scope

Policy `assembly_independent_support_v1` replaces the Nanopore exception that
allowed a strong competitive read-mapping signal alone to imply contamination.
Both Illumina and Nanopore now require independent assembly evidence. Existing
species ANI/AF thresholds, clinically relevant subspecies, reporting groups,
cohort cross-talk rules and measured QC cutoffs are unchanged.

This is a conservative software policy, **not a validated diagnostic limit of
detection**. It does not establish whether contamination arose in culture,
extraction, library preparation or sequencing. The historical internal identifier
`suspected_culture_contamination` is retained for compatibility; the UI describes
supported contamination without asserting its source.

## Evidence sequence

1. SpeciesFinder nominates from assembly fragments. Whole-assembly ANI evaluates
   species identity. Competitive filtered-read mapping screens secondary taxa.
   Mapping percentages are assignments among nominated references, not organism
   abundance.
2. Only a flagged secondary triggers the bounded assembly check. Reuse extracted
   candidate references; use the best whole-assembly ANI reference per taxon,
   falling back to the existing canonical reference when no ANI is measurable.
   Other competitively assessed nominees remain competitors.
3. Divide a derived assembly copy into nonoverlapping 10 kb windows; retain
   tails/contigs at least 2 kb for this check. The original assembly is unchanged.
   Align against the candidate chromosome proxies with minimap2 asm10, retaining
   secondary alignments. Shorter fragments remain in the original QC assembly
   and species nomination, but do not corroborate contamination.
4. Chromosome proxies are reference contigs at least 500 kb whose headers do not
   identify plasmids. This deliberately conservative proxy is not a complete
   mobile-element annotation. Organisms with smaller/fragmented references may
   therefore remain Review.
5. Count a window only at >=95% alignment identity and >=90% window coverage.
   Distinct assignment requires a >=1 percentage-point advantage in matching
   bases per window base over every competing taxon. Near-ties support neither
   genome independently. Union target intervals to avoid repeated-locus inflation.
6. A genome has independent support only with >=100 kb distinct sequence,
   >=10% of its eligible reference sequence, and alignment midpoints distributed
   across >=5 of 20 reference bins. Both primary and at least one retained
   secondary must pass. These safeguards reject isolated conserved/mobile loci;
   they cannot guarantee exclusion of every possible horizontally transferred region.
7. Direct reference-versus-reference skani ANI is retained separately. A secondary
   is reference ambiguity when >=90% of at least 10 aligned windows are also
   adequately explained by the primary, independent secondary support is absent,
   and the two references reach >=95% ANI and either AF >=50%.
   Similar reference labels alone never merge species, rename a result, or dismiss
   independent evidence of a mixture.
8. Cohort review excludes qualifying carryover before the final decision. Reapply
   the saved assembly assessment to the retained signals; never restore an
   excluded secondary or reuse it after the primary taxon changes.

## Decisions and reports

| Evidence | Interpretation | Isolate QC |
| --- | --- | --- |
| Shared assembly + overlapping references; no independent secondary support | Reference ambiguity | No mixture-related failure |
| Weak signal meeting donor/recipient cross-talk gates | Possible cross-talk | No mixture-related failure |
| Secondary signal without sufficient independent evidence | Secondary-signal review | Review |
| Strong residual screening signal + distributed independent support for both genomes | Contamination supported | Fail |

A confident primary identity and its ANI are retained in Review and Fail cases.
If identity itself is ambiguous or unsupported, it remains unconfirmed.
Low depth, low quality and CheckM failures remain independent: resolving reference
ambiguity does **not** waive them. Missing evidence never becomes a measured Pass.

The organism-input warning is independent of secondary-signal warnings, so a typo
and a contamination/review indicator can both be shown. The sample page presents
distinct assembly support and explains that mapping percentages are not abundance.

## Resources, audit and limitations

The check uses at most four threads per flagged sample; there is no extra
reference download, read remapping, new database, or investigation for clean
single-signal samples. Each minimap2 call is bounded to 600 seconds and each
reference ANI comparison to 120 seconds. Tool/reference failures leave the
secondary assessment unresolved, not clean and not an automatic contamination call.

Per-sample `assembly_mixture/` retains the PAF, reference-pair ANI tables, logs,
and `assembly_mixture.json`. The JSON records policy, assembly/reference hashes,
reference names/accessions and counts. The result embeds this evidence under
`mixed_signal.assembly_evidence` so cohort reconciliation preserves it.
Generated query/reference FASTA copies are removed, not the source assemblies.

No new discriminating-site/pileup analysis is automatically launched. Cases that
still need read-level allele discrimination remain Review. Low-level mixtures,
very closely related strain mixtures and mixtures not reconstructed in the
assembly can be missed by these assembly gates; do not interpret their absence as
proof of purity. CheckM remains complementary evidence and its QC limit still
applies.

This change does not edit historical published outputs or launch reanalyses.
New executions use the new evidence policy. A source deployment is required
before production executions can use it.

## Focused saved-data validation (2026-11-04)

Read-only S3 staging and four-thread, network-disabled diagnostics used unchanged
assemblies and saved mapping results, without rerunning QC.

| Run / sample | Previous signal verdict | Revised signal verdict |
| --- | --- | --- |
| 261103_17_SEQ10_LMO / C20_P | Suspected contamination | Reference ambiguity; K. quasipneumoniae ANI 99.45% retained |
| 261103_17_SEQ10_LMO / C20_Z | Confirmed species | Unchanged |
| 251015_16_SEQ10_LMO / P2539688-4 | Suspected contamination | Supported contamination retained; primary Salmonella retained |

C20_P: CP077279.1 and AP023264.1 are 99.08% ANI, with AF 89.65%/89.13%.
366/378 secondary-aligned windows (96.83%) are also explained by the primary.
Secondary-exclusive sequence is 129,185 bp, only 2.45% of its reference, versus
1,007,341 bp (19.21%) for the primary. This is not evidence sufficient for two
independent genomes under this policy. The references' differing species labels
remain unchanged; no taxonomic misannotation is assumed proven.

P2539688-4: Salmonella has 4,396,389 bp (92.24%) of distinct support and
Bordetella hinzii has 3,615,785 bp (71.66%), each distributed across 20 bins.
That genuine positive case is not cleared by the new conservative rule.

Artifacts: `/home/ngslab/temp_data/assembly-mixture-20261104/`;
source staging script and execution log:
`/home/ngslab/status/assembly-mixture-20261104-stage.py`,
`/home/ngslab/status/assembly-mixture-20261104.log`.
The final `comparison.json` and `results/*/*/decision.json` contain the reviewed
decisions; the first execution log predates the final saved-evidence replay.

The real-data check above is Nanopore-focused. Synthetic regression tests exercise
both platforms, shared reference windows, true independent support, repeated
loci, sparse accessory sequence, incomplete alignments, missing evidence, cohort
exclusions and preservation of independent QC failures. Broader prospective
Illumina and defined-mixture validation remain necessary before treating the
provisional gates as laboratory-validated cutoffs.
