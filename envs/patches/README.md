# KMA reassignment correction

Validated on the captured stall and four successful Salmonella genomes: all
four complete allele profiles and cgST summaries were identical before/after.
The affected sample completed cgMLST in 208 seconds with reassignment enabled.
The installer builds the correction from pinned source and records checksums.
The separate build-label patch identifies it as
`1.6.13+adl-reassign-progress-1`; it makes no algorithm change.

`kma-1.6.13-reassign-progress.patch` applies to KMA 1.6.13 source commit
`5cfd66b310437035d2f866d97dfba8b85dc7e408`. In `reassign_matchseqs`, a failed
search for the candidate's seed incorrectly resets the search offset to the
candidate length. A sufficiently long consensus can then revisit the same
offset indefinitely. The correction returns "no match" for that candidate
when no seed remains, allowing the outer reassignment candidate search to
continue. It does not disable `-reassign`, alter scoring, filter loci, or change
the cgMLST database.

Diagnostic evidence captured a 1,290-base consensus and 168-base candidate
`STMMW_11301_299`, repeatedly visiting offset 169. The candidate's first 32 bases
are absent from the captured consensus. The consensus's repetitive character
does not itself establish poor assembly quality or contamination.

The probe is linked against original and corrected KMA libraries to demonstrate
the original loop and test successful matches and rejected candidates. Validation
must also run cgMLST with reassignment enabled and compare complete allele
profiles/cgST on successful controls before this correction is promoted.

Do not describe a corrected binary as an unchanged upstream build. Retain the
upstream commit, patch identity (`adl-reassign-progress-1`), patch checksum and
binary checksum in execution provenance if deployed later.

Installation: clone the pinned KMA commit, then run
`bash envs/install_kma_reassign_fix.sh /path/to/kma /path/to/bactipipe-kma-env`.
Requires a C compiler, make, patch, zlib headers, git and Python 3. Web builds
compile in an isolated stage and copy only the binary and provenance into the
runtime environment. The SpeciesFinder KMA environment is not changed.

Known separate upstream issue: four synthetic tests with candidate lengths
divisible by 32 accepted mismatches in both original and corrected KMA. This
patch does not address that independent comparison defect. Do not claim all
matching behavior is validated or silently include a broader correction here.
