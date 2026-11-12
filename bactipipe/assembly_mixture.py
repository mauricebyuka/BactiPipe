"""Bounded, assembly-first corroboration of secondary isolate signals.

Competitive read assignments are a screening test, not evidence of two genomes.
Compare non-overlapping assembly windows against *all* nominated chromosomes,
retaining alternative alignments. Shared windows cannot support both organisms.
The provisional safeguards below intentionally leave uncertain cases in Review.
"""
from __future__ import annotations

import gzip
import csv
import hashlib
import json
import subprocess
from dataclasses import asdict, dataclass, replace
from pathlib import Path


@dataclass(frozen=True)
class AssemblyMixturePolicy:
    version: str = "assembly_independent_support_v1"
    window_bases: int = 10000
    minimum_window_bases: int = 2000
    minimum_reference_contig_bases: int = 500000
    minimum_identity: float = 95.0
    minimum_window_coverage: float = 90.0
    minimum_score_margin: float = 1.0
    minimum_exclusive_bases: int = 100000
    minimum_exclusive_reference_percent: float = 10.0
    reference_bins: int = 20
    minimum_distributed_bins: int = 5
    minimum_shared_percent: float = 90.0
    minimum_shared_windows: int = 10
    minimum_reference_pair_ani: float = 95.0
    minimum_reference_pair_af: float = 50.0


POLICY = AssemblyMixturePolicy()


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        header, sequence = None, []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:].strip(), []
            else:
                sequence.append(line.strip())
        if header is not None:
            yield header, "".join(sequence)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def interval_bases(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def assess_paf(path, *, targets, primary_taxid, policy=POLICY):
    """Summarize unique sequence, not mapping count, with explicit alternatives.

Target coordinates are unioned so repeats cannot inflate support. Bins use the
alignment midpoint: one long/shared element cannot occupy many bins by itself.
Only a near-full-window alignment can corroborate a genome. Short alternative
alignments still compete on matching bases per query base.
"""
    windows = {}
    sizes = {}
    for target, metadata in targets.items():
        sizes[metadata["taxid"]] = sizes.get(metadata["taxid"], 0) + metadata["length"]
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 12 or fields[5] not in targets:
                continue
            name, qlen, start, end, _, target, _, tstart, tend, matches, length = fields[:11]
            qlen, start, end, tstart, tend, matches, length = map(
                int, (qlen, start, end, tstart, tend, matches, length))
            if qlen < policy.minimum_window_bases or not length:
                continue
            coverage = 100 * (end - start) / qlen
            identity = 100 * matches / length
            score = 100 * matches / qlen
            taxid = targets[target]["taxid"]
            entry = (score, coverage, identity, target, tstart, tend, qlen)
            previous = windows.setdefault(name, {}).get(taxid)
            if previous is None or entry[0] > previous[0]:
                windows[name][taxid] = entry
    support = {taxid: {"aligned_windows": 0, "shared_with_primary_windows": 0,
                       "exclusive_windows": 0, "intervals": {}, "bins": set()}
               for taxid in sizes}
    for hits in windows.values():
        primary = hits.get(primary_taxid)
        for taxid, entry in hits.items():
            score, coverage, identity, target, start, end, _ = entry
            if coverage < policy.minimum_window_coverage or identity < policy.minimum_identity:
                continue
            row = support[taxid]
            row["aligned_windows"] += 1
            if primary and primary[1] >= policy.minimum_window_coverage and primary[2] >= policy.minimum_identity and primary[0] >= score - policy.minimum_score_margin:
                row["shared_with_primary_windows"] += 1
            competitors = [other[0] for key, other in hits.items() if key != taxid]
            if competitors and score - max(competitors) < policy.minimum_score_margin:
                continue
            row["exclusive_windows"] += 1
            row["intervals"].setdefault(target, []).append((start, end))
            row["bins"].add((target, min(policy.reference_bins - 1,
                int(((start + end) / 2) * policy.reference_bins / targets[target]["length"]))))
    result = {}
    for taxid, row in support.items():
        bases = sum(interval_bases(intervals) for intervals in row.pop("intervals").values())
        bins = len(row.pop("bins"))
        percent = 100 * bases / sizes[taxid]
        shared = 100 * row["shared_with_primary_windows"] / max(row["aligned_windows"], 1)
        independent = (bases >= policy.minimum_exclusive_bases
            and percent >= policy.minimum_exclusive_reference_percent
            and bins >= policy.minimum_distributed_bins)
        result[taxid] = dict(row, exclusive_bases=bases, exclusive_reference_percent=percent,
            distributed_bins=bins, shared_with_primary_percent=shared,
            reference_bases=sizes[taxid], independent_support=independent)
    return result


def investigate_assembly(*, assembly, catalog, mixed, output_directory, threads=1,
                         competitors=(), ani_hits=(), policy=POLICY):
    """Run only for flagged secondary signals; never download a new database.

Long, non-plasmid-labelled reference contigs are a conservative chromosome
proxy, not a claim that every accessory/mobile sequence has been annotated.
Unavailable tools/references are Review, never evidence of a clean isolate.
"""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    evidence = {"policy": asdict(policy), "status": "not_assessed", "species": {},
                "primary_taxid": mixed.primary.species_taxid if mixed.primary else "",
                "assembly_sha256": sha256(assembly), "references": []}
    if not mixed.primary or not mixed.qualifying_secondary:
        return evidence
    primary_taxid = mixed.primary.species_taxid
    selected = (mixed.primary, *mixed.qualifying_secondary)
    targets = {}
    reference_paths = {}
    reference = output / "candidate_chromosomes.fna"
    query = output / "assembly_windows.fna"
    try:
        # Include other nominees as competitors, not just the flagged pair.
        wanted = {signal.species_taxid for signal in selected}
        groups = catalog.mixed_signal_groups_by_taxid
        selected_groups = {groups[taxid] for taxid in wanted if groups.get(taxid)}
        # Respect the established taxid-based reporting groups (e.g. E. coli /
        # Shigella). A sibling label must not steal all primary support and hide
        # a genuinely separate Salmonella genome. This creates no new groups.
        candidate_taxids = wanted | {s.species_taxid for s in competitors
            if s.competitive_mapped_reads is not None
            and groups.get(s.species_taxid) not in selected_groups}
        evidence["existing_reporting_groups"] = {taxid: groups[taxid]
            for taxid in candidate_taxids if groups.get(taxid)}
        records = []
        for taxid in sorted(candidate_taxids):
            hits = [hit for hit in ani_hits if hit.species_taxid == taxid]
            best = max(hits, key=lambda hit: (hit.ani, max(hit.aligned_fraction_query,
                       hit.aligned_fraction_reference))) if hits else None
            records.append(catalog.resolve(best.reference_path) if best
                           else catalog.canonical_for_species(taxid))
        with reference.open("w") as destination:
            for record in records:
                if record is None:
                    continue
                source = catalog.fasta_for(record)
                reference_paths[record.species_taxid] = source
                evidence["references"].append({"taxid": record.species_taxid,
                    "accession": record.accession, "species": record.species_name,
                    "sha256": sha256(source)})
                for index, (header, sequence) in enumerate(fasta_records(source)):
                    if "plasmid" in header.casefold() or len(sequence) < policy.minimum_reference_contig_bases:
                        continue
                    target = f"t{record.species_taxid}_{index}"
                    targets[target] = {"taxid": record.species_taxid, "length": len(sequence)}
                    destination.write(f">{target}\n{sequence}\n")
        if wanted - {metadata["taxid"] for metadata in targets.values()}:
            raise ValueError("Candidate chromosome references are incomplete; mixture evidence remains unresolved.")
        with query.open("w") as destination:
            for contig, (_, sequence) in enumerate(fasta_records(assembly)):
                for offset in range(0, len(sequence), policy.window_bases):
                    window = sequence[offset:offset + policy.window_bases]
                    if len(window) >= policy.minimum_window_bases:
                        destination.write(f">q{contig}_{offset}\n{window}\n")
        paf = output / "assembly_candidates.paf"
        command = ["minimap2", "-x", "asm10", "-c", "--secondary=yes", "-N", "50",
                   "-p", "0.5", "-t", str(max(1, min(threads, 4))), str(reference), str(query)]
        evidence["command"] = command
        with paf.open("w") as stdout, (output / "alignment.log").open("w") as stderr:
            subprocess.run(command, stdout=stdout, stderr=stderr, check=True, timeout=600)
        evidence["species"] = assess_paf(paf, targets=targets, primary_taxid=primary_taxid, policy=policy)
        for record in records:
            if record and record.species_taxid in evidence["species"]:
                evidence["species"][record.species_taxid]["species_name"] = record.species_name
        # Pairwise ANI is descriptive reference QC, not an automatic taxon merge.
        comparisons = []
        for signal in mixed.qualifying_secondary:
            pair_path = output / f"reference_pair_{signal.species_taxid}.tsv"
            command = ["skani", "dist", "-r", str(reference_paths[primary_taxid]),
                       "-q", str(reference_paths[signal.species_taxid]),
                       "-t", str(max(1, min(threads, 4))), "-o", str(pair_path)]
            try:
                with (output / "reference_comparison.log").open("a") as log:
                    subprocess.run(command, stdout=log, stderr=log, check=True, timeout=120)
                with pair_path.open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                comparisons.append({"secondary_taxid": signal.species_taxid, "measurements": rows})
            except (OSError, subprocess.SubprocessError):
                comparisons.append({"secondary_taxid": signal.species_taxid, "status": "not_measurable"})
        evidence.update(status="assessed", reference_comparisons=comparisons)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        evidence.update(status="not_assessed", reason=f"{type(exc).__name__}: {exc}")
    finally:
        # These exact generated FASTAs are reproducible; retain compact evidence.
        reference.unlink(missing_ok=True)
        query.unlink(missing_ok=True)
    (output / "assembly_mixture.json").write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    return evidence


def apply_assembly_evidence(mixed, evidence):
    """Reapply saved evidence after cohort exclusions without reviving a donor hit."""
    if not evidence or evidence.get("policy", {}).get("version") != POLICY.version:
        return mixed
    mixed = replace(mixed, assembly_evidence=evidence)
    if (evidence.get("status") != "assessed" or not mixed.primary
            or mixed.primary.species_taxid != evidence.get("primary_taxid")):
        return mixed
    species = evidence.get("species", {})
    overlapping_references = set()
    for comparison in evidence.get("reference_comparisons", ()):
        for row in comparison.get("measurements", ()):
            try:
                if (float(row["ANI"]) >= POLICY.minimum_reference_pair_ani
                        and max(float(row["Align_fraction_ref"]), float(row["Align_fraction_query"]))
                        >= POLICY.minimum_reference_pair_af):
                    overlapping_references.add(comparison["secondary_taxid"])
            except (KeyError, TypeError, ValueError):
                continue
    remaining, shared = [], []
    for signal in mixed.qualifying_secondary:
        row = species.get(signal.species_taxid, {})
        # Similarity of the references alone cannot rule out a real mixture.
        if (row.get("aligned_windows", 0) >= POLICY.minimum_shared_windows
                and row.get("shared_with_primary_percent", 0) >= POLICY.minimum_shared_percent
                and not row.get("independent_support", False)
                and signal.species_taxid in overlapping_references):
            shared.append(signal.species_name)
        else:
            remaining.append(signal)
    if shared and not remaining:
        return replace(mixed, status="reference_ambiguity", qualifying_secondary=(),
            reason="Secondary reference matches mostly explain the same assembled sequence; independent second-genome support was not established. This is not a contamination verdict.")
    return replace(mixed, qualifying_secondary=tuple(remaining))


def independent_mixture_supported(mixed):
    evidence = mixed.assembly_evidence
    if (evidence.get("policy", {}).get("version") != POLICY.version
            or evidence.get("status") != "assessed" or not mixed.primary
            or mixed.primary.species_taxid != evidence.get("primary_taxid")):
        return False
    species = evidence.get("species", {})
    return (species.get(mixed.primary.species_taxid, {}).get("independent_support", False)
        and any(species.get(signal.species_taxid, {}).get("independent_support", False)
                for signal in mixed.qualifying_secondary))
