"""Reusable orchestration for one bacterial species-identification sample."""

from __future__ import annotations

import json
import hashlib
import os
import re
import shutil
import subprocess
import uuid
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Sequence

from bactipipe.kraken_evidence import (
    parse_kraken2_report,
    run_kraken2,
    write_taxonomic_signals,
)
from bactipipe.mixed_signal_mapping import confirm_signal_breadth
from bactipipe.phylogenetic_placement import run_bac120_placement
from bactipipe.species_database import validate_database
from bactipipe.speciesfinder_database import (
    DATABASE_KIND as SPECIESFINDER_DATABASE_KIND,
    extract_candidate_references,
    validate_speciesfinder_database,
)
from bactipipe.speciesfinder_evidence import (
    run_speciesfinder,
    select_candidate_hits,
    speciesfinder_hits_to_signals,
    write_speciesfinder_hits,
)
from bactipipe.speciesfinder_query import prepare_nomination_query
from bactipipe.__version__ import __version__
from bactipipe.species_identification import (
    IdentificationPolicy,
    IdentificationResult,
    ReferenceCatalog,
    MixedSignalAssessment,
    MixedSignalPolicy,
    assess_mixed_signals,
    default_kraken_confidence,
    default_mixed_signal_policy,
    decide_identification,
    parse_kmerfinder_output,
    parse_skani_output,
    parse_taxonomic_signals,
    run_skani_search,
    write_identification_outputs,
)


@dataclass(frozen=True)
class LegacySummaryFields:
    """Values used by the stable BactiPipe QC summary columns."""

    identified_organism: str
    match_percent: str
    taxonomic_confirmation: str
    coverage_organism: str | None


def command_version(command: Sequence[str]) -> str:
    """Return a report-friendly version for a tool that was actually invoked."""

    try:
        version_flag = "-v" if Path(command[0]).name == "kma" else "--version"
        completed = subprocess.run(
            [*command, version_flag],
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "v?"
    if completed.returncode != 0:
        return "v?"
    raw = (completed.stdout or completed.stderr or "").splitlines()
    match = re.search(r"\bv?(\d+(?:\.\d+){1,3})\b", raw[0] if raw else "")
    return f"v{match.group(1)}" if match else "v?"


def identification_database_version(database_root: str | Path) -> str:
    identity_path = Path(database_root) / "database_identity.json"
    try:
        identity = json.loads(identity_path.read_text())
    except (OSError, ValueError, TypeError):
        return "unknown"
    return str(identity.get("version") or "unknown")


def identification_database_kind(database_root: str | Path) -> str:
    """Return the declared backend without accepting an unverified database."""

    identity_path = Path(database_root) / "database_identity.json"
    try:
        identity = json.loads(identity_path.read_text())
    except (OSError, ValueError, TypeError):
        return ""
    return str(identity.get("database_kind") or "reviewed_reference_collection")


def legacy_summary_fields(result: IdentificationResult) -> LegacySummaryFields:
    """Adapt the richer result without inventing a numeric match value."""

    if result.status == "confirmed_species" or (
        result.status in {"suspected_culture_contamination", "mixed_signal_review"}
        and result.match_percent is not None
        and result.reported_organism != "No confident identification"
    ):
        top_hit = result.closest_species[0] if result.closest_species else None
        coverage_organism = (
            top_hit.species_name
            if top_hit is not None and top_hit.reporting_group
            else result.reported_organism
        )
        return LegacySummaryFields(
            identified_organism=result.reported_organism,
            match_percent=(
                f"{result.match_percent:.2f}"
                if result.match_percent is not None
                else "N/A"
            ),
            taxonomic_confirmation=("MixedSignalReview" if result.status == "mixed_signal_review"
                                    else result.taxonomic_confirmation),
            # Genome-size lookup remains species-specific even when the
            # diagnostic report deliberately uses a broader complex label.
            coverage_organism=coverage_organism,
        )
    # Run summaries report identification, not nomination/mixed-signal evidence.
    # Preserve the detailed analytical status and closest hits in sample sidecars.
    label = "No confident identification"
    coverage_organism = (
        result.closest_species[0].species_name if result.closest_species else None
    )
    return LegacySummaryFields(
        identified_organism=label,
        match_percent="N/A",
        taxonomic_confirmation=("MixedSignalReview" if result.status == "mixed_signal_review"
                                else result.taxonomic_confirmation),
        coverage_organism=coverage_organism,
    )


def identify_assembly(*, output_directory: str | Path, **kwargs) -> IdentificationResult:
    """Run into a sibling directory and atomically publish complete evidence."""

    final = Path(output_directory).resolve()
    final.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    staging = final.parent / f".{final.name}.building-{token}"
    backup = final.parent / f".{final.name}.backup-{token}"
    try:
        result = _identify_assembly_once(
            output_directory=staging,
            **kwargs,
        )
        if final.exists():
            os.replace(final, backup)
        try:
            os.replace(staging, final)
        except Exception:
            if backup.exists() and not final.exists():
                os.replace(backup, final)
            raise
        shutil.rmtree(backup, ignore_errors=True)
        return result
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _identify_assembly_once(
    *,
    sample_id: str,
    assembly: str | Path,
    database_root: str | Path,
    output_directory: str | Path,
    expected_organism: str = "Unknown",
    reads: Sequence[str | Path] = (),
    platform: str = "",
    kraken_database: str | Path | None = None,
    taxonomic_signals_tsv: str | Path | None = None,
    kmerfinder_results: str | Path | None = None,
    gtdb_metadata_tsv: str | Path | None = None,
    gtdbtk_command_prefix: Sequence[str] = ("gtdbtk",),
    gtdbtk_use_scratch: bool = True,
    threads: int = 1,
    kraken_confidence: float | None = None,
    policy: IdentificationPolicy = IdentificationPolicy(),
    mixed_signal_policy: MixedSignalPolicy | None = None,
    allow_candidate_database: bool = False,
    require_mixed_signal_assessment: bool = False,
    speciesfinder_executable: str = "speciesfinder",
    kma_executable: str = "kma",
    maximum_candidate_species: int = 20,
    candidate_references_per_species: int = 3,
) -> IdentificationResult:
    if kraken_database and taxonomic_signals_tsv:
        raise ValueError("Supply Kraken2 input or a taxonomic-signal TSV, not both")
    if mixed_signal_policy is None:
        mixed_signal_policy = default_mixed_signal_policy(platform)
    if kraken_confidence is None:
        kraken_confidence = default_kraken_confidence(platform)
    database_root = Path(database_root).resolve()
    database_kind = identification_database_kind(database_root)
    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)

    raw_skani = output / "skani.tsv"
    catalog: ReferenceCatalog | None
    speciesfinder_invoked = database_kind == SPECIESFINDER_DATABASE_KIND
    nomination_query = None
    if speciesfinder_invoked:
        if kraken_database:
            raise ValueError(
                "A SpeciesFinder candidate database cannot be combined with a Kraken2 database"
            )
        prepared_database = validate_speciesfinder_database(
            database_root, allow_candidate=allow_candidate_database
        )
        database_identity = prepared_database.identity
        speciesfinder_output = output / "speciesfinder"
        # Isolate nomination uses assembled sequence only, on both platforms.
        # Reads remain available for confirmation of assembly-nominated taxa;
        # do not introduce a FASTQ fallback or an extra read-based screen here.
        # Bound KMA query lengths to avoid origin-sensitive whole-contig
        # alignment failures. ANI, placement and downstream retain the original.
        query_directory = output / "speciesfinder_query"
        nomination_query = prepare_nomination_query(assembly, query_directory)
        query_path = query_directory / nomination_query["query_filename"]
        candidate_inputs = (query_path,)
        speciesfinder_hits = run_speciesfinder(
            inputs=candidate_inputs,
            output_directory=speciesfinder_output,
            database_prefix=prepared_database.database_prefix,
            taxonomy_file=prepared_database.taxonomy_file,
            executable=speciesfinder_executable,
        )
        # Only our derived query is disposable; keep its checksums/coordinates.
        query_path.unlink()
        write_speciesfinder_hits(
            speciesfinder_hits, output / "speciesfinder_candidates.tsv"
        )
        selected_speciesfinder_hits = select_candidate_hits(
            speciesfinder_hits,
            maximum_species=maximum_candidate_species,
            references_per_species=candidate_references_per_species,
        )
        candidate_root = output / "speciesfinder_candidate_references"
        references = extract_candidate_references(
            database=prepared_database,
            hits=selected_speciesfinder_hits,
            output_directory=candidate_root,
            maximum_species=maximum_candidate_species,
            references_per_species=candidate_references_per_species,
            kma_executable=kma_executable,
        )
        if references:
            catalog = ReferenceCatalog.from_tsv(
                candidate_root / "reference_manifest.tsv"
            )
            run_skani_candidate_dist(
                assembly=assembly,
                references=[record.fasta_path for record in references],
                output_path=raw_skani,
                threads=threads,
            )
            ani_hits = parse_skani_output(raw_skani, catalog)
        else:
            catalog = None
            raw_skani.write_text(
                "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\t"
                "Align_fraction_query\tRef_name\tQuery_name\n"
            )
            ani_hits = ()
        speciesfinder_signals = (
            speciesfinder_hits_to_signals(selected_speciesfinder_hits)
            if reads
            else ()
        )
        if speciesfinder_signals:
            write_taxonomic_signals(
                speciesfinder_signals,
                output / "taxonomic_signals.preliminary.tsv",
            )
    else:
        database_identity = validate_database(
            database_root, allow_candidate=allow_candidate_database
        )
        catalog = ReferenceCatalog.from_tsv(database_root / "reference_manifest.tsv")
        run_skani_search(
            assembly=assembly,
            database=database_root / "skani",
            output_path=raw_skani,
            threads=threads,
        )
        ani_hits = parse_skani_output(raw_skani, catalog)
        speciesfinder_signals = ()
    kmer = parse_kmerfinder_output(kmerfinder_results) if kmerfinder_results else ()

    if kraken_database:
        if not reads or not platform:
            raise ValueError("Kraken2 evidence requires reads and platform")
        kraken_output = output / "kraken2"
        kraken_report = run_kraken2(
            reads=reads,
            database=kraken_database,
            output_directory=kraken_output,
            sample_id=sample_id,
            threads=threads,
            confidence=kraken_confidence,
        )
        signals = parse_kraken2_report(
            kraken_report,
            identification_names_by_taxid=catalog.identification_names_by_taxid,
        )
        write_taxonomic_signals(
            signals, kraken_output / "taxonomic_signals.preliminary.tsv"
        )
    else:
        if taxonomic_signals_tsv:
            signals = parse_taxonomic_signals(taxonomic_signals_tsv)
        else:
            signals = speciesfinder_signals
    if signals and reads:
        if not platform:
            raise ValueError("Read-based breadth confirmation requires platform")
        if catalog is None:
            raise ValueError(
                "Read evidence was produced without extractable candidate references"
            )
        signals = confirm_signal_breadth(
            signals=signals,
            catalog=catalog,
            reads=reads,
            platform=platform,
            output_directory=output / "mixed_signal_mapping",
            threads=threads,
            preferred_species_taxids=(
                (ani_hits[0].species_taxid,) if ani_hits else ()
            ),
            maximum_species=(
                maximum_candidate_species if speciesfinder_invoked else 5
            ),
        )
    # Publish the final signal contract even when no candidate was nominated.
    # Batch reconciliation and unresolved-sample placement must not silently
    # skip a valid no-hit result merely because the table has zero data rows.
    write_taxonomic_signals(signals, output / "taxonomic_signals.tsv")

    reporting_groups = (
        catalog.mixed_signal_groups_by_taxid if catalog is not None else {}
    )
    reporting_group_members = (
        catalog.reporting_group_member_names if catalog is not None else {}
    )
    mixed_signal = assess_mixed_signals(
        signals,
        policy=mixed_signal_policy,
        reporting_groups=reporting_groups,
    )
    if mixed_signal.qualifying_secondary and catalog is not None:
        from .assembly_mixture import investigate_assembly, apply_assembly_evidence
        evidence = investigate_assembly(assembly=assembly, catalog=catalog,
            mixed=mixed_signal, competitors=signals, ani_hits=ani_hits,
            output_directory=output / "assembly_mixture", threads=threads)
        mixed_signal = apply_assembly_evidence(mixed_signal, evidence)
    if require_mixed_signal_assessment and mixed_signal.status == "not_assessed":
        mixed_signal = MixedSignalAssessment(
            status="indeterminate",
            primary=None,
            reason=(
                "No species-level read evidence was available to exclude a mixed "
                "bacterial signal."
            ),
        )

    result = decide_identification(
        sample_id=sample_id,
        ani_hits=ani_hits,
        expected_organism=expected_organism,
        kmerfinder_support=kmer,
        mixed_signal=mixed_signal,
        policy=policy,
        reporting_group_members=reporting_group_members,
    )
    result = replace(
        result,
        mixed_signal_thresholds=asdict(mixed_signal_policy),
    )
    write_identification_outputs(result, output)
    gtdbtk_invoked = False
    if gtdb_metadata_tsv:
        placement_output = output / "phylogenetic_placement"
        if result.status in {
            "mixed",
            "mixed_signal_review",
            "suspected_culture_contamination",
        }:
            placement_output.mkdir(parents=True, exist_ok=True)
            (placement_output / "phylogenetic_placement.json").write_text(
                json.dumps(
                    {
                        "status": "suppressed",
                        "reason": (
                            "A combined marker tree is not reliable when mixed species "
                            "signals are present."
                        ),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
        elif result.status != "confirmed_species":
            gtdbtk_invoked = True
            run_bac120_placement(
                assembly=assembly,
                sample_id=sample_id,
                output_directory=placement_output,
                gtdb_metadata_tsv=gtdb_metadata_tsv,
                cpus=threads,
                pplacer_cpus=max(1, min(threads, 4)),
                command_prefix=gtdbtk_command_prefix,
                use_pplacer_scratch=gtdbtk_use_scratch,
            )
    tools = {"skani": command_version(("skani",))}
    if speciesfinder_invoked:
        tools["speciesfinder"] = python_package_version("SpeciesFinder")
        tools["kma"] = command_version((kma_executable,))
    if kraken_database:
        tools["kraken2"] = command_version(("kraken2",))
    if gtdbtk_invoked:
        tools["gtdbtk"] = command_version(gtdbtk_command_prefix)
    assembly_path = Path(assembly).resolve()
    assembly_sha256 = _sha256(assembly_path)
    if nomination_query and assembly_sha256 != nomination_query["source_sha256"]:
        raise RuntimeError("Assembly changed during species identification; retry with stable inputs")
    provenance = {
        "schema_version": 4,
        "sample_id": sample_id,
        "bactipipe_version": __version__,
        "database": database_identity,
        "tools": tools,
        "parameters": {
            "platform": platform,
            "expected_organism": expected_organism,
            "kraken_confidence": kraken_confidence,
            "identification_thresholds": asdict(policy),
            "mixed_signal_thresholds": asdict(mixed_signal_policy),
            "assembly_mixture_policy": mixed_signal.assembly_evidence.get("policy", {}),
            "mixed_signal_assessment_required": require_mixed_signal_assessment,
            "kraken_role": "candidate_nomination_only",
            "candidate_provider": (
                "speciesfinder" if speciesfinder_invoked else "reviewed_reference_collection"
            ),
            "speciesfinder_input_kind": "fragmented_assembly_fasta" if speciesfinder_invoked else None,
            "speciesfinder_query": nomination_query,
            "maximum_candidate_species": maximum_candidate_species,
            "candidate_references_per_species": candidate_references_per_species,
        },
        "assembly": {
            "filename": assembly_path.name,
            "size_bytes": assembly_path.stat().st_size,
            "sha256": assembly_sha256,
        },
        "reads": [
            {
                "filename": path.name,
                "size_bytes": path.stat().st_size,
                "mtime_ns": path.stat().st_mtime_ns,
            }
            for value in reads
            for path in [Path(value).resolve()]
        ],
        "result_status": result.status,
        "contamination_evidence": asdict(result.contamination_evidence),
    }
    temporary_provenance = output / "analysis_provenance.json.tmp"
    temporary_provenance.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    os.replace(temporary_provenance, output / "analysis_provenance.json")
    return result


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def python_package_version(distribution: str) -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version

        return f"v{version(distribution)}"
    except PackageNotFoundError:
        return "v?"


def run_skani_candidate_dist(
    *,
    assembly: str | Path,
    references: Sequence[str | Path],
    output_path: str | Path,
    threads: int = 1,
    executable: str = "skani",
) -> None:
    """Calculate ANI only against SpeciesFinder-nominated exact templates."""

    assembly_path = Path(assembly).resolve()
    reference_paths = [Path(value).resolve() for value in references]
    if not assembly_path.is_file() or assembly_path.stat().st_size == 0:
        raise ValueError(f"Assembly is missing or empty: {assembly_path}")
    if not reference_paths:
        raise ValueError("At least one SpeciesFinder candidate reference is required")
    for reference in reference_paths:
        if not reference.is_file() or reference.stat().st_size == 0:
            raise ValueError(f"Candidate reference is missing or empty: {reference}")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    reference_list = output.parent / "speciesfinder_reference_paths.txt"
    reference_list.write_text("".join(f"{path}\n" for path in reference_paths))
    temporary = output.with_suffix(output.suffix + ".tmp")
    completed = subprocess.run(
        [
            executable,
            "dist",
            str(assembly_path),
            "--rl",
            str(reference_list),
            "-o",
            str(temporary),
            "-t",
            str(max(1, int(threads))),
        ],
        text=True,
        capture_output=True,
    )
    (output.parent / "skani.stderr.log").write_text(completed.stderr or "")
    if completed.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"skani candidate ANI failed with exit {completed.returncode}; see "
            f"{output.parent / 'skani.stderr.log'}"
        )
    if not temporary.is_file() or temporary.stat().st_size == 0:
        raise RuntimeError("skani candidate ANI completed without a result table")
    os.replace(temporary, output)
