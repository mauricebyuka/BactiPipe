"""Run-wide review of weak bacterial signals in multiplexed cohorts."""

from __future__ import annotations

import csv
import json
import os
import shutil
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence
from bactipipe.lambda_control import is_lambda_control

from bactipipe.phylogenetic_placement import run_bac120_placement

from bactipipe.species_identification import (
    CrossSampleEvidence,
    IdentificationPolicy,
    IdentificationResult,
    KmerSpeciesEvidence,
    MixedSignalPolicy,
    ReferenceCatalog,
    TaxonomicSignal,
    assess_mixed_signals,
    checkm_contamination_evidence,
    decide_identification,
    default_mixed_signal_policy,
    parse_skani_output,
    parse_taxonomic_signals,
    write_identification_outputs,
)
from bactipipe.species_workflow import legacy_summary_fields


@dataclass(frozen=True)
class CrossSamplePolicy:
    """Provisional two-sided carryover gates pending cohort validation."""

    maximum_source_fraction: float = 0.10
    maximum_sample_fraction: float = 0.20
    minimum_source_mapped_reads: int = 100

    def __post_init__(self) -> None:
        for name in ("maximum_source_fraction", "maximum_sample_fraction"):
            value = getattr(self, name)
            if not 0 <= value <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.minimum_source_mapped_reads < 1:
            raise ValueError("minimum_source_mapped_reads must be at least 1")


@dataclass(frozen=True)
class _GroupedSignal:
    sample_id: str
    key: str
    species_taxid: str
    species_name: str
    reporting_group: str
    mapped_reads: int


def _signal_key(signal: TaxonomicSignal, reporting_groups: Mapping[str, str]) -> str:
    group = reporting_groups.get(signal.species_taxid, "").strip()
    return f"group:{group.casefold()}" if group else f"taxid:{signal.species_taxid}"


def _group_signals(
    sample_id: str,
    signals: Sequence[TaxonomicSignal],
    signal_groups: Mapping[str, str],
    reporting_groups: Mapping[str, str],
) -> list[_GroupedSignal]:
    grouped: dict[str, list[TaxonomicSignal]] = {}
    for signal in signals:
        if signal.competitive_mapped_reads is None:
            continue
        grouped.setdefault(_signal_key(signal, signal_groups), []).append(signal)
    output = []
    for key, members in grouped.items():
        representative = max(
            members,
            key=lambda item: (
                item.confirmed_read_count,
                item.confirmed_read_percent,
            ),
        )
        output.append(
            _GroupedSignal(
                sample_id=sample_id,
                key=key,
                species_taxid=representative.species_taxid,
                species_name=representative.species_name,
                reporting_group=reporting_groups.get(
                    representative.species_taxid, ""
                ),
                mapped_reads=sum(item.confirmed_read_count for item in members),
            )
        )
    return sorted(output, key=lambda item: item.mapped_reads, reverse=True)


def classify_cross_sample_carryover(
    signals_by_sample: Mapping[str, Sequence[TaxonomicSignal]],
    *,
    reporting_groups: Mapping[str, str] | None = None,
    signal_groups: Mapping[str, str] | None = None,
    policy: CrossSamplePolicy = CrossSamplePolicy(),
) -> dict[str, tuple[CrossSampleEvidence, ...]]:
    """Classify low signals only when a credible multiplex donor is present."""

    reporting_groups = reporting_groups or {}
    signal_groups = signal_groups or reporting_groups
    grouped_by_sample = {
        sample_id: _group_signals(
            sample_id, signals, signal_groups, reporting_groups
        )
        for sample_id, signals in signals_by_sample.items()
    }
    primary_by_sample = {
        sample_id: rows[0] for sample_id, rows in grouped_by_sample.items() if rows
    }
    by_key: dict[str, list[_GroupedSignal]] = {}
    for rows in grouped_by_sample.values():
        for row in rows:
            by_key.setdefault(row.key, []).append(row)

    evidence_by_sample: dict[str, list[CrossSampleEvidence]] = {
        sample_id: [] for sample_id in signals_by_sample
    }
    for rows in by_key.values():
        primary_sources = [
            row
            for row in rows
            if primary_by_sample.get(row.sample_id) is not None
            and primary_by_sample[row.sample_id].key == row.key
        ]
        source = max(primary_sources, key=lambda item: item.mapped_reads) if primary_sources else None
        credible_source = bool(
            source is not None
            and source.mapped_reads >= policy.minimum_source_mapped_reads
        )
        for row in rows:
            sample_primary = primary_by_sample.get(row.sample_id)
            source_fraction = (
                row.mapped_reads / source.mapped_reads
                if source is not None and source.mapped_reads
                else 0.0
            )
            sample_fraction = (
                row.mapped_reads / sample_primary.mapped_reads
                if sample_primary and sample_primary.mapped_reads
                else 0.0
            )
            status = "observed"
            if source is not None and row == source and credible_source:
                status = "source"
            elif (
                credible_source
                and source is not None
                and row.sample_id != source.sample_id
                and source_fraction <= policy.maximum_source_fraction
                and sample_fraction <= policy.maximum_sample_fraction
            ):
                status = "possible_cross_talk"
            evidence_by_sample[row.sample_id].append(
                CrossSampleEvidence(
                    species_taxid=row.species_taxid,
                    species_name=row.species_name,
                    reporting_group=row.reporting_group,
                    status=status,
                    source_sample=(source.sample_id if credible_source and source else ""),
                    mapped_reads=row.mapped_reads,
                    source_mapped_reads=(
                        source.mapped_reads if credible_source and source else 0
                    ),
                    fraction_of_source=source_fraction if credible_source else 0.0,
                    fraction_of_sample_primary=sample_fraction,
                )
            )
    return {
        sample_id: tuple(
            sorted(rows, key=lambda item: item.mapped_reads, reverse=True)
        )
        for sample_id, rows in evidence_by_sample.items()
    }


def reconcile_batch_identifications(
    *,
    species_output_root: str | Path,
    database_root: str | Path,
    expected_organisms: Mapping[str, str],
    platform: str,
    identification_policy: IdentificationPolicy = IdentificationPolicy(),
    mixed_signal_policy: MixedSignalPolicy | None = None,
    cross_sample_policy: CrossSamplePolicy = CrossSamplePolicy(),
    assemblies: Mapping[str, str | Path] | None = None,
    gtdb_metadata_tsv: str | Path | None = None,
    gtdbtk_command_prefix: Sequence[str] = ("gtdbtk",),
    placement_threads: int = 1,
    checkm_contamination_by_sample: Mapping[str, object] | None = None,
    maximum_checkm_contamination: float = 10.0,
) -> dict[str, IdentificationResult]:
    """Apply cohort context and republish per-sample decisions atomically."""

    root = Path(species_output_root)
    global_manifest = Path(database_root) / "reference_manifest.tsv"
    global_catalog = (
        ReferenceCatalog.from_tsv(global_manifest)
        if global_manifest.is_file()
        else None
    )
    mixed_signal_policy = mixed_signal_policy or default_mixed_signal_policy(platform)
    checkm_contamination_by_sample = checkm_contamination_by_sample or {}
    signals_by_sample: dict[str, tuple[TaxonomicSignal, ...]] = {}
    sample_directories: dict[str, Path] = {}
    catalogs_by_sample: dict[str, ReferenceCatalog] = {}
    for sample_dir in sorted(root.iterdir() if root.is_dir() else ()):
        if is_lambda_control(expected_organisms.get(sample_dir.name)):
            continue  # Ignore even stale bacterial identification outputs for a control.
        if not sample_dir.is_dir():
            continue
        signals_path = sample_dir / "taxonomic_signals.tsv"
        skani_path = sample_dir / "skani.tsv"
        result_path = sample_dir / "species_identification.json"
        candidate_manifest = (
            sample_dir / "speciesfinder_candidate_references" / "reference_manifest.tsv"
        )
        sample_catalog = global_catalog
        if sample_catalog is None and candidate_manifest.is_file():
            try:
                sample_catalog = ReferenceCatalog.from_tsv(candidate_manifest)
            except ValueError as exc:
                if "contains no records" in str(exc):
                    sample_catalog = ReferenceCatalog(
                        (), root=candidate_manifest.parent
                    )
                else:
                    raise
        if (
            signals_path.is_file()
            and skani_path.is_file()
            and result_path.is_file()
            and sample_catalog is not None
        ):
            signals_by_sample[sample_dir.name] = parse_taxonomic_signals(signals_path)
            sample_directories[sample_dir.name] = sample_dir
            catalogs_by_sample[sample_dir.name] = sample_catalog
    reporting_groups, signal_groups, reporting_group_members = _merge_catalog_groups(
        catalogs_by_sample.values()
    )
    evidence = classify_cross_sample_carryover(
        signals_by_sample,
        reporting_groups=reporting_groups,
        signal_groups=signal_groups,
        policy=cross_sample_policy,
    )
    results: dict[str, IdentificationResult] = {}
    for sample_id, signals in signals_by_sample.items():
        sample_dir = sample_directories[sample_id]
        catalog = catalogs_by_sample[sample_id]
        sample_evidence = evidence.get(sample_id, ())
        excluded_keys = {
            _signal_key(
                TaxonomicSignal(
                    species_taxid=item.species_taxid,
                    species_name=item.species_name,
                    read_count=0,
                    read_percent=0.0,
                ),
                signal_groups,
            )
            for item in sample_evidence
            if item.status == "possible_cross_talk"
        }
        reviewed_signals = tuple(
            signal
            for signal in signals
            if _signal_key(signal, signal_groups)
            not in excluded_keys
        )
        mixed = assess_mixed_signals(
            reviewed_signals,
            policy=mixed_signal_policy,
            reporting_groups=signal_groups,
        )
        prior_payload = json.loads(
            (sample_dir / "species_identification.json").read_text()
        )
        from .assembly_mixture import apply_assembly_evidence
        mixed = apply_assembly_evidence(mixed,
            prior_payload.get("mixed_signal", {}).get("assembly_evidence", {}))
        kmer_support = tuple(
            KmerSpeciesEvidence(**row)
            for row in prior_payload.get("kmerfinder_support", ())
        )
        result = decide_identification(
            sample_id=sample_id,
            ani_hits=parse_skani_output(sample_dir / "skani.tsv", catalog),
            expected_organism=expected_organisms.get(sample_id, "Unknown"),
            kmerfinder_support=kmer_support,
            mixed_signal=mixed,
            policy=identification_policy,
            reporting_group_members=reporting_group_members,
            cohort_review_complete=True,
            platform=platform,
            contamination_evidence=checkm_contamination_evidence(
                checkm_contamination_by_sample.get(sample_id),
                maximum_contamination=maximum_checkm_contamination,
            ),
        )
        possible = [
            item for item in sample_evidence if item.status == "possible_cross_talk"
        ]
        if possible and result.status == "confirmed_species":
            names = ", ".join(
                item.reporting_group or item.species_name for item in possible
            )
            result = replace(
                result,
                reason=(
                    f"{result.reason} Weak secondary signal(s) ({names}) were "
                    "retained as possible cross-sample carryover."
                ),
            )
        result = replace(
            result,
            mixed_signal_thresholds=asdict(mixed_signal_policy),
            cross_sample_evidence=sample_evidence,
            cross_sample_thresholds=asdict(cross_sample_policy),
        )
        results[sample_id] = result

    _write_batch_audit(root, evidence, cross_sample_policy)
    for sample_id, result in results.items():
        sample_directory = sample_directories[sample_id]
        write_identification_outputs(result, sample_directory)
        placement = _publish_reconciled_placement(
            sample_id=sample_id,
            result=result,
            sample_directory=sample_directory,
            assembly=(assemblies or {}).get(sample_id),
            gtdb_metadata_tsv=gtdb_metadata_tsv,
            gtdbtk_command_prefix=gtdbtk_command_prefix,
            placement_threads=placement_threads,
        )
        _update_sample_provenance(
            sample_directory,
            result,
            cross_sample_policy,
            placement=placement,
        )
    return results


def _merge_catalog_groups(
    catalogs: Sequence[ReferenceCatalog],
) -> tuple[dict[str, str], dict[str, str], dict[str, tuple[str, ...]]]:
    reporting: dict[str, str] = {}
    signals: dict[str, str] = {}
    members: dict[str, set[str]] = {}
    for catalog in catalogs:
        for target, source, label in (
            (reporting, catalog.reporting_groups_by_taxid, "reporting"),
            (signals, catalog.mixed_signal_groups_by_taxid, "mixed-signal"),
        ):
            for taxid, group in source.items():
                previous = target.setdefault(taxid, group)
                if previous != group:
                    raise ValueError(
                        f"Conflicting {label} groups for taxid {taxid}: "
                        f"{previous!r} and {group!r}"
                    )
        for group, names in catalog.reporting_group_member_names.items():
            members.setdefault(group, set()).update(names)
    return (
        reporting,
        signals,
        {
            group: tuple(sorted(names, key=str.casefold))
            for group, names in members.items()
        },
    )


def _publish_reconciled_placement(
    *,
    sample_id: str,
    result: IdentificationResult,
    sample_directory: Path,
    assembly: str | Path | None,
    gtdb_metadata_tsv: str | Path | None,
    gtdbtk_command_prefix: Sequence[str],
    placement_threads: int,
) -> dict[str, object] | None:
    """Publish a tree only after final cohort-level mixture interpretation."""

    if not gtdb_metadata_tsv:
        return None
    output = sample_directory / "phylogenetic_placement"
    if output.exists():
        if output.is_symlink() or output.is_file():
            output.unlink()
        else:
            shutil.rmtree(output)
    if result.status == "confirmed_species":
        return {"status": "not_required", "reason": "Species was confirmed by ANI."}
    if result.status in {
        "mixed",
        "mixed_signal_review",
        "suspected_culture_contamination",
    }:
        payload: dict[str, object] = {
            "status": "suppressed",
            "reason": (
                "A combined marker tree is not reliable when mixed species "
                "signals are present."
            ),
        }
        output.mkdir(parents=True)
        _write_json_atomic(output / "phylogenetic_placement.json", payload)
        return payload
    if not assembly:
        payload = {
            "status": "unavailable",
            "reason": "The assembled genome was not available for placement.",
        }
        output.mkdir(parents=True)
        _write_json_atomic(output / "phylogenetic_placement.json", payload)
        return payload
    try:
        placement = run_bac120_placement(
            assembly=assembly,
            sample_id=sample_id,
            output_directory=output,
            gtdb_metadata_tsv=gtdb_metadata_tsv,
            cpus=max(1, int(placement_threads)),
            pplacer_cpus=max(1, min(int(placement_threads), 4)),
            command_prefix=gtdbtk_command_prefix,
        )
        return asdict(placement)
    except (OSError, RuntimeError, ValueError) as exc:
        payload = {"status": "unavailable", "reason": str(exc)}
        output.mkdir(parents=True, exist_ok=True)
        _write_json_atomic(output / "phylogenetic_placement.json", payload)
        return payload


def _write_json_atomic(path: Path, payload: Mapping[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def expected_organisms_from_summary(path: str | Path) -> dict[str, str]:
    _prefix, fieldnames, rows = _read_summary(path)
    sample_field = "Sample" if "Sample" in fieldnames else "Sample ID"
    expected_field = "Expected organism"
    if sample_field not in fieldnames or expected_field not in fieldnames:
        raise ValueError("QC summary lacks sample or expected-organism columns")
    return {row[sample_field]: row[expected_field] for row in rows}


def checkm_contamination_from_summary(path: str | Path) -> dict[str, str]:
    """Return available per-sample CheckM contamination values from a QC summary."""

    _prefix, fieldnames, rows = _read_summary(path)
    sample_field = "Sample" if "Sample" in fieldnames else "Sample ID"
    contamination_field = "CheckM contamination"
    if sample_field not in fieldnames or contamination_field not in fieldnames:
        return {}
    return {
        row[sample_field]: row[contamination_field]
        for row in rows
        if (row.get(contamination_field) or "").strip().casefold()
        not in {"", "n/a", "na", "none"}
    }


def update_summary_identifications(
    path: str | Path,
    results: Mapping[str, IdentificationResult],
) -> None:
    summary = Path(path)
    prefix, fieldnames, rows = _read_summary(summary)
    sample_field = "Sample" if "Sample" in fieldnames else "Sample ID"
    for row in rows:
        if is_lambda_control(row.get("Expected organism")):
            continue
        result = results.get(row.get(sample_field, ""))
        if result is None:
            continue
        fields = legacy_summary_fields(result)
        from bactipipe.illumina_filtering import refresh_depth_row
        refresh_depth_row(row, fields, summary)
        if "Identified organism" in fieldnames:
            row["Identified organism"] = fields.identified_organism
        if "% Match" in fieldnames:
            row["% Match"] = fields.match_percent
        if "tax_confirm" in fieldnames:
            row["tax_confirm"] = (
                "Contamination" if result.status == "suspected_culture_contamination"
                else fields.taxonomic_confirmation
            )
        if "Overall Quality" in fieldnames:
            if "Isolate QC" not in fieldnames:
                fieldnames.append("Isolate QC")
            row["Isolate QC"] = ("Fail" if result.status == "suspected_culture_contamination"
                                 else "Review" if result.status == "mixed_signal_review" else "Pass")
    temporary = summary.with_suffix(summary.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        if prefix:
            handle.write("\n".join(prefix) + "\n")
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, summary)


def update_final_summary_verdicts(
    final_summary_path: str | Path,
    temporary_summary_path: str | Path,
    *,
    minimum_completeness: float = 90.0,
    maximum_contamination: float = 10.0,
) -> None:
    """Recalculate merged final verdicts after a cohort-context update."""

    final_path = Path(final_summary_path)
    prefix, fieldnames, rows = _read_summary(final_path)
    _temp_prefix, temp_fields, temp_rows = _read_summary(temporary_summary_path)
    temp_sample_field = "Sample" if "Sample" in temp_fields else "Sample ID"
    temp_by_sample = {row[temp_sample_field]: row for row in temp_rows}
    sample_field = "Sample" if "Sample" in fieldnames else "Sample ID"
    required = {
        "CheckM completeness",
        "CheckM contamination",
        "Overall Quality",
    }
    if not required.issubset(fieldnames):
        return
    for row in rows:
        temporary_row = temp_by_sample.get(row.get(sample_field, ""))
        if temporary_row is None:
            continue
        from bactipipe.qc_status import threshold_verdict
        completeness_verdict = threshold_verdict(row["CheckM completeness"], minimum_completeness)
        contamination_verdict = threshold_verdict(row["CheckM contamination"], maximum_contamination, maximum=True)
        from bactipipe.qc_status import interpret
        row["Overall Quality"] = interpret(
            row.get("Expected organism", ""), row.get("Identified organism", ""),
            [temporary_row.get("qc_verdict", "N/A"), temporary_row.get("cov_verdict", "N/A"),
             completeness_verdict, contamination_verdict],
            isolate_qc_failure=row.get("Isolate QC") == "Fail",
            isolate_qc_review=row.get("Isolate QC") == "Review",
            control_qc=row.get("Control QC", "Not assessed"),
        )["overall"]
    temporary = final_path.with_suffix(final_path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        if prefix:
            handle.write("\n".join(prefix) + "\n")
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, final_path)


def _read_summary(
    path: str | Path,
) -> tuple[list[str], list[str], list[dict[str, str]]]:
    summary = Path(path)
    lines = summary.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = next(
        (
            index
            for index, line in enumerate(lines)
            if line.startswith("Sample\t") or line.startswith("Sample ID\t")
        ),
        None,
    )
    if header_index is None:
        raise ValueError(f"QC summary header was not found: {summary}")
    reader = csv.DictReader(lines[header_index:], delimiter="\t")
    return lines[:header_index], list(reader.fieldnames or ()), list(reader)


def _write_batch_audit(
    root: Path,
    evidence: Mapping[str, Sequence[CrossSampleEvidence]],
    policy: CrossSamplePolicy,
) -> None:
    # A control-only batch (including batches whose bacteria failed QC) has no
    # species directory yet. An empty audit is valid, not a run-wide failure.
    root.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample_id",
        "species_taxid",
        "species_name",
        "reporting_group",
        "status",
        "source_sample",
        "mapped_reads",
        "source_mapped_reads",
        "fraction_of_source",
        "fraction_of_sample_primary",
    ]
    rows = [
        {
            "sample_id": sample_id,
            **asdict(item),
            "fraction_of_source": f"{item.fraction_of_source:.6f}",
            "fraction_of_sample_primary": f"{item.fraction_of_sample_primary:.6f}",
        }
        for sample_id, items in sorted(evidence.items())
        for item in items
    ]
    temporary = root / "cross_sample_signal_review.tsv.tmp"
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, root / "cross_sample_signal_review.tsv")
    payload = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "policy": asdict(policy),
        "sample_count": len(evidence),
        "possible_cross_talk_count": sum(
            item.status == "possible_cross_talk"
            for items in evidence.values()
            for item in items
        ),
    }
    temporary_json = root / "cross_sample_signal_review.json.tmp"
    temporary_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary_json, root / "cross_sample_signal_review.json")


def _update_sample_provenance(
    sample_directory: Path,
    result: IdentificationResult,
    policy: CrossSamplePolicy,
    *,
    placement: Mapping[str, object] | None = None,
) -> None:
    path = sample_directory / "analysis_provenance.json"
    try:
        payload = json.loads(path.read_text()) if path.is_file() else {}
    except (OSError, ValueError, TypeError):
        payload = {}
    payload.update(
        {
            "schema_version": max(3, int(payload.get("schema_version") or 1)),
            "result_status": result.status,
            "contamination_evidence": asdict(result.contamination_evidence),
            "cohort_review": {
                "performed": True,
                "completed_at": datetime.now(timezone.utc).isoformat(),
                "policy": asdict(policy),
                "possible_cross_talk": [
                    asdict(item)
                    for item in result.cross_sample_evidence
                    if item.status == "possible_cross_talk"
                ],
            },
        }
    )
    if placement is not None:
        payload["phylogenetic_placement"] = dict(placement)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)
