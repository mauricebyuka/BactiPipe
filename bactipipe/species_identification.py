"""Evidence-based bacterial species identification.

This module deliberately separates evidence parsing from the diagnostic
decision. Whole-genome ANI is the primary species evidence. Read-classifier
signals are candidate nominations only and are independently confirmed before
they can stop a mixed isolate from being forced into a single-species result.
KmerFinder parsing remains solely for optional legacy-result comparison.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from .organism_names import expected_organism_name


SCHEMA_VERSION = 3
UNKNOWN_NAMES = {"", "unknown", "n/a", "na", "none"}


def _normalise_name(value: str | None) -> str:
    return " ".join((value or "").strip().split()).casefold()


def _finite_float(value: object, *, field_name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be numeric; received {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{field_name} must be finite; received {value!r}")
    return result


def _optional_float(value: object) -> float | None:
    if value is None or str(value).strip().casefold() in UNKNOWN_NAMES:
        return None
    return _finite_float(value, field_name="value")


def _optional_int(value: object) -> int | None:
    if value is None or str(value).strip().casefold() in UNKNOWN_NAMES:
        return None
    try:
        return int(str(value).strip())
    except (TypeError, ValueError) as exc:
        raise ValueError(f"value must be an integer; received {value!r}") from exc


@dataclass(frozen=True)
class IdentificationPolicy:
    """Candidate thresholds pending laboratory validation."""

    minimum_ani: float = 95.0
    minimum_aligned_fraction: float = 50.0
    minimum_candidate_aligned_fraction: float = 30.0
    minimum_ani_margin: float = 0.5
    minimum_aligned_fraction_margin: float = 2.0
    maximum_closest_species: int = 10
    aligned_fraction_mode: str = "either"

    def aligned_fraction(self, hit: AniHit) -> float:
        """GTDB-Tk uses max(query AF, reference AF); retain legacy replay."""
        fractions = (hit.aligned_fraction_query, hit.aligned_fraction_reference)
        return max(fractions) if self.aligned_fraction_mode == "either" else min(fractions)

    @classmethod
    def from_saved_thresholds(cls, values):
        if not values:
            return cls()
        # Existing evidence predates the explicit mode and used both AFs >=65%.
        saved = dict(values)
        if 'aligned_fraction_mode' not in saved:
            saved['aligned_fraction_mode'] = 'both'
            saved.setdefault('minimum_aligned_fraction', 65.0)
        return cls(**saved)

    def __post_init__(self) -> None:
        if self.aligned_fraction_mode not in {"either", "both"}:
            raise ValueError("aligned_fraction_mode must be either or both")
        for name in (
            "minimum_ani",
            "minimum_aligned_fraction",
            "minimum_candidate_aligned_fraction",
            "minimum_ani_margin",
            "minimum_aligned_fraction_margin",
        ):
            value = getattr(self, name)
            if not 0 <= value <= 100:
                raise ValueError(f"{name} must be between 0 and 100")
        if self.minimum_candidate_aligned_fraction > self.minimum_aligned_fraction:
            raise ValueError(
                "minimum_candidate_aligned_fraction cannot exceed "
                "minimum_aligned_fraction"
            )
        if self.maximum_closest_species < 1:
            raise ValueError("maximum_closest_species must be at least 1")


@dataclass(frozen=True)
class MixedSignalPolicy:
    """Conservative gates for a reportable secondary bacterial signal.

    These are deliberately configurable and are not a claim of validated
    diagnostic cutoffs.  A secondary signal requires read count, abundance,
    and breadth evidence; a few shared or horizontally transferred loci are
    insufficient.
    """

    minimum_secondary_reads: int = 100
    minimum_secondary_distinct_minimizers: int = 100
    minimum_secondary_read_percent: float = 5.0
    minimum_secondary_reference_breadth: float = 20.0
    comparable_to_primary_ratio: float = 0.25
    minimum_contamination_read_percent: float = 10.0
    minimum_contamination_reference_breadth: float = 50.0
    minimum_contamination_primary_breadth: float = 70.0

    def __post_init__(self) -> None:
        if self.minimum_secondary_reads < 1:
            raise ValueError("minimum_secondary_reads must be at least 1")
        if self.minimum_secondary_distinct_minimizers < 1:
            raise ValueError(
                "minimum_secondary_distinct_minimizers must be at least 1"
            )
        for name in (
            "minimum_secondary_read_percent",
            "minimum_secondary_reference_breadth",
            "minimum_contamination_read_percent",
            "minimum_contamination_reference_breadth",
            "minimum_contamination_primary_breadth",
        ):
            if not 0 <= getattr(self, name) <= 100:
                raise ValueError(f"{name} must be between 0 and 100")
        if self.minimum_contamination_read_percent < self.minimum_secondary_read_percent:
            raise ValueError(
                "minimum_contamination_read_percent cannot be lower than "
                "minimum_secondary_read_percent"
            )
        if (
            self.minimum_contamination_reference_breadth
            < self.minimum_secondary_reference_breadth
        ):
            raise ValueError(
                "minimum_contamination_reference_breadth cannot be lower than "
                "minimum_secondary_reference_breadth"
            )
        if not 0 <= self.comparable_to_primary_ratio <= 1:
            raise ValueError("comparable_to_primary_ratio must be between 0 and 1")


def default_mixed_signal_policy(platform: str) -> MixedSignalPolicy:
    """Return provisional platform-specific gates for validation.

    Illumina's short reads receive stricter read, distinct-minimizer, and
    competitive-breadth gates because closely related bacteria can share many
    short exact matches. These are candidate validation settings, not approved
    diagnostic cutoffs.
    """

    normalised = platform.strip().casefold()
    if normalised == "illumina":
        return MixedSignalPolicy(
            minimum_secondary_reads=250,
            minimum_secondary_distinct_minimizers=250,
            minimum_secondary_read_percent=5.0,
            minimum_secondary_reference_breadth=25.0,
            comparable_to_primary_ratio=0.25,
            minimum_contamination_read_percent=10.0,
            minimum_contamination_reference_breadth=50.0,
            minimum_contamination_primary_breadth=70.0,
        )
    if normalised == "nanopore":
        return MixedSignalPolicy()
    raise ValueError("platform must be Illumina or Nanopore")


@dataclass(frozen=True)
class ContaminationEvidence:
    """Assembly-level evidence used to corroborate read-level mixture signals."""

    status: str = "not_assessed"
    checkm_contamination: float | None = None
    maximum_checkm_contamination: float = 10.0
    reason: str = "Assembly-level contamination evidence was not supplied."

    def __post_init__(self) -> None:
        if self.status not in {
            "not_assessed",
            "supports_contamination",
            "does_not_support_contamination",
        }:
            raise ValueError(f"Unsupported contamination-evidence status: {self.status}")
        if self.checkm_contamination is not None:
            _finite_float(
                self.checkm_contamination,
                field_name="CheckM contamination",
            )
        if not 0 <= self.maximum_checkm_contamination <= 100:
            raise ValueError(
                "maximum_checkm_contamination must be between 0 and 100"
            )


def checkm_contamination_evidence(
    value: object,
    *,
    maximum_contamination: float = 10.0,
) -> ContaminationEvidence:
    """Interpret CheckM as corroborating QC evidence, never species evidence."""

    threshold = _finite_float(
        maximum_contamination,
        field_name="maximum CheckM contamination",
    )
    if not 0 <= threshold <= 100:
        raise ValueError("maximum CheckM contamination must be between 0 and 100")
    contamination = _optional_float(value)
    if contamination is None:
        return ContaminationEvidence(maximum_checkm_contamination=threshold)
    if contamination > threshold:
        return ContaminationEvidence(
            status="supports_contamination",
            checkm_contamination=contamination,
            maximum_checkm_contamination=threshold,
            reason=(
                f"CheckM contamination ({contamination:.2f}%) exceeded the "
                f"{threshold:.2f}% isolate-QC limit."
            ),
        )
    return ContaminationEvidence(
        status="does_not_support_contamination",
        checkm_contamination=contamination,
        maximum_checkm_contamination=threshold,
        reason=(
            f"CheckM contamination ({contamination:.2f}%) did not exceed the "
            f"{threshold:.2f}% isolate-QC limit."
        ),
    )


def default_kraken_confidence(platform: str) -> float:
    """Return the provisional classifier threshold used only for nomination."""

    normalised = platform.strip().casefold()
    if normalised == "illumina":
        return 0.10
    if normalised == "nanopore":
        return 0.05
    raise ValueError("platform must be Illumina or Nanopore")


@dataclass(frozen=True)
class ReferenceRecord:
    accession: str
    species_taxid: str
    species_name: str
    fasta_path: str
    source: str
    reference_role: str
    assembly_level: str
    sequence_origin: str
    type_material: bool
    canonical_for_tree: bool
    callable_for_identification: bool
    include_in_kraken: bool = False
    reporting_group: str = ""
    ncbi_taxonomy_name: str = ""
    gtdb_taxonomy: str = ""
    sha256: str = ""


class ReferenceCatalog:
    REQUIRED_COLUMNS = {
        "accession",
        "species_taxid",
        "species_name",
        "fasta_path",
        "source",
        "reference_role",
        "assembly_level",
        "sequence_origin",
        "type_material",
        "canonical_for_tree",
        "callable_for_identification",
    }

    def __init__(self, records: Sequence[ReferenceRecord], *, root: Path) -> None:
        self.records = tuple(records)
        self.root = root.resolve()
        self._by_accession = {record.accession: record for record in records}
        self._by_path: dict[str, ReferenceRecord] = {}
        for record in records:
            path = Path(record.fasta_path)
            resolved = path if path.is_absolute() else self.root / path
            self._by_path[str(resolved.resolve())] = record
            self._by_path[path.name] = record
        self.reporting_groups_by_taxid = {
            record.species_taxid: record.reporting_group
            for record in records
            if record.reporting_group
        }
        self.mixed_signal_groups_by_taxid = dict(self.reporting_groups_by_taxid)
        taxonomy_name_taxids: dict[str, set[str]] = {}
        for record in records:
            if not record.callable_for_identification:
                continue
            taxonomy_name = _normalise_name(record.ncbi_taxonomy_name)
            if taxonomy_name:
                taxonomy_name_taxids.setdefault(taxonomy_name, set()).add(
                    record.species_taxid
                )
        for taxonomy_name, taxids in taxonomy_name_taxids.items():
            if len(taxids) < 2:
                continue
            group = f"ncbi_species:{taxonomy_name}"
            for taxid in taxids:
                self.mixed_signal_groups_by_taxid.setdefault(taxid, group)
        self.identification_names_by_taxid: dict[str, str] = {}
        for record in records:
            if not record.callable_for_identification:
                continue
            previous = self.identification_names_by_taxid.setdefault(
                record.species_taxid, record.species_name
            )
            if previous != record.species_name:
                raise ValueError(
                    "Reference manifest gives conflicting identification names "
                    f"for taxid {record.species_taxid}: {previous!r} and "
                    f"{record.species_name!r}"
                )
        grouped_names: dict[str, set[str]] = {}
        for record in records:
            if record.reporting_group:
                grouped_names.setdefault(record.reporting_group, set()).add(
                    record.species_name
                )
        self.reporting_group_member_names = {
            group: tuple(sorted(names, key=str.casefold))
            for group, names in grouped_names.items()
        }

    @classmethod
    def from_tsv(cls, path: str | Path) -> "ReferenceCatalog":
        manifest_path = Path(path)
        with manifest_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing = cls.REQUIRED_COLUMNS.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    "Reference manifest is missing columns: " + ", ".join(sorted(missing))
                )
            records: list[ReferenceRecord] = []
            seen_accessions: set[str] = set()
            for line_number, row in enumerate(reader, start=2):
                accession = (row.get("accession") or "").strip()
                species_taxid = (row.get("species_taxid") or "").strip()
                species_name = " ".join((row.get("species_name") or "").split())
                fasta_path = (row.get("fasta_path") or "").strip()
                if not accession or accession in seen_accessions:
                    raise ValueError(
                        f"Invalid or duplicate accession on manifest line {line_number}: "
                        f"{accession!r}"
                    )
                if not species_taxid.isdigit() or not species_name or not fasta_path:
                    raise ValueError(
                        f"Incomplete species metadata on manifest line {line_number}"
                    )
                records.append(
                    ReferenceRecord(
                        accession=accession,
                        species_taxid=species_taxid,
                        species_name=species_name,
                        fasta_path=fasta_path,
                        source=(row.get("source") or "").strip(),
                        reference_role=(row.get("reference_role") or "").strip(),
                        assembly_level=(row.get("assembly_level") or "").strip(),
                        sequence_origin=(row.get("sequence_origin") or "").strip(),
                        type_material=_parse_bool(
                            row.get("type_material"), line_number=line_number
                        ),
                        canonical_for_tree=_parse_bool(
                            row.get("canonical_for_tree"), line_number=line_number
                        ),
                        callable_for_identification=_parse_bool(
                            row.get("callable_for_identification"), line_number=line_number
                        ),
                        include_in_kraken=_parse_bool(
                            row.get("include_in_kraken")
                            if "include_in_kraken" in row
                            else row.get("callable_for_identification"),
                            line_number=line_number,
                        ),
                        reporting_group=(row.get("reporting_group") or "").strip(),
                        ncbi_taxonomy_name=(row.get("ncbi_taxonomy_name") or "").strip(),
                        gtdb_taxonomy=(row.get("gtdb_taxonomy") or "").strip(),
                        sha256=(row.get("sha256") or "").strip().lower(),
                    )
                )
                seen_accessions.add(accession)
        if not records:
            raise ValueError("Reference manifest contains no records")
        return cls(records, root=manifest_path.parent)

    def resolve(self, reference_path: str) -> ReferenceRecord:
        candidate = Path(reference_path)
        for key in (str(candidate.resolve()), candidate.name):
            if key in self._by_path:
                return self._by_path[key]
        accession = next(
            (
                value
                for value in sorted(self._by_accession, key=len, reverse=True)
                if candidate.name == value
                or candidate.name.startswith(value + ".")
                or candidate.name.startswith(value + "_")
            ),
            None,
        )
        if accession is not None:
            return self._by_accession[accession]
        raise ValueError(
            f"skani reference {reference_path!r} is absent from the reviewed manifest"
        )

    def canonical_for_species(self, species_taxid: str) -> ReferenceRecord | None:
        matches = [
            record
            for record in self.records
            if record.species_taxid == str(species_taxid)
            and record.canonical_for_tree
            and record.callable_for_identification
        ]
        if len(matches) > 1:
            raise ValueError(
                f"Species taxid {species_taxid} has more than one canonical reference"
            )
        return matches[0] if matches else None

    def fasta_for(self, record: ReferenceRecord) -> Path:
        path = Path(record.fasta_path)
        return path if path.is_absolute() else self.root / path

    def validate_files(self, *, require_checksums: bool = True) -> None:
        canonical_counts: dict[str, int] = {}
        callable_species: set[str] = set()
        kraken_species: set[str] = set()
        reporting_group_by_species: dict[str, str] = {}
        allowed_roles = {"type", "reference", "representative", "diversity", "context"}
        for record in self.records:
            previous_group = reporting_group_by_species.setdefault(
                record.species_taxid, record.reporting_group
            )
            if previous_group != record.reporting_group:
                raise ValueError(
                    f"Species taxid {record.species_taxid} has inconsistent reporting groups"
                )
            path = Path(record.fasta_path)
            path = path if path.is_absolute() else self.root / path
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Reference FASTA is missing or empty: {path}")
            if require_checksums and not record.sha256:
                raise ValueError(f"Reference {record.accession} has no SHA-256 checksum")
            if record.sha256:
                observed = _sha256(path)
                if observed != record.sha256:
                    raise ValueError(
                        f"Reference checksum mismatch for {record.accession}: "
                        f"expected {record.sha256}, observed {observed}"
                    )
            if record.reference_role.casefold() not in allowed_roles:
                raise ValueError(
                    f"Reference {record.accession} has unsupported role "
                    f"{record.reference_role!r}"
                )
            if record.callable_for_identification:
                callable_species.add(record.species_taxid)
                if record.sequence_origin.casefold() != "isolate":
                    raise ValueError(
                        f"Callable reference {record.accession} is not isolate-derived"
                    )
            if record.include_in_kraken:
                kraken_species.add(record.species_taxid)
                if not record.callable_for_identification:
                    raise ValueError(
                        f"Kraken reference {record.accession} is not callable for "
                        "identification"
                    )
            if record.canonical_for_tree and not record.callable_for_identification:
                raise ValueError(
                    f"Canonical reference {record.accession} is not callable for "
                    "identification"
                )
            if record.canonical_for_tree:
                canonical_counts[record.species_taxid] = (
                    canonical_counts.get(record.species_taxid, 0) + 1
                )
        invalid = [
            taxid for taxid in callable_species if canonical_counts.get(taxid, 0) != 1
        ]
        if invalid:
            raise ValueError(
                "Each represented tree species must have one canonical reference: "
                + ", ".join(sorted(invalid))
            )
        missing_kraken = callable_species.difference(kraken_species)
        if missing_kraken:
            raise ValueError(
                "Each callable species must have at least one Kraken reference: "
                + ", ".join(sorted(missing_kraken))
            )


def _parse_bool(value: object, *, line_number: int) -> bool:
    normalised = str(value or "").strip().casefold()
    if normalised in {"1", "true", "yes", "y"}:
        return True
    if normalised in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean {value!r} on manifest line {line_number}")


@dataclass(frozen=True)
class AniHit:
    accession: str
    species_taxid: str
    species_name: str
    reference_path: str
    ani: float
    aligned_fraction_reference: float
    aligned_fraction_query: float
    reporting_group: str = ""

    @property
    def minimum_aligned_fraction(self) -> float:
        return min(self.aligned_fraction_reference, self.aligned_fraction_query)


@dataclass(frozen=True)
class KmerSpeciesEvidence:
    species_taxid: str
    species_name: str
    query_coverage: float
    template_coverage: float
    score: float
    hit_count: int


@dataclass(frozen=True)
class TaxonomicSignal:
    species_taxid: str
    species_name: str
    read_count: int
    read_percent: float
    reference_breadth: float | None = None
    distinct_minimizers: int | None = None
    competitive_mapped_reads: int | None = None
    competitive_mapped_percent: float | None = None
    nomination_score: float | None = None
    nomination_percent: float | None = None

    @property
    def confirmed_read_count(self) -> int:
        if self.competitive_mapped_reads is not None:
            return self.competitive_mapped_reads
        return self.read_count

    @property
    def confirmed_read_percent(self) -> float:
        if self.competitive_mapped_percent is not None:
            return self.competitive_mapped_percent
        return self.read_percent


@dataclass(frozen=True)
class MixedSignalAssessment:
    status: str
    primary: TaxonomicSignal | None
    qualifying_secondary: tuple[TaxonomicSignal, ...] = ()
    reason: str = ""


@dataclass(frozen=True)
class CrossSampleEvidence:
    species_taxid: str
    species_name: str
    reporting_group: str
    status: str
    source_sample: str
    mapped_reads: int
    source_mapped_reads: int
    fraction_of_source: float
    fraction_of_sample_primary: float


@dataclass(frozen=True)
class _SignalGroupEvidence:
    key: str
    primary: TaxonomicSignal
    confirmed_read_count: int
    confirmed_read_percent: float
    distinct_minimizers: int | None
    reference_breadth: float | None


def _signal_group_evidence(
    signals: Sequence[TaxonomicSignal],
    reporting_groups: Mapping[str, str],
) -> list[_SignalGroupEvidence]:
    grouped: dict[str, list[TaxonomicSignal]] = {}
    for signal in signals:
        reporting_group = reporting_groups.get(signal.species_taxid, "").strip()
        key = (
            f"reporting_group:{reporting_group.casefold()}"
            if reporting_group
            else f"species_taxid:{signal.species_taxid}"
        )
        grouped.setdefault(key, []).append(signal)
    evidence = []
    for key, members in grouped.items():
        primary = max(
            members,
            key=lambda item: (item.confirmed_read_count, item.confirmed_read_percent),
        )
        minimizers = [
            item.distinct_minimizers
            for item in members
            if item.distinct_minimizers is not None
        ]
        breadths = [
            item.reference_breadth
            for item in members
            if item.reference_breadth is not None
        ]
        evidence.append(
            _SignalGroupEvidence(
                key=key,
                primary=primary,
                confirmed_read_count=sum(
                    item.confirmed_read_count for item in members
                ),
                confirmed_read_percent=sum(
                    item.confirmed_read_percent for item in members
                ),
                distinct_minimizers=max(minimizers) if minimizers else None,
                reference_breadth=max(breadths) if breadths else None,
            )
        )
    return sorted(
        evidence,
        key=lambda item: (item.confirmed_read_count, item.confirmed_read_percent),
        reverse=True,
    )


@dataclass(frozen=True)
class IdentificationResult:
    sample_id: str
    status: str
    reported_organism: str
    match_percent: float | None
    taxonomic_confirmation: str
    reason: str
    closest_species: tuple[AniHit, ...] = ()
    kmerfinder_support: tuple[KmerSpeciesEvidence, ...] = ()
    mixed_signal: MixedSignalAssessment = field(
        default_factory=lambda: MixedSignalAssessment(
            status="not_assessed", primary=None, reason="Secondary-signal assessment was not performed."
        )
    )
    identification_thresholds: dict[str, object] = field(default_factory=dict)
    mixed_signal_thresholds: dict[str, object] = field(default_factory=dict)
    cross_sample_evidence: tuple[CrossSampleEvidence, ...] = ()
    cross_sample_thresholds: dict[str, object] = field(default_factory=dict)
    contamination_evidence: ContaminationEvidence = field(
        default_factory=ContaminationEvidence
    )
    schema_version: int = SCHEMA_VERSION

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def parse_skani_output(
    path: str | Path, catalog: ReferenceCatalog
) -> tuple[AniHit, ...]:
    required = {
        "Ref_file",
        "ANI",
        "Align_fraction_ref",
        "Align_fraction_query",
    }
    hits: list[AniHit] = []
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError("skani output is missing columns: " + ", ".join(sorted(missing)))
        for line_number, row in enumerate(reader, start=2):
            record = catalog.resolve(row["Ref_file"])
            if not record.callable_for_identification:
                continue
            ani = _finite_float(row.get("ANI"), field_name=f"ANI line {line_number}")
            af_ref = _finite_float(
                row.get("Align_fraction_ref"),
                field_name=f"reference aligned fraction line {line_number}",
            )
            af_query = _finite_float(
                row.get("Align_fraction_query"),
                field_name=f"query aligned fraction line {line_number}",
            )
            if not all(0 <= value <= 100 for value in (ani, af_ref, af_query)):
                raise ValueError(f"Out-of-range skani value on line {line_number}")
            hits.append(
                AniHit(
                    accession=record.accession,
                    species_taxid=record.species_taxid,
                    species_name=record.species_name,
                    reference_path=record.fasta_path,
                    ani=ani,
                    aligned_fraction_reference=af_ref,
                    aligned_fraction_query=af_query,
                    reporting_group=record.reporting_group,
                )
            )
    return tuple(hits)


def parse_kmerfinder_output(path: str | Path) -> tuple[KmerSpeciesEvidence, ...]:
    """Parse KmerFinder without treating query coverage as percent identity."""

    by_species: dict[tuple[str, str], dict[str, float | int]] = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Species", "TAXID Species", "Query_Coverage", "Template_Coverage", "Score"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "KmerFinder output is missing columns: " + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            species = " ".join((row.get("Species") or "").split())
            taxid = (row.get("TAXID Species") or "").strip()
            if not species:
                continue
            query_value = row.get("tot_query_Coverage") or row.get("Query_Coverage")
            template_value = row.get("tot_template_Coverage") or row.get("Template_Coverage")
            query_coverage = _finite_float(
                query_value, field_name=f"KmerFinder query coverage line {line_number}"
            )
            template_coverage = _finite_float(
                template_value, field_name=f"KmerFinder template coverage line {line_number}"
            )
            score = _finite_float(
                row.get("Score"), field_name=f"KmerFinder score line {line_number}"
            )
            entry = by_species.setdefault(
                (taxid, species),
                {"query": 0.0, "template": 0.0, "score": 0.0, "count": 0},
            )
            entry["query"] = max(float(entry["query"]), min(100.0, query_coverage))
            entry["template"] = max(
                float(entry["template"]), min(100.0, template_coverage)
            )
            entry["score"] = float(entry["score"]) + score
            entry["count"] = int(entry["count"]) + 1
    evidence = [
        KmerSpeciesEvidence(
            species_taxid=taxid,
            species_name=species,
            query_coverage=float(values["query"]),
            template_coverage=float(values["template"]),
            score=float(values["score"]),
            hit_count=int(values["count"]),
        )
        for (taxid, species), values in by_species.items()
    ]
    evidence.sort(
        key=lambda item: (item.query_coverage, item.template_coverage, item.score),
        reverse=True,
    )
    return tuple(evidence)


def parse_taxonomic_signals(path: str | Path) -> tuple[TaxonomicSignal, ...]:
    """Load classifier-independent, species-level read evidence.

    A Kraken adapter can create this compact table, but the decision layer is
    intentionally not coupled to Kraken's report format.
    """

    required = {"species_taxid", "species_name", "read_count", "read_percent"}
    signals: list[TaxonomicSignal] = []
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Taxonomic signal table is missing columns: " + ", ".join(sorted(missing))
            )
        for row in reader:
            signals.append(
                TaxonomicSignal(
                    species_taxid=(row.get("species_taxid") or "").strip(),
                    species_name=" ".join((row.get("species_name") or "").split()),
                    read_count=_optional_int(row.get("read_count")) or 0,
                    read_percent=_optional_float(row.get("read_percent")) or 0.0,
                    reference_breadth=_optional_float(row.get("reference_breadth")),
                    distinct_minimizers=_optional_int(row.get("distinct_minimizers")),
                    competitive_mapped_reads=_optional_int(
                        row.get("competitive_mapped_reads")
                    ),
                    competitive_mapped_percent=_optional_float(
                        row.get("competitive_mapped_percent")
                    ),
                    nomination_score=_optional_float(row.get("nomination_score")),
                    nomination_percent=_optional_float(
                        row.get("nomination_percent")
                    ),
                )
            )
    return tuple(signals)


def assess_mixed_signals(
    signals: Sequence[TaxonomicSignal],
    policy: MixedSignalPolicy = MixedSignalPolicy(),
    reporting_groups: Mapping[str, str] | None = None,
) -> MixedSignalAssessment:
    if not signals:
        return MixedSignalAssessment(
            status="not_assessed",
            primary=None,
            reason="Secondary-signal assessment was not performed.",
        )
    reporting_groups = reporting_groups or {}
    competitive_mapping_available = any(
        signal.competitive_mapped_reads is not None for signal in signals
    )
    if competitive_mapping_available:
        mapped_signals = [
            signal
            for signal in signals
            if signal.competitive_mapped_reads is not None
        ]
        mapped_group_keys = {
            item.key for item in _signal_group_evidence(mapped_signals, reporting_groups)
        }
        unmapped_groups = _signal_group_evidence(
            [
                signal
                for signal in signals
                if signal.competitive_mapped_reads is None
            ],
            reporting_groups,
        )
        unmapped_relevant_signal = any(
            group.key not in mapped_group_keys
            and group.confirmed_read_count >= policy.minimum_secondary_reads
            and group.confirmed_read_percent >= policy.minimum_secondary_read_percent
            for group in unmapped_groups
        )
    else:
        mapped_signals = list(signals)
        unmapped_relevant_signal = False
    ordered = _signal_group_evidence(mapped_signals, reporting_groups)
    if not ordered:
        return MixedSignalAssessment(
            status="indeterminate",
            primary=None,
            reason="No classifier candidate received competitive mapping evidence.",
        )
    primary_group = ordered[0]
    primary = primary_group.primary
    qualifying: list[_SignalGroupEvidence] = []
    incomplete_evidence = unmapped_relevant_signal
    for signal_group in ordered[1:]:
        if (
            signal_group.confirmed_read_count < policy.minimum_secondary_reads
            or signal_group.confirmed_read_percent
            < policy.minimum_secondary_read_percent
        ):
            continue
        if signal_group.distinct_minimizers is None:
            incomplete_evidence = True
            continue
        if (
            signal_group.distinct_minimizers
            < policy.minimum_secondary_distinct_minimizers
        ):
            continue
        if signal_group.reference_breadth is None:
            incomplete_evidence = True
            continue
        if (
            signal_group.reference_breadth
            >= policy.minimum_secondary_reference_breadth
        ):
            qualifying.append(signal_group)
    if qualifying:
        strongest = qualifying[0]
        ratio = strongest.confirmed_read_count / max(
            primary_group.confirmed_read_count, 1
        )
        strong_contamination_signal = bool(
            ratio >= policy.comparable_to_primary_ratio
            and strongest.confirmed_read_percent
            >= policy.minimum_contamination_read_percent
            and strongest.reference_breadth is not None
            and strongest.reference_breadth
            >= policy.minimum_contamination_reference_breadth
            and primary_group.reference_breadth is not None
            and primary_group.reference_breadth
            >= policy.minimum_contamination_primary_breadth
        )
        status = (
            "strong_secondary"
            if strong_contamination_signal
            else "minor_secondary"
        )
        return MixedSignalAssessment(
            status=status,
            primary=primary,
            qualifying_secondary=tuple(item.primary for item in qualifying),
            reason=(
                "At least one secondary species or reporting group has sufficient "
                "competitive reads, meaningful abundance, and distributed genome "
                "evidence."
            ),
        )
    if incomplete_evidence:
        return MixedSignalAssessment(
            status="indeterminate",
            primary=primary,
            reason=(
                "Secondary taxonomic signals were present, but distinct-minimizer "
                "or reference-breadth evidence was unavailable."
            ),
        )
    return MixedSignalAssessment(
        status="single_signal",
        primary=primary,
        reason="No secondary species passed all mixed-signal evidence gates.",
    )


def decide_identification(
    *,
    sample_id: str,
    ani_hits: Sequence[AniHit],
    expected_organism: str = "Unknown",
    kmerfinder_support: Sequence[KmerSpeciesEvidence] = (),
    mixed_signal: MixedSignalAssessment | None = None,
    policy: IdentificationPolicy = IdentificationPolicy(),
    reporting_group_members: Mapping[str, Sequence[str]] | None = None,
    cohort_review_complete: bool = False,
    platform: str = "",
    contamination_evidence: ContaminationEvidence | None = None,
) -> IdentificationResult:
    mixed_signal = mixed_signal or MixedSignalAssessment(
        status="not_assessed", primary=None, reason="Secondary-signal assessment was not performed."
    )
    thresholds = asdict(policy)
    contamination_evidence = contamination_evidence or ContaminationEvidence()
    normalised_platform = platform.strip().casefold()
    if normalised_platform and normalised_platform not in {"illumina", "nanopore"}:
        raise ValueError("platform must be Illumina or Nanopore")
    species_hits = _best_hit_per_species(ani_hits, policy=policy)
    closest = tuple(
        hit
        for hit in species_hits
        if policy.aligned_fraction(hit) >= policy.minimum_candidate_aligned_fraction
    )[: policy.maximum_closest_species]

    passing = sorted(
        [
            hit
            for hit in species_hits
            if hit.ani >= policy.minimum_ani
            and policy.aligned_fraction(hit) >= policy.minimum_aligned_fraction
        ],
        key=lambda hit: (hit.ani, policy.aligned_fraction(hit)),
        reverse=True,
    )
    top = passing[0] if passing else None
    competitor = (
        next(
            (
                hit
                for hit in passing[1:]
                if top is not None
                and top.ani - hit.ani < policy.minimum_ani_margin
                and (
                    policy.aligned_fraction(top) - policy.aligned_fraction(hit)
                    < policy.minimum_aligned_fraction_margin
                )
            ),
            None,
        )
        if top is not None
        else None
    )

    group_members = reporting_group_members or {}
    accepted_names: set[str] = set()
    report_hit: AniHit | None = None
    if top is not None and competitor is None:
        accepted_names = set(group_members.get(top.reporting_group, ()))
        if top.reporting_group:
            accepted_names.update(
                hit.species_name
                for hit in ani_hits
                if hit.reporting_group == top.reporting_group
            )
        report_hit = _select_user_facing_hit(
            top=top,
            ani_hits=ani_hits,
            expected_organism=expected_organism,
            accepted_group_members=accepted_names,
            policy=policy,
        )
        if report_hit is not top and closest:
            closest = (report_hit, *closest[1:])

    if mixed_signal.status in {"mixed", "strong_secondary"}:
        if not cohort_review_complete:
            return _result(
                sample_id=sample_id,
                status="mixed_signal_review",
                expected_organism=expected_organism,
                reason=(
                    "A strong secondary bacterial signal requires completed "
                    "cohort-level cross-talk review before isolate contamination "
                    "can be reported."
                ),
                closest=closest,
                kmerfinder_support=kmerfinder_support,
                mixed_signal=mixed_signal,
                identification_thresholds=thresholds,
                contamination_evidence=contamination_evidence,
            )
        if (
            normalised_platform == "illumina"
            and contamination_evidence.status != "supports_contamination"
        ):
            return _result(
                sample_id=sample_id,
                status="mixed_signal_review",
                expected_organism=expected_organism,
                reason=(
                    "A strong secondary short-read mapping signal remained after "
                    "cohort cross-talk review, but Illumina competitive mapping "
                    "alone is insufficient to report culture contamination. "
                    f"{contamination_evidence.reason}"
                ),
                closest=closest,
                kmerfinder_support=kmerfinder_support,
                mixed_signal=mixed_signal,
                identification_thresholds=thresholds,
                contamination_evidence=contamination_evidence,
            )
        contamination_reason = (
            "A strong secondary bacterial signal remained after cohort-level "
            "cross-talk review, consistent with possible culture contamination."
        )
        if contamination_evidence.status == "supports_contamination":
            contamination_reason += f" {contamination_evidence.reason}"
        if report_hit is not None:
            reported_organism = report_hit.species_name
            return IdentificationResult(
                sample_id=sample_id,
                status="suspected_culture_contamination",
                reported_organism=reported_organism,
                match_percent=report_hit.ani,
                taxonomic_confirmation=_taxonomic_confirmation(
                    expected_organism,
                    reported_organism,
                    accepted_names=accepted_names,
                ),
                reason=(
                    f"{contamination_reason} The primary assembly signal still "
                    "passed the species ANI and aligned-fraction thresholds."
                ),
                closest_species=closest,
                kmerfinder_support=tuple(kmerfinder_support),
                mixed_signal=mixed_signal,
                identification_thresholds=thresholds,
                contamination_evidence=contamination_evidence,
            )
        return _result(
            sample_id=sample_id,
            status="suspected_culture_contamination",
            expected_organism=expected_organism,
            reason=contamination_reason,
            closest=closest,
            kmerfinder_support=kmerfinder_support,
            mixed_signal=mixed_signal,
            identification_thresholds=thresholds,
            contamination_evidence=contamination_evidence,
        )
    if mixed_signal.status in {"indeterminate", "minor_secondary"}:
        reason = (
            "A supported secondary bacterial signal did not meet the stricter "
            "criteria for suspected culture contamination."
            if mixed_signal.status == "minor_secondary"
            else "Possible secondary species signals require review before a species call."
        )
        return _result(
            sample_id=sample_id,
            status="mixed_signal_review",
            expected_organism=expected_organism,
            reason=reason,
            closest=closest,
            kmerfinder_support=kmerfinder_support,
            mixed_signal=mixed_signal,
            identification_thresholds=thresholds,
            contamination_evidence=contamination_evidence,
        )
    if not passing:
        reason = (
            "No reference passed both the whole-genome ANI and bidirectional "
            "aligned-fraction thresholds."
        )
        if not closest:
            reason = "No reference covered enough of the assembly to report a closest species."
        return _result(
            sample_id=sample_id,
            status="no_confident_identification",
            expected_organism=expected_organism,
            reason=reason,
            closest=closest,
            kmerfinder_support=kmerfinder_support,
            mixed_signal=mixed_signal,
            identification_thresholds=thresholds,
            contamination_evidence=contamination_evidence,
        )
    if competitor is not None:
        return _result(
            sample_id=sample_id,
            status="ambiguous_species",
            expected_organism=expected_organism,
            reason=(
                f"The top two species differ by less than {policy.minimum_ani_margin:.2f} "
                "ANI percentage points without a sufficient bidirectional "
                "aligned-fraction advantage."
            ),
            closest=closest,
            kmerfinder_support=kmerfinder_support,
            mixed_signal=mixed_signal,
            identification_thresholds=thresholds,
            contamination_evidence=contamination_evidence,
        )
    assert report_hit is not None
    reported_organism = report_hit.species_name
    taxonomic_confirmation = _taxonomic_confirmation(
        expected_organism,
        reported_organism,
        accepted_names=accepted_names,
    )
    return IdentificationResult(
        sample_id=sample_id,
        status="confirmed_species",
        reported_organism=reported_organism,
        match_percent=report_hit.ani,
        taxonomic_confirmation=taxonomic_confirmation,
        reason=(
            "The top reference passed the ANI and aligned-fraction thresholds "
            "without an equally plausible independent competing organism."
        ),
        closest_species=closest,
        kmerfinder_support=tuple(kmerfinder_support),
        mixed_signal=mixed_signal,
        identification_thresholds=thresholds,
        contamination_evidence=contamination_evidence,
    )


def _select_user_facing_hit(
    *,
    top: AniHit,
    ani_hits: Sequence[AniHit],
    expected_organism: str,
    accepted_group_members: Iterable[str],
    policy: IdentificationPolicy,
) -> AniHit:
    """Keep reporting groups internal while preserving clinical terminology."""

    if not top.reporting_group:
        return top
    expected = expected_organism_name(expected_organism)
    accepted = {expected_organism_name(name) for name in accepted_group_members}
    if expected in UNKNOWN_NAMES or expected not in accepted:
        return top
    expected_hits = [
        hit
        for hit in ani_hits
        if hit.reporting_group == top.reporting_group
        and expected_organism_name(hit.species_name) == expected
        and hit.ani >= policy.minimum_ani
        and policy.aligned_fraction(hit) >= policy.minimum_aligned_fraction
    ]
    if not expected_hits:
        return top
    return max(
        expected_hits,
        key=lambda hit: (hit.ani, policy.aligned_fraction(hit)),
    )


def _result(
    *,
    sample_id: str,
    status: str,
    expected_organism: str,
    reason: str,
    closest: tuple[AniHit, ...],
    kmerfinder_support: Sequence[KmerSpeciesEvidence],
    mixed_signal: MixedSignalAssessment,
    identification_thresholds: dict[str, object],
    contamination_evidence: ContaminationEvidence,
) -> IdentificationResult:
    expected_known = _normalise_name(expected_organism) not in UNKNOWN_NAMES
    return IdentificationResult(
        sample_id=sample_id,
        status=status,
        reported_organism="No confident identification",
        match_percent=None,
        taxonomic_confirmation=(
            "Fail"
            if expected_known or status == "suspected_culture_contamination"
            else "N/A"
        ),
        reason=reason,
        closest_species=closest,
        kmerfinder_support=tuple(kmerfinder_support),
        mixed_signal=mixed_signal,
        identification_thresholds=identification_thresholds,
        contamination_evidence=contamination_evidence,
    )


def _best_hit_per_species(
    hits: Sequence[AniHit], *, policy: IdentificationPolicy
) -> list[AniHit]:
    grouped: dict[str, list[AniHit]] = {}
    for hit in hits:
        key = (
            f"reporting_group:{_normalise_name(hit.reporting_group)}"
            if hit.reporting_group
            else f"species_taxid:{hit.species_taxid}"
        )
        grouped.setdefault(key, []).append(hit)
    by_species: dict[str, AniHit] = {}
    for group_key, species_hits in grouped.items():
        adequately_aligned = [
            hit
            for hit in species_hits
            if policy.aligned_fraction(hit) >= policy.minimum_aligned_fraction
        ]
        if adequately_aligned:
            selected = max(
                adequately_aligned,
                key=lambda hit: (hit.ani, policy.aligned_fraction(hit)),
            )
        else:
            selected = max(
                species_hits,
                key=lambda hit: (policy.aligned_fraction(hit), hit.ani),
            )
        by_species[group_key] = selected
    def rank(hit: AniHit) -> tuple[float, float, float]:
        if (
            hit.ani >= policy.minimum_ani
            and policy.aligned_fraction(hit) >= policy.minimum_aligned_fraction
        ):
            return (2.0, hit.ani, policy.aligned_fraction(hit))
        return (1.0, policy.aligned_fraction(hit), hit.ani)

    return sorted(by_species.values(), key=rank, reverse=True)


def _taxonomic_confirmation(
    expected: str,
    observed: str,
    *,
    accepted_names: Iterable[str] = (),
) -> str:
    if _normalise_name(expected) in UNKNOWN_NAMES:
        return "N/A"
    accepted = {_species_level_expectation(observed)}
    accepted.update(_species_level_expectation(name) for name in accepted_names)
    return "Pass" if _species_level_expectation(expected) in accepted else "Fail"


def _species_level_expectation(value: object) -> str:
    """Use the same conservative comparison as QC input and dashboard checks."""
    return expected_organism_name(value)


def run_skani_search(
    *,
    assembly: str | Path,
    database: str | Path,
    output_path: str | Path,
    threads: int = 1,
    executable: str = "skani",
) -> None:
    assembly_path = Path(assembly)
    database_path = Path(database)
    output = Path(output_path)
    if not assembly_path.is_file() or assembly_path.stat().st_size == 0:
        raise ValueError(f"Assembly is missing or empty: {assembly_path}")
    if not database_path.exists():
        raise ValueError(f"skani database does not exist: {database_path}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    command = [
        executable,
        "search",
        str(assembly_path),
        "-d",
        str(database_path),
        "-o",
        str(temporary),
        "-t",
        str(max(1, int(threads))),
    ]
    try:
        completed = subprocess.run(command, text=True, capture_output=True)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"skani search could not be started: {exc}") from exc
    if completed.returncode != 0:
        temporary.unlink(missing_ok=True)
        detail = (completed.stderr or completed.stdout or "unknown error").strip()
        raise RuntimeError(f"skani search failed with exit {completed.returncode}: {detail}")
    if not temporary.is_file():
        raise RuntimeError("skani search completed without producing its result table")
    os.replace(temporary, output)


def write_identification_outputs(
    result: IdentificationResult, output_directory: str | Path
) -> None:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    _atomic_text(
        output / "species_identification.json",
        json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
    )
    summary_fields = [
        "sample_id",
        "status",
        "reported_organism",
        "match_percent",
        "taxonomic_confirmation",
        "reason",
    ]
    _write_tsv(
        output / "species_identification.tsv",
        summary_fields,
        [
            {
                "sample_id": result.sample_id,
                "status": result.status,
                "reported_organism": result.reported_organism,
                "match_percent": (
                    f"{result.match_percent:.2f}" if result.match_percent is not None else "N/A"
                ),
                "taxonomic_confirmation": result.taxonomic_confirmation,
                "reason": result.reason,
            }
        ],
    )
    cross_sample_path = output / "cross_sample_signals.tsv"
    if result.cross_sample_evidence:
        cross_sample_fields = [
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
        _write_tsv(
            cross_sample_path,
            cross_sample_fields,
            [
                {
                    **asdict(evidence),
                    "fraction_of_source": f"{evidence.fraction_of_source:.6f}",
                    "fraction_of_sample_primary": (
                        f"{evidence.fraction_of_sample_primary:.6f}"
                    ),
                }
                for evidence in result.cross_sample_evidence
            ],
        )
    else:
        cross_sample_path.unlink(missing_ok=True)
    closest_fields = [
        "rank",
        "species_taxid",
        "species_name",
        "ani",
        "aligned_fraction_reference",
        "aligned_fraction_query",
        "reference_accession",
        "reporting_group",
    ]
    _write_tsv(
        output / "closest_species.tsv",
        closest_fields,
        [
            {
                "rank": rank,
                "species_taxid": hit.species_taxid,
                "species_name": hit.species_name,
                "ani": f"{hit.ani:.2f}",
                "aligned_fraction_reference": f"{hit.aligned_fraction_reference:.2f}",
                "aligned_fraction_query": f"{hit.aligned_fraction_query:.2f}",
                "reference_accession": hit.accession,
                "reporting_group": hit.reporting_group,
            }
            for rank, hit in enumerate(result.closest_species, start=1)
        ],
    )


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _atomic_text(path: Path, text: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
