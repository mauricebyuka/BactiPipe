"""Controlled organism-name and genome-size resolution for BactiPipe QC.

The controlled vocabulary remains ``data/pathogenic_bacteria.txt``.  Its size
column is a reviewed reference-record size, not an estimate supplied by the
operator and not a species-wide statistical median.  Keeping resolution in one
small module prevents runners and multiprocessing workers from silently using
different fallback sizes.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from importlib.resources import files as resource_files
from pathlib import Path
from .organism_names import expected_organism_name, normalize_organism_name


CATALOG_VERSION = "BactiPipe-genome-sizes-v1"
CATALOG_RESOURCE = "pathogenic_bacteria.txt"
UNKNOWN_NAME = "Unknown"


@dataclass(frozen=True)
class GenomeSizeRecord:
    taxid: str
    genome_size_bp: int
    reference_accession: str
    canonical_name: str
    source: str = "controlled_reference_record"


@dataclass(frozen=True)
class GenomeSizeResolution:
    entered_name: str
    status: str
    canonical_name: str | None = None
    taxid: str | None = None
    genome_size_bp: int | None = None
    reference_accession: str | None = None
    source: str | None = None
    catalog_version: str = CATALOG_VERSION

    @property
    def is_known(self) -> bool:
        return self.status == "known"

    @property
    def is_unknown(self) -> bool:
        return self.status in {"unknown", "unsupported", "missing"}


class GenomeSizeCatalog:
    """Case-insensitive resolver over the controlled organism vocabulary."""

    def __init__(self, records: list[GenomeSizeRecord], checksum: str):
        self.records = tuple(records)
        self.checksum = checksum
        self._by_normalized_name: dict[str, GenomeSizeRecord] = {}
        for record in self.records:
            key = normalize_organism_name(record.canonical_name)
            if not key:
                raise ValueError("Controlled organism names cannot be blank.")
            if key in self._by_normalized_name:
                raise ValueError(
                    f"Duplicate normalized controlled organism name: {record.canonical_name}"
                )
            if record.genome_size_bp <= 0:
                raise ValueError(
                    f"Genome size must be positive for {record.canonical_name}."
                )
            self._by_normalized_name[key] = record

        # Lambda is an internal process-control genome and is intentionally not
        # part of the veterinary bacterial organism list.
        self._lambda = GenomeSizeRecord(
            taxid="10710",
            genome_size_bp=48_502,
            reference_accession="bundled:Lambda",
            canonical_name="Lambda",
            source="bundled_process_control",
        )

    @classmethod
    def from_path(cls, path: str | Path) -> "GenomeSizeCatalog":
        path = Path(path)
        raw = path.read_bytes()
        records: list[GenomeSizeRecord] = []
        for line_number, raw_line in enumerate(
            raw.decode("utf-8").splitlines(), start=1
        ):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"Invalid controlled organism row {line_number}: expected four columns."
                )
            taxid, size, accession = parts[:3]
            name = parts[-1].strip()
            records.append(
                GenomeSizeRecord(
                    taxid=taxid.strip(),
                    genome_size_bp=int(size),
                    reference_accession=accession.strip(),
                    canonical_name=name,
                )
            )
        return cls(records, sha256(raw).hexdigest())

    @classmethod
    def from_package(cls) -> "GenomeSizeCatalog":
        resource = resource_files("bactipipe.data").joinpath(CATALOG_RESOURCE)
        with resource.open("rb") as handle:
            raw = handle.read()
        records: list[GenomeSizeRecord] = []
        for line_number, raw_line in enumerate(
            raw.decode("utf-8").splitlines(), start=1
        ):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                raise ValueError(
                    f"Invalid controlled organism row {line_number}: expected four columns."
                )
            taxid, size, accession = parts[:3]
            records.append(
                GenomeSizeRecord(
                    taxid=taxid.strip(),
                    genome_size_bp=int(size),
                    reference_accession=accession.strip(),
                    canonical_name=parts[-1].strip(),
                )
            )
        return cls(records, sha256(raw).hexdigest())

    def resolve(self, entered_name: object) -> GenomeSizeResolution:
        entered = str(entered_name or "").strip()
        normalized = normalize_organism_name(entered)
        if normalized in {"unknown", "organism unknown"}:
            return GenomeSizeResolution(
                entered_name=entered,
                status="unknown",
                canonical_name=UNKNOWN_NAME,
                source="operator_unspecified",
            )
        if normalized == "lambda":
            record = self._lambda
        else:
            record = self._by_normalized_name.get(normalized)
            if record is None:
                record = self._by_normalized_name.get(expected_organism_name(entered))
        if record is None:
            return GenomeSizeResolution(
                entered_name=entered,
                status="unsupported" if normalized else "missing",
                canonical_name=UNKNOWN_NAME,
                source="unrecognized_input" if normalized else "missing_input",
            )
        return GenomeSizeResolution(
            entered_name=entered,
            status="known",
            canonical_name=record.canonical_name,
            taxid=record.taxid,
            genome_size_bp=record.genome_size_bp,
            reference_accession=record.reference_accession,
            source=record.source,
        )


def coverage_and_verdict(
    total_bases: int,
    resolution: GenomeSizeResolution,
    minimum_depth: float,
) -> tuple[float | None, str]:
    """Calculate coverage without treating unavailable coverage as a pass."""

    if not resolution.genome_size_bp:
        return None, "N/A"
    coverage = total_bases / resolution.genome_size_bp
    return coverage, "Pass" if coverage >= minimum_depth else "Fail"
