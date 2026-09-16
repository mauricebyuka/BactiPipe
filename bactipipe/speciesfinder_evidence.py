"""SpeciesFinder candidate nomination with fail-closed result validation."""

from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Sequence

from bactipipe.species_identification import TaxonomicSignal


@dataclass(frozen=True)
class SpeciesFinderHit:
    """One SpeciesFinder/KMA template hit with explicit metric semantics."""

    assembly_accession: str
    template_accession: str
    description: str
    organism_taxid: str
    taxonomy: str
    species_taxid: str
    species_name: str
    query_count: int
    score: float
    query_identity: float
    query_coverage: float
    template_coverage: float
    depth: float


def run_speciesfinder(
    *,
    inputs: Sequence[str | Path],
    output_directory: str | Path,
    database_prefix: str | Path,
    taxonomy_file: str | Path,
    executable: str = "speciesfinder",
) -> tuple[SpeciesFinderHit, ...]:
    """Run SpeciesFinder and reject silent KMA failures or malformed outputs."""

    if len(inputs) not in {1, 2}:
        raise ValueError("SpeciesFinder requires one input or one paired-end read pair")
    input_paths = [Path(value).resolve() for value in inputs]
    for path in input_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"SpeciesFinder input is missing or empty: {path}")
    prefix = Path(database_prefix).resolve()
    taxonomy = Path(taxonomy_file).resolve()
    for suffix in (".comp.b", ".seq.b", ".length.b", ".name"):
        candidate = Path(f"{prefix}{suffix}")
        if not candidate.is_file() or candidate.stat().st_size == 0:
            raise ValueError(f"SpeciesFinder KMA database component is missing: {candidate}")
    if not taxonomy.is_file() or taxonomy.stat().st_size == 0:
        raise ValueError(f"SpeciesFinder taxonomy file is missing or empty: {taxonomy}")

    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # A retry may target an existing failed work directory. Remove only the
    # exact files owned by SpeciesFinder so stale success artifacts cannot make
    # a later incomplete invocation appear valid.
    for name in (
        "results.txt",
        "results.res",
        "output.json",
        "speciesfinder.log",
        "speciesfinder.stdout.log",
        "speciesfinder.stderr.log",
    ):
        (output / name).unlink(missing_ok=True)
    command = [
        executable,
        "-i",
        *[str(path) for path in input_paths],
        "-o",
        str(output),
        "-db",
        str(prefix),
        "-x",
        "-tax",
        str(taxonomy),
    ]
    try:
        completed = subprocess.run(command, text=True, capture_output=True)
    except OSError as exc:
        raise RuntimeError(f"SpeciesFinder could not be started: {exc}") from exc
    (output / "speciesfinder.stdout.log").write_text(completed.stdout or "")
    (output / "speciesfinder.stderr.log").write_text(completed.stderr or "")
    if completed.returncode != 0:
        raise RuntimeError(
            f"SpeciesFinder failed with exit {completed.returncode}; see "
            f"{output / 'speciesfinder.stderr.log'}"
        )

    internal_log = output / "speciesfinder.log"
    log_text = internal_log.read_text(errors="replace") if internal_log.is_file() else ""
    if "KMA failed with return code" in log_text:
        raise RuntimeError(
            "SpeciesFinder reported an internal KMA failure despite returning success; "
            f"see {internal_log}"
        )
    results = output / "results.txt"
    raw_results = output / "results.res"
    json_result = output / "output.json"
    for required in (results, raw_results, json_result):
        if not required.is_file() or required.stat().st_size == 0:
            raise RuntimeError(
                f"SpeciesFinder completed without a required result: {required}"
            )
    try:
        json.loads(json_result.read_text())
    except (OSError, ValueError, TypeError) as exc:
        raise RuntimeError(f"SpeciesFinder produced invalid JSON: {json_result}") from exc
    return parse_speciesfinder_results(results)


def parse_speciesfinder_results(path: str | Path) -> tuple[SpeciesFinderHit, ...]:
    """Parse extended TSV output without trusting its repackaged JSON aliases."""

    source = Path(path)
    with source.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            headers = [_normalise_header(value) for value in next(reader)]
        except StopIteration as exc:
            raise ValueError("SpeciesFinder output is empty") from exc
        # SpeciesFinder 4.0.4 hard-codes both the template-organism and species
        # taxid headers as "TAXID". Preserve their positional semantics instead
        # of allowing csv.DictReader to overwrite the first value.
        taxid_positions = [
            index for index, header in enumerate(headers) if header == "taxid"
        ]
        if "taxid species" not in headers and len(taxid_positions) >= 2:
            headers[taxid_positions[-1]] = "taxid species"
        fields = set(headers)
        required = {
            "assembly",
            "accession number",
            "species",
            "taxid species",
            "score",
            "query coverage",
            "template coverage",
        }
        missing = required.difference(fields)
        if missing:
            raise ValueError(
                "SpeciesFinder output is missing columns: " + ", ".join(sorted(missing))
            )
        parsed: list[SpeciesFinderHit] = []
        for line_number, values in enumerate(reader, start=2):
            if not values or not any(value.strip() for value in values):
                continue
            if len(values) != len(headers):
                raise ValueError(
                    f"SpeciesFinder output line {line_number} has {len(values)} "
                    f"columns; expected {len(headers)}"
                )
            row = dict(zip(headers, values))
            description = _clean(row.get("description"))
            if is_excluded_reference_source(description):
                continue
            species_name, species_taxid = _clinically_relevant_taxon(
                species_name=_clean(row.get("species")),
                species_taxid=_clean(row.get("taxid species")),
                organism_taxid=_clean(row.get("taxid")),
                taxonomy=_clean(row.get("taxonomy")),
                description=description,
            )
            template_accession = _clean(row.get("accession number")).split(" ", 1)[0]
            if (
                not species_name
                or not species_taxid.isdigit()
                or not template_accession
                or not is_reportable_species_name(species_name)
            ):
                continue
            parsed.append(
                SpeciesFinderHit(
                    assembly_accession=_clean(row.get("assembly")),
                    template_accession=template_accession,
                    description=description,
                    organism_taxid=_clean(row.get("taxid")),
                    taxonomy=_clean(row.get("taxonomy")),
                    species_taxid=species_taxid,
                    species_name=species_name,
                    query_count=_integer(
                        row.get("num") or row.get("number"),
                        default=0,
                        field=f"query count line {line_number}",
                    ),
                    score=_number(
                        row.get("score"), field=f"score line {line_number}"
                    ),
                    query_identity=_bounded_percent(
                        row.get("query identity"),
                        default=0.0,
                        field=f"query identity line {line_number}",
                    ),
                    # KMA's query coverage is a query/template length ratio and
                    # can exceed 100 for a multi-contig assembly. Preserve the
                    # value under its upstream name; it is nomination evidence,
                    # not BactiPipe's whole-genome aligned fraction.
                    query_coverage=_nonnegative_number(
                        row.get("tot query coverage") or row.get("query coverage"),
                        field=f"query coverage line {line_number}",
                    ),
                    # KMA defines this as 100 * aligned_length / template_length.
                    # Insertions can therefore make it exceed 100. Like its
                    # reciprocal Query_Coverage, it is nomination evidence and
                    # not a conventional reference-breadth percentage.
                    template_coverage=_nonnegative_number(
                        row.get("tot template coverage")
                        or row.get("template coverage"),
                        field=f"template coverage line {line_number}",
                    ),
                    depth=_number(
                        row.get("tot depth") or row.get("depth"),
                        default=0.0,
                        field=f"depth line {line_number}",
                    ),
                )
            )
    parsed.sort(
        key=lambda hit: (
            hit.query_count,
            hit.score,
            hit.template_coverage,
            hit.query_identity,
        ),
        reverse=True,
    )
    return tuple(parsed)


def select_candidate_hits(
    hits: Sequence[SpeciesFinderHit],
    *,
    maximum_species: int = 20,
    references_per_species: int = 3,
) -> tuple[SpeciesFinderHit, ...]:
    """Keep a bounded, diverse candidate set without making an identification."""

    if maximum_species < 1 or references_per_species < 1:
        raise ValueError("SpeciesFinder candidate limits must be positive")
    selected: list[SpeciesFinderHit] = []
    species_order: list[str] = []
    counts: dict[str, int] = {}
    seen_templates: set[str] = set()
    for hit in hits:
        if hit.template_accession in seen_templates:
            continue
        if hit.species_taxid not in counts:
            if len(species_order) >= maximum_species:
                continue
            species_order.append(hit.species_taxid)
            counts[hit.species_taxid] = 0
        if counts[hit.species_taxid] >= references_per_species:
            continue
        selected.append(hit)
        seen_templates.add(hit.template_accession)
        counts[hit.species_taxid] += 1
    return tuple(selected)


def write_speciesfinder_hits(
    hits: Sequence[SpeciesFinderHit], path: str | Path
) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    fieldnames = [field.name for field in fields(SpeciesFinderHit)]
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(hit.__dict__)
    os.replace(temporary, destination)


def speciesfinder_hits_to_signals(
    hits: Sequence[SpeciesFinderHit],
) -> tuple[TaxonomicSignal, ...]:
    """Convert nominations to preliminary signals for competitive remapping.

    SpeciesFinder/KMA does not expose Kraken-style distinct minimizers. Its
    score is used only to rank candidates before competitive mapping; final
    mixture decisions use mapped-read counts and reference breadth.
    """

    grouped: dict[tuple[str, str], dict[str, float | int]] = {}
    for hit in hits:
        key = (hit.species_taxid, hit.species_name)
        item = grouped.setdefault(
            key,
            {"support": 0, "score": 0.0},
        )
        item["support"] = int(item["support"]) + max(hit.query_count, 0)
        item["score"] = float(item["score"]) + max(hit.score, 0.0)
    total_support = sum(int(item["support"]) for item in grouped.values())
    total_score = sum(max(0.0, float(item["score"])) for item in grouped.values())
    signals = []
    for (taxid, name), item in grouped.items():
        support = int(item["support"])
        score = max(0.0, float(item["score"]))
        signals.append(
            TaxonomicSignal(
                species_taxid=taxid,
                species_name=name,
                # SpeciesFinder/KMA does not report a read count in its real
                # extended output. Keep an optional upstream count when a
                # compatible table supplies one, but never relabel Score as
                # reads. Competitive mapping supplies the diagnostic counts.
                read_count=support,
                read_percent=(
                    100.0 * support / total_support if total_support else 0.0
                ),
                distinct_minimizers=None,
                nomination_score=score,
                nomination_percent=(
                    100.0 * score / total_score if total_score else 0.0
                ),
            )
        )
    signals.sort(
        key=lambda signal: (
            signal.nomination_score or 0.0,
            signal.read_count,
        ),
        reverse=True,
    )
    return tuple(signals)


def _clinically_relevant_taxon(
    *,
    species_name: str,
    species_taxid: str,
    organism_taxid: str,
    taxonomy: str,
    description: str,
) -> tuple[str, str]:
    """Retain explicitly controlled infraspecific diagnostic taxa."""

    evidence = " ".join((taxonomy, description)).casefold()
    controlled = (
        (
            "streptococcus equi subsp. zooepidemicus",
            "40041",
            "Streptococcus equi subsp. zooepidemicus",
        ),
        (
            "streptococcus equi subsp. equi",
            "148942",
            "Streptococcus equi subsp. equi",
        ),
    )
    for phrase, taxid, name in controlled:
        if phrase in evidence:
            return name, taxid
    return species_name, species_taxid or organism_taxid


def is_excluded_reference_source(description: str) -> bool:
    """Reject non-isolate templates that cannot support a diagnostic call."""

    text = description.casefold()
    return bool(
        "metagenome-assembled" in text
        or "metagenomic assembled" in text
        or "single-cell" in text
        or re.search(r"(?:^|[^a-z0-9])(?:mag|sag)(?:[^a-z0-9]|$)", text)
    )


def is_reportable_species_name(species_name: str) -> bool:
    """Return whether a taxonomy label can support a species-level call."""

    text = " ".join(species_name.casefold().replace(".", "").split())
    if not text or text in {"na", "n/a", "unknown"}:
        return False
    if any(
        phrase in text
        for phrase in (
            "uncultured",
            "unclassified",
            "candidate division",
        )
    ):
        return False
    if re.search(
        r"(?:^|\s)(?:sp|bacterium|archaeon|symbiont|endosymbiont)(?:\s|$)",
        text,
    ):
        return False
    tokens = text.split()
    minimum_tokens = 3 if tokens[0] == "candidatus" else 2
    return len(tokens) >= minimum_tokens


def _normalise_header(value: object) -> str:
    return " ".join(str(value or "").lstrip("#").replace("_", " ").split()).casefold()


def _clean(value: object) -> str:
    return " ".join(str(value or "").split())


def _number(value: object, *, field: str, default: float | None = None) -> float:
    text = str(value or "").strip()
    if not text and default is not None:
        return default
    try:
        return float(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid SpeciesFinder {field}: {value!r}") from exc


def _integer(value: object, *, field: str, default: int | None = None) -> int:
    text = str(value or "").strip()
    if not text and default is not None:
        return default
    try:
        return int(float(text))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid SpeciesFinder {field}: {value!r}") from exc


def _bounded_percent(
    value: object, *, field: str, default: float | None = None
) -> float:
    result = _number(value, field=field, default=default)
    if not 0.0 <= result <= 100.0:
        raise ValueError(f"Out-of-range SpeciesFinder {field}: {result}")
    return result


def _nonnegative_number(value: object, *, field: str) -> float:
    result = _number(value, field=field)
    if result < 0.0:
        raise ValueError(f"Negative SpeciesFinder {field}: {result}")
    return result
