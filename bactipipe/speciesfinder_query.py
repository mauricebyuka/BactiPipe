"""Lossless assembly fragments for SpeciesFinder nomination, never for ANI."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


FRAGMENT_LIMIT_BP = 20_000
QUERY_METHOD = "nonoverlapping_assembly_fragments_v1"


def prepare_nomination_query(
    assembly: str | Path, output_directory: str | Path,
    *, limit_bp: int = FRAGMENT_LIMIT_BP,
) -> dict:
    """Write a deterministic query plus compact audit in a NEW directory.

    Every contig is represented, including short contigs and final tails. No
    sequence trimming, overlap, rotation or reverse complementation is applied.
    Coordinates are zero-based, half-open; numeric IDs tolerate duplicate source
    headers. The caller may remove query.fasta after successful nomination, but
    must retain query.json and fragments.tsv. A failed caller owns this directory.
    """
    if isinstance(limit_bp, bool) or not isinstance(limit_bp, int) or limit_bp < 1:
        raise ValueError("Fragment limit must be a positive integer")
    source = Path(assembly).resolve(strict=True)
    output = Path(output_directory).resolve()
    # Never overwrite an assembly, previous audit or another concurrent attempt.
    output.mkdir(parents=True, exist_ok=False)
    query = output / "query.fasta"
    source_hash = hashlib.sha256()
    sequence_hash = hashlib.sha256()
    query_hash = hashlib.sha256()
    contigs = fragments = bases = 0
    header = None
    pieces = []

    with query.open("xb") as fasta, (output / "fragments.tsv").open("x", newline="") as audit:
        writer = csv.writer(audit, delimiter="\t", lineterminator="\n")
        writer.writerow(("fragment_id", "contig_number", "source_header", "contig_length_bp",
                         "start_0based", "end_exclusive", "length_bp"))

        def emit():
            nonlocal contigs, fragments, bases
            if header is None:
                return
            sequence = "".join(pieces)
            if not sequence:
                raise ValueError(f"Empty assembly FASTA record: {header}")
            contigs += 1
            encoded = sequence.encode("ascii")
            sequence_hash.update(encoded)
            bases += len(sequence)
            for start in range(0, len(sequence), limit_bp):
                stop = min(start + limit_bp, len(sequence))
                name = f"c{contigs}_p{start}"
                record = f">{name}\n{sequence[start:stop]}\n".encode("ascii")
                fasta.write(record)
                query_hash.update(record)
                writer.writerow((name, contigs, header, len(sequence), start, stop, stop - start))
                fragments += 1

        with source.open("rb") as handle:
            for raw in handle:
                source_hash.update(raw)
                line = raw.decode("utf-8").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    emit()
                    header = line[1:].strip()
                    if not header:
                        raise ValueError("Assembly FASTA has an empty header")
                    pieces = []
                else:
                    if header is None:
                        raise ValueError("Assembly FASTA sequence precedes its header")
                    sequence = "".join(line.split())
                    if set(sequence.upper()) - set("ACGTRYSWKMBDHVNU"):
                        raise ValueError(f"Invalid assembly nucleotide sequence in {header}")
                    pieces.append(sequence)
        emit()
        if not contigs:
            raise ValueError("Assembly FASTA contains no sequences")

    evidence = {
        "schema_version": 1,
        "method": QUERY_METHOD,
        "limit_bp": limit_bp,
        "source_filename": source.name,
        "source_sha256": source_hash.hexdigest(),
        "query_filename": query.name,
        "query_sha256": query_hash.hexdigest(),
        "concatenated_sequence_sha256": sequence_hash.hexdigest(),
        "contig_count": contigs,
        "fragment_count": fragments,
        "total_bases": bases,
        "coordinate_system": "zero_based_half_open",
        "coordinates_file": "fragments.tsv",
        "ani_input": "original_assembly",
    }
    (output / "query.json").write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    return evidence
