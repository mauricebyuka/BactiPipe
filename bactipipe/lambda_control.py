"""Explicit Lambda control handling, separate from bacterial identification."""

import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from datetime import datetime, timezone

REFERENCE_SEQUENCE_SHA256 = "135c9042edbda8729c0bb550e2747fe968383904344123fad3e71a17b8559383"
MIN_IDENTITY = 99.0
MIN_COVERAGE = 90.0

CONTROL_IDENTIFICATION = "Control confirmation not assessed"


def is_lambda_control(expected_organism):
    """Use the declared organism, never infer a control from its sample filename."""
    return str(expected_organism or "").strip().casefold() == "lambda"


def control_summary_row(sample, quality, quality_verdict, organism, depth, minimum_depth, depth_verdict, evidence=None):
    evidence = evidence or {}
    status = evidence.get("status")
    hit = "Lambda" if status == "Pass" else "Lambda not confirmed" if status == "Fail" else CONTROL_IDENTIFICATION
    identity = evidence.get("identity_percent")
    return [sample, f"{quality:.2f}", quality_verdict, organism, hit,
            f"{identity:.6f}" if identity is not None else "N/A", depth, minimum_depth,
            depth_verdict, "Control" + status if status in {"Pass", "Fail", "Error"} else "Review"]


def percent_threshold(value):
    value = float(value)
    if not math.isfinite(value) or not 0 < value <= 100:
        raise ValueError("Control thresholds must be greater than zero and at most 100.")
    return value


def add_control_arguments(parser):
    parser.add_argument("--lambda-min-identity", type=percent_threshold, default=MIN_IDENTITY,
                        help="Minimum Lambda reference-alignment identity percent (default: 99).")
    parser.add_argument("--lambda-min-coverage", type=percent_threshold, default=MIN_COVERAGE,
                        help="Minimum Lambda reference genome coverage percent (default: 90).")


def alignment_metrics(lines, reference_length):
    """Count each reference position once, assigning overlapping HSPs by best score.

    Identity includes mismatches and alignment gaps. Coverage counts reference bases
    paired to unambiguous query bases, not deletions. Insertions are assigned to a
    neighbouring reference position and counted only with that position's HSP.
    """
    hsps = []
    for line in lines:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 10:
            raise ValueError("Invalid Lambda BLAST alignment record.")
        qid, sid, qs, qe, ss, se, evalue, score, qseq, sseq = fields
        start, end = int(ss), int(se)
        if not 1 <= start <= reference_length or not 1 <= end <= reference_length or len(qseq) != len(sseq):
            raise ValueError("Invalid Lambda alignment coordinates.")
        if sum(c != '-' for c in sseq) != abs(end-start)+1:
            raise ValueError("Lambda alignment reference span does not match its sequence.")
        score = float(score)
        if not math.isfinite(score):
            raise ValueError("Invalid Lambda alignment score.")
        hsps.append((score, qid, start, end, qseq.upper(), sseq.upper()))
    assigned = {}
    for score, qid, start, end, query, subject in sorted(hsps, key=lambda x: (-x[0], x[1:])):
        step = 1 if end >= start else -1
        pos = start
        columns = {}
        insertions = {}
        for q, s in zip(query, subject):
            if s == '-':
                anchor = pos if min(start,end) <= pos <= max(start,end) else end
                insertions[anchor] = insertions.get(anchor, 0) + 1
                continue
            callable_base = q in 'ACGT' and s in 'ACGT'
            columns[pos] = (int(callable_base and q == s), int(callable_base))
            pos += step
        for pos, (match, covered) in columns.items():
            if pos not in assigned:
                assigned[pos] = (match, covered, 1 + insertions.get(pos, 0))
    matches = sum(x[0] for x in assigned.values())
    covered = sum(x[1] for x in assigned.values())
    length = sum(x[2] for x in assigned.values())
    return {"identity_percent": 100 * matches / length if length else None,
            "reference_coverage_percent": 100 * covered / reference_length,
            "covered_reference_bases": covered, "matching_bases": matches,
            "alignment_columns": length, "alignment_count": len(hsps)}


def assess_control(assembly, output_dir, *, min_identity=MIN_IDENTITY, min_coverage=MIN_COVERAGE, report=lambda message: None):
    """Align an assembly to the pinned bundled control; retain small evidence files."""
    min_identity, min_coverage = percent_threshold(min_identity), percent_threshold(min_coverage)
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 1, "assessed_at": datetime.now(timezone.utc).isoformat(),
              "method": "BLASTn assembly-to-control; unique reference-position HSP selection",
              "minimum_identity_percent": min_identity, "minimum_coverage_percent": min_coverage,
              "status": "Error", "identity_percent": None, "reference_coverage_percent": None}
    try:
        reference = Path(__file__).parent / "data" / "lambda.fasta"
        payload = reference.read_bytes()
        sequence = ''.join(line.strip() for line in payload.decode().splitlines() if not line.startswith('>')).upper()
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if len(sequence) != 48502 or digest != REFERENCE_SEQUENCE_SHA256:
            raise ValueError("Bundled Lambda reference failed length/checksum validation.")
        result["reference"] = {"name": "Bundled Lambda control", "accession": None, "length_bp": len(sequence),
                               "sequence_sha256": digest, "fasta_sha256": hashlib.sha256(payload).hexdigest(),
                               "ncbi_comparison_accession": "NC_001416.1", "ncbi_comparison_differences": 4}
        assembly = Path(assembly)
        result["assembly_sha256"] = hashlib.sha256(assembly.read_bytes()).hexdigest()
        version = subprocess.run(["blastn", "-version"], capture_output=True, text=True, timeout=30, check=True)
        result["tool_version"] = version.stdout.splitlines()[0]
        command = ["blastn", "-query", str(assembly), "-subject", str(reference), "-task", "blastn",
                   "-dust", "no", "-soft_masking", "false", "-evalue", "1e-10",
                   "-outfmt", "6 qseqid sseqid qstart qend sstart send evalue bitscore qseq sseq"]
        result["command"] = command
        with (root / "alignments.tsv").open('w') as out, (root / "alignment.log").open('w') as err:
            subprocess.run(command, stdout=out, stderr=err, timeout=120, check=True)
        with (root / "alignments.tsv").open() as handle:
            result.update(alignment_metrics(handle, len(sequence)))
        identity = result["identity_percent"]
        result["identity_pass"] = identity is not None and identity >= min_identity
        result["coverage_pass"] = result["reference_coverage_percent"] >= min_coverage
        result["status"] = "Pass" if result["identity_pass"] and result["coverage_pass"] else "Fail"
        result["reason"] = "Lambda confirmed against the bundled control reference." if result["status"] == "Pass" else "Lambda reference identity and/or genome coverage did not meet the control thresholds."
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        result["reason"] = f"Lambda reference confirmation could not be completed: {exc}"
    temporary = root / "result.json.tmp"
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    os.replace(temporary, root / "result.json")
    with (root / "metrics.tsv").open('w') as handle:
        keys = ['status','identity_percent','reference_coverage_percent','covered_reference_bases',
                'minimum_identity_percent','minimum_coverage_percent','reason']
        writer = csv.DictWriter(handle,fieldnames=keys,delimiter='\t',extrasaction='ignore')
        writer.writeheader(); writer.writerow(result)
    report("Lambda control evidence: " + json.dumps(result, sort_keys=True))
    return result
