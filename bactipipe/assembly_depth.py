"""Post-assembly yield depth, independent of organism identification."""
import hashlib
import json
import math
from pathlib import Path


def assembly_depth(total_bases, assembly, minimum_depth, *, audit_path=None):
    """Return depth/verdict using all assembled contigs, including ambiguous bases.

    This is yield divided by assembly length, not alignment-measured mean depth.
    Call with QC-retained bases, never an assembly-only subset's base count.
    """
    length = 0
    contigs = 0
    digest = hashlib.sha256()
    reason = ""
    try:
        with Path(assembly).open("rb") as handle:
            for raw in handle:
                digest.update(raw)
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(b">"):
                    contigs += 1
                elif contigs:
                    length += len(b"".join(line.split()))
                else:
                    raise ValueError("Assembly is not FASTA")
        if not length or not contigs:
            raise ValueError("Assembly contains no sequence")
    except (OSError, ValueError) as error:
        reason = str(error)
        length = 0
    try:
        bases = float(total_bases)
        if not math.isfinite(bases) or bases < 0:
            raise ValueError("Invalid retained base count")
    except (TypeError, ValueError):
        bases = None
        reason = reason or "QC-retained base count unavailable"
    depth = bases / length if length and bases is not None else None
    try:
        threshold = float(minimum_depth)
        if not math.isfinite(threshold) or threshold < 0:
            raise ValueError("Invalid depth threshold")
    except (TypeError, ValueError):
        threshold = None
    verdict = "N/A" if depth is None or threshold is None else ("Pass" if depth >= threshold else "Fail")
    if audit_path is not None:
        path = Path(audit_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": 1, "method": "qc_retained_bases / assembly_length_bp",
            "retained_bases": bases, "assembly_length_bp": length or None,
            "assembly_contigs": contigs, "assembly_sha256": digest.hexdigest() if length else None,
            "depth": depth, "minimum_depth": threshold, "verdict": verdict, "reason": reason,
        }
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(payload, indent=2) + "\n")
        temporary.replace(path)
    return depth, verdict
