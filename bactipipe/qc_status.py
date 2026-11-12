"""Shared interpretation; execution failure, measured QC and identity are distinct."""
from .genome_sizes import GenomeSizeCatalog, normalize_organism_name
from .organism_names import expected_organism_name
from functools import lru_cache
import math
from .lambda_control import is_lambda_control


@lru_cache(maxsize=1)
def catalog():
    return GenomeSizeCatalog.from_package()


def organism_input(entered):
    resolution = catalog().resolve(entered)
    warning = resolution.status in {"missing", "unsupported"}
    message = ""
    if warning:
        message = (f"Entered organism '{resolution.entered_name}' was not recognized; analyzed as Unknown."
                   if resolution.entered_name else "No organism was supplied; analyzed as Unknown.")
    return {"status": resolution.status, "warning": warning, "message": message,
            "effective_organism": resolution.canonical_name or "Unknown"}


def identification_state(expected, identified):
    hit = normalize_organism_name(identified)
    if hit.startswith("analysis failed"):
        return "analysis_failed"
    if (not hit or hit in {"unknown", "n/a", "na", "none", "not identified", "unidentified", "no hit"}
            or hit.startswith(("closest:", "ambiguous:", "no confident", "mixed bacterial", "possible mixed", "no hit", "secondary bacterial", "suspected culture"))):
        return "unresolved"
    resolution = catalog().resolve(expected)
    if resolution.is_unknown:
        return "identified"
    return "matched" if expected_organism_name(hit) == expected_organism_name(resolution.canonical_name) else "mismatch"


def quality_status(verdicts):
    values = [str(v).strip().casefold() for v in verdicts]
    if "fail" in values:
        return "Fail"
    if any(v != "pass" for v in values):
        return "Not assessed"
    return "Pass"


def threshold_verdict(value, threshold, *, maximum=False):
    """A missing/non-finite measurement is not a measured threshold failure."""
    try:
        value, threshold = float(value), float(threshold)
    except (TypeError, ValueError):
        return "N/A"
    if not math.isfinite(value) or not math.isfinite(threshold):
        return "N/A"
    return "Pass" if (value <= threshold if maximum else value >= threshold) else "Fail"


def interpret(expected, identified, verdicts, *, analysis_failed=False, isolate_qc_failure=False, isolate_qc_review=False, control_qc="Not assessed"):
    input_info = organism_input(expected)
    identification = identification_state(expected, identified)
    if is_lambda_control(expected):
        # Only quality/depth apply here; bacterial marker/identity verdicts do not.
        quality = quality_status(list(verdicts)[:2])
        failed = analysis_failed or identification == "analysis_failed" or quality == "Fail" or control_qc in {"Fail", "Error"}
        confirmed = control_qc == "Pass" and str(identified).strip().casefold() == "lambda"
        return {"overall": "Fail" if failed else "Pass" if confirmed and quality == "Pass" else "Review", "quality": quality,
                "identification": "analysis_failed" if analysis_failed or identification == "analysis_failed"
                else "matched" if confirmed else "control_not_confirmed" if control_qc == "Fail"
                else "control_not_assessed", "input": input_info}
    quality = quality_status(verdicts)
    if analysis_failed or isolate_qc_failure or identification == "analysis_failed" or quality == "Fail":
        overall = "Fail"
    elif quality != "Pass" or identification == "unresolved" or isolate_qc_review:
        overall = "Review"
    else:
        overall = "Pass"
    return {"overall": overall, "quality": quality, "identification": identification, "input": input_info}
