from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import re
import json
import math
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple, Union
from xml.sax.saxutils import escape

import pandas as pd

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle


# ──────────────────────────────────────────────────────────────────────────────
# Layout configuration
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class LayoutConfig:
    # Page & margins
    page_size: Tuple[float, float] = landscape(letter)
    left_margin: float = 36
    right_margin: float = 36
    top_margin: float = 48
    bottom_margin: float = 36

    # Dual header band split (portion of usable width for left vs right panels)
    left_panel_ratio: float = 0.36
    right_panel_ratio: float = 0.64

    # Column width strategy for main QC table
    reserved_cols: Tuple[int, int] = (2, 3)     # Expected & Identified organism (0-based)
    reserved_width: float = 2.0 * inch
    last3_narrow_width: float = 0.9 * inch

    # Typography (separate sizes for header blocks vs the table)
    header_block_font_size: int = 10            # used for Run Metadata + Tools blocks
    table_header_font_size: int = 8             # table header row
    table_body_font_size: int = 8               # table data rows

    # Header wrapping threshold for table columns
    wrap_threshold: int = 12

    # Body alignment: columns to LEFT-align (0-based). Others are centered.
    left_align_cols: Tuple[int, ...] = (0, 2, 3)

    # Conditional coloring
    pass_color: colors.Color = colors.HexColor("#1b5e20")
    fail_color: colors.Color = colors.HexColor("#b71c1c")


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

def render_pdf_from_rows(
    rows: Sequence[Sequence[object]],
    header_labels: Sequence[str],
    pdf_path: str,
    run_metadata: Optional[Iterable[str]] = None,
    legend: Optional[Iterable[str]] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    tools_header: Optional[Sequence[Tuple[str, str]]] = None,  # list of (tool, version/cutoff)
    layout: LayoutConfig = LayoutConfig(),
    run_name: Optional[str] = None,  # NEW: explicit run name override
    lambda_evidence: Optional[dict] = None,
) -> str:
    """Render a QC report PDF from in-memory rows."""
    df = pd.DataFrame(rows, columns=list(header_labels))
    # If no explicit run_name, try to infer from run_metadata
    if not run_name and run_metadata:
        run_name = _infer_run_name_from_metadata(list(run_metadata))
    return _render_pdf_common(
        df=df,
        pdf_path=pdf_path,
        run_metadata=list(run_metadata) if run_metadata else None,
        legend=list(legend) if legend else None,
        title=title,
        subtitle=subtitle,
        tools_header=list(tools_header) if tools_header else None,
        layout=layout,
        run_name=run_name,
        lambda_evidence=lambda_evidence,
    )


def render_pdf_from_tsv(
    tsv_path: str,
    pdf_path: str,
    run_metadata: Optional[Iterable[str]] = None,
    legend: Optional[Iterable[str]] = None,
    title: Optional[str] = None,
    subtitle: Optional[str] = None,
    layout: LayoutConfig = LayoutConfig(),
    run_name: Optional[str] = None,  # NEW: explicit run name override
    lambda_evidence: Optional[dict] = None,
) -> str:
    """Render a QC report PDF directly from an already written TSV file."""
    lines = _read_lines(tsv_path)
    header_idx = _find_table_header_idx(lines)
    if header_idx is None:
        raise RuntimeError("Could not locate the table header line containing 'Sample ID'.")

    # Parse region above the table: left metadata lines (with any '====' preserved)
    # and a tools span (tool >> value) surrounded by separators if present.
    top_block, tools_pairs = _parse_header_blocks(lines, header_idx)

    table_text = "\n".join(lines[header_idx:])
    df = pd.read_csv(StringIO(table_text), sep="\t", keep_default_na=False)
    if lambda_evidence is None:
        lambda_evidence = _load_lambda_evidence(df, Path(tsv_path).parent)

    run_meta = list(run_metadata) if run_metadata else (top_block or None)

    # Infer run name if not provided
    if not run_name:
        run_name = _infer_run_name_from_metadata(run_meta)

    return _render_pdf_common(
        df=df,
        pdf_path=pdf_path,
        run_metadata=run_meta,
        legend=list(legend) if legend else None,
        title=title,
        subtitle=subtitle,
        tools_header=tools_pairs or None,
        layout=layout,
        run_name=run_name,
        lambda_evidence=lambda_evidence,
    )


# ──────────────────────────────────────────────────────────────────────────────
# Parsing helpers
# ──────────────────────────────────────────────────────────────────────────────

def _read_lines(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().splitlines()


def _find_table_header_idx(lines: List[str]) -> Optional[int]:
    for i, line in enumerate(lines):
        if "Sample ID" in line and "\t" in line:
            return i
    return None


def _is_separator(line: str) -> bool:
    # Be permissive: allow lines that are mostly '=' after stripping whitespace.
    s = line.strip()
    if len(s) < 3:
        return False
    eq_count = sum(1 for ch in s if ch == '=')
    return eq_count >= max(3, int(0.8 * len(s)))


def _parse_header_blocks(lines: List[str], header_idx: int) -> Tuple[List[str], List[Tuple[str, str]]]:
    """
    Parse pre-table region (lines[:header_idx]) into:
      - top_block: metadata lines (preserve as-is, including any '====' lines)
      - tools_pairs: list of (tool, value) parsed from lines containing '>>'
                     that are inside the last [==== ... ====] span in the region.
    If no such tools span is found, everything is metadata and tools_pairs is empty.
    """
    region = lines[:header_idx]

    # Identify all separator positions
    sep_positions = [i for i, ln in enumerate(region) if _is_separator(ln)]

    # Collect candidate spans [start, end] bounded by separators
    spans: List[Tuple[int, int]] = []
    for i in range(len(sep_positions) - 1):
        start = sep_positions[i]
        end = sep_positions[i + 1]
        if end > start + 1:
            spans.append((start, end))

    # Choose the last span that contains at least one '>>' line as the tools block
    tools_span: Optional[Tuple[int, int]] = None
    for start, end in spans:
        if any(">>" in ln for ln in region[start + 1:end]):
            tools_span = (start, end)

    if tools_span:
        start, end = tools_span
        # Metadata: everything BEFORE the tools span (preserve including any '====')
        top_block = [ln for ln in region[:start] if ln.strip()]
        # Tools: parse '>>' lines only
        tools_pairs = []
        for ln in region[start + 1:end]:
            if ">>" in ln:
                k, v = ln.split(">>", 1)
                tools_pairs.append((k.strip(), v.strip()))
        return top_block, tools_pairs

    # No tools span found → treat all as metadata
    top_block = [ln for ln in region if ln.strip()]
    return top_block, []


def _infer_run_name_from_metadata(run_meta: Optional[List[str]]) -> Optional[str]:
    if not run_meta:
        return None
    for ln in run_meta:
        if ln.lower().startswith("run:"):
            return ln.split(":", 1)[1].strip()
    return None


# ──────────────────────────────────────────────────────────────────────────────
# Styles & layout helpers
# ──────────────────────────────────────────────────────────────────────────────

def _styles(layout: LayoutConfig):
    from pathlib import Path
    import reportlab
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont

    heading_font = "BactiPipeVeraBold"
    if heading_font not in pdfmetrics.getRegisteredFontNames():
        font_path = Path(reportlab.__file__).resolve().parent / "fonts" / "VeraBd.ttf"
        pdfmetrics.registerFont(TTFont(heading_font, str(font_path)))
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="TitleCenter", parent=styles["Title"], alignment=TA_CENTER))
    styles.add(ParagraphStyle(
        name="HeaderBlock",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=layout.header_block_font_size,
        leading=layout.header_block_font_size + 2,
    ))
    styles.add(ParagraphStyle(name="Header", parent=styles["Heading2"], spaceAfter=4))
    styles.add(ParagraphStyle(name="HeaderRight", parent=styles["Heading2"],
                              fontName=heading_font, spaceAfter=4))
    styles.add(ParagraphStyle(
        name="ToolLabel", parent=styles["HeaderBlock"], fontSize=9, leading=11,
        splitLongWords=False,
    ))
    styles.add(ParagraphStyle(
        name="ToolValue", parent=styles["HeaderBlock"], fontSize=9, leading=11,
    ))
    styles.add(ParagraphStyle(
        name="DatasetValue", parent=styles["HeaderBlock"], fontSize=8.5, leading=10.5,
        spaceAfter=3,
    ))
    styles.add(ParagraphStyle(
        name="HeaderSection", parent=styles["HeaderBlock"], fontName="Helvetica-Bold",
        fontSize=10, leading=12, spaceBefore=5, spaceAfter=3,
    ))
    styles.add(ParagraphStyle(
        name="Separator",
        parent=styles["Normal"],
        fontName="Courier",
        fontSize=layout.header_block_font_size,
        leading=layout.header_block_font_size + 2,
    ))
    # NEW: style for wrapped organism names in table body
    styles.add(ParagraphStyle(
        name="TableBodyPara",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=layout.table_body_font_size,
        leading=layout.table_body_font_size + 2,
    ))
    return styles


def _page_footer_factory(run_name: Optional[str]):
    """Returns a ReportLab onPage callback that draws the desired footer."""
    prefix = "BactiPipe QC Report"
    if run_name:
        label = f"{prefix} - Run {run_name}"
    else:
        label = prefix

    def _footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 8)
        text = f"{label}    |    Page {doc.page}"
        # Right-align to the page's right margin
        canvas.drawRightString(doc.pagesize[0] - 36, 20, text)
        canvas.restoreState()

    return _footer


def _wrap_header(text: str, threshold: int = 12) -> str:
    if len(text) <= threshold:
        return text
    parts = text.split()
    if len(parts) > 1:
        mid = len(parts) // 2
        return "\n".join([" ".join(parts[:mid]), " ".join(parts[mid:])])
    mid = len(text) // 2
    return text[:mid] + "\n" + text[mid:]


def _compute_col_widths(df: pd.DataFrame, layout: LayoutConfig) -> List[float]:
    page_w, _ = layout.page_size
    usable = page_w - (layout.left_margin + layout.right_margin)

    n_cols = df.shape[1]
    last3_total = 3 * layout.last3_narrow_width
    reserved_total = 2 * layout.reserved_width
    remaining = max(usable - (last3_total + reserved_total), 1.0)

    other_cols = n_cols - 5  # all except last 3 and reserved 2
    other_width = remaining / other_cols if other_cols > 0 else remaining

    widths: List[float] = []
    for i in range(n_cols):
        if i in layout.reserved_cols:
            widths.append(layout.reserved_width)
        elif i >= n_cols - 3:
            widths.append(layout.last3_narrow_width)
        else:
            widths.append(other_width)
    return widths


# ──────────────────────────────────────────────────────────────────────────────
# Building sections
# ──────────────────────────────────────────────────────────────────────────────

def _report_header_groups(tools_header):
    """Presentation only: preserve original identifiers/filter settings in TSVs."""
    software, datasets, thresholds = [], [], []
    seen = set()
    has_source_db = any(str(k).strip().casefold() == "speciesfinder db" for k, v in tools_header or [])
    threshold_labels = {
        "required mean quality": "Minimum mean quality",
        "completeness cutoff": "Minimum completeness",
        "contamination cutoff": "Maximum contamination",
        "read filtering trigger": "Read filtering trigger",
    }
    for name, value in tools_header or []:
        name, value = str(name).strip(), str(value).strip()
        key = " ".join(name.casefold().split())
        if key in {"minimum read length", "min read length", "minimum required length"}:
            continue
        if key == "speciesfinder db":
            pair = ("SpeciesFinder db", value)
            target = datasets
        elif key in {"species identification database", "species id references", "speciesfinder database"} and (
            value.startswith("BactiPipe-SpeciesFinder-") or key == "speciesfinder database"
        ):
            if has_source_db:
                continue
            pair = ("SpeciesFinder db", value.removeprefix("BactiPipe-SpeciesFinder-"))
            target = datasets
        elif "database" in key or key in {"reference bundle", "checkm references", "species id references"}:
            continue  # Other database identities stay in logs/TSVs.
        elif key in threshold_labels or "cutoff" in key or key.startswith(("minimum ", "maximum ", "required ")):
            pair = (threshold_labels.get(key, name), value)
            target = thresholds
        else:
            # Normalize only software versions, never database release/build identifiers.
            value = re.sub(r"^" + re.escape(name) + r"\s+", "", value, flags=re.I)
            value = re.sub(r"^version\s*:?\s*", "", value, flags=re.I)
            if re.fullmatch(r"v?\d+(?:\.\d+)+(?:[+-][\w.-]+)?", value, flags=re.I):
                value = "v" + value.lstrip("vV")
            if value.casefold() in {"v?", "?", "unknown"}:
                value = "Not recorded"
            pair = (name, value)
            target = software
        identity = (pair[0].casefold(), pair[1])
        if identity not in seen:
            target.append(pair)
            seen.add(identity)
    # Keep the depth policy immediately before quality, regardless of where
    # tool/control entries were recorded in the source header.
    depth = [pair for pair in thresholds if pair[0].casefold() == "minimum depth"]
    if depth:
        thresholds = [pair for pair in thresholds if pair not in depth]
        quality_index = next((i for i, pair in enumerate(thresholds)
                              if pair[0] == "Minimum mean quality"), 0)
        thresholds[quality_index:quality_index] = depth
    return software, datasets, thresholds


def _build_dual_header_band(
    run_metadata_lines: Optional[List[str]],
    tools_header: Optional[Sequence[Tuple[str, str]]],
    layout: LayoutConfig,
    styles,
) -> Optional[Table]:
    """
    Run metadata on the left; separate tools and threshold columns on the right.
    Restore the original double-dashed separators without narrowing tool columns.
    """
    if not run_metadata_lines and not tools_header:
        return None

    page_w, _ = layout.page_size
    usable = page_w - (layout.left_margin + layout.right_margin)
    left_ratio = max(0.05, min(0.95, layout.left_panel_ratio))
    right_ratio = max(0.05, min(0.95, layout.right_panel_ratio))
    total = left_ratio + right_ratio
    right_width = usable * right_ratio / total
    col_widths = [usable * left_ratio / total, right_width / 2, right_width / 2]

    def text_width(text, style_name):
        from reportlab.pdfbase.pdfmetrics import stringWidth
        style = styles[style_name]
        return stringWidth(str(text), style.fontName, style.fontSize)

    def separator(width):
        from reportlab.pdfbase.pdfmetrics import stringWidth
        style = styles["Separator"]
        count = max(3, int((width - 1) / stringWidth("=", style.fontName, style.fontSize)))
        return Paragraph("=" * count, style)

    # LEFT: Run Metadata. Long ASCII separators must not force extra wrapped lines.
    left_flow: List[Union[Paragraph, Table]] = []
    if run_metadata_lines:
        metadata = [str(ln) for ln in run_metadata_lines if str(ln).strip() and not _is_separator(ln)]
        metadata_width = min(col_widths[0] - 12, max(
            [text_width(ln, "HeaderBlock") for ln in metadata] or [0]
        ))
        left_flow.append(Paragraph("Run Metadata", styles["HeaderRight"]))
        left_flow.append(separator(metadata_width))
        for ln in metadata:
            left_flow.append(Paragraph(escape(ln), styles["HeaderBlock"]))
        left_flow.append(separator(metadata_width))

    # RIGHT: tools and thresholds occupy independent, top-aligned columns.
    software, datasets, thresholds = _report_header_groups(tools_header)

    def pair_panel(title, pairs, available_width):
        panel: List[Union[Paragraph, Table]] = []
        if not pairs:
            return panel
        label_width = min(available_width * .80, max(text_width(k, "ToolLabel") for k, v in pairs) + 12)
        value_width = min(available_width - label_width, max(text_width(v, "ToolValue") for k, v in pairs) + 1)
        section_width = label_width + value_width
        panel.append(Paragraph(title, styles["HeaderRight"]))
        panel.append(separator(section_width))
        tool_data = [
            [Paragraph(escape(k), styles["ToolLabel"]),
             Paragraph(escape(v), styles["ToolValue"])]
            for k, v in pairs
        ]
        wide_policy_rows = []
        for index, (label, value) in enumerate(pairs):
            if label.casefold() == "minimum depth" and ";" in value:
                # A multi-tier policy needs the panel width, not a narrow
                # version-number cell. Keep all other label/value rows intact.
                tool_data[index] = [Paragraph(escape(label) + ": " + escape(value),
                                             styles["ToolValue"]), ""]
                wide_policy_rows.append(index)
        tool_table = Table(tool_data, colWidths=[label_width, value_width],
                           rowHeights=[styles["ToolValue"].leading if not k and not v else None for k, v in pairs],
                           hAlign="LEFT")
        tool_table.setStyle(TableStyle([
            ("ALIGN", (0,0), (-1,-1), "LEFT"),
            ("VALIGN", (0,0), (-1,-1), "TOP"),
            ("LEFTPADDING", (0,0), (-1,-1), 0),
            ("RIGHTPADDING", (0,0), (-1,-1), 0),
            ("TOPPADDING", (0,0), (-1,-1), 0),
            ("BOTTOMPADDING", (0,0), (-1,-1), 0),  # compact vertically
        ]))
        for index in wide_policy_rows:
            tool_table.setStyle(TableStyle([("SPAN", (0, index), (1, index))]))
        panel.append(tool_table)
        panel.append(separator(section_width))
        return panel

    versions = software + datasets
    tools_flow = pair_panel("Tool Versions", versions, col_widths[1] - 12)
    thresholds_flow = pair_panel("QC Thresholds", thresholds, col_widths[2])

    parent = Table([[left_flow, tools_flow, thresholds_flow]], colWidths=col_widths)
    parent.setStyle(TableStyle([
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("ALIGN", (0,0), (-1,-1), "LEFT"),
        ("LEFTPADDING", (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (-1,-1), 0),
        ("RIGHTPADDING", (0,0), (1,-1), 12),
        ("TOPPADDING", (0,0), (-1,-1), 0),
        ("BOTTOMPADDING", (0,0), (-1,-1), 0),
    ]))
    return parent


def _load_lambda_evidence(df, result_root):
    """Read only retained per-control metrics; never infer coverage from identity."""
    from bactipipe.lambda_control import is_lambda_control
    evidence = {}
    for _, row in df.iterrows():
        sample = str(row.get("Sample ID", ""))
        if not is_lambda_control(row.get("Expected organism")):
            continue
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sample) or sample in {".", ".."}:
            continue
        path = result_root / "qc_out" / sample / "lambda_control" / "result.json"
        try:
            payload = json.loads(path.read_text())
            if isinstance(payload, dict):
                evidence[sample] = payload
        except (OSError, ValueError):
            # Historical standalone TSVs may not have their companion metrics.
            # Report unavailable coverage instead of inventing a measurement.
            continue
    return evidence


def _lambda_match_cell(identity, evidence, styles, layout):
    def numeric(value):
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return number if math.isfinite(number) and 0 <= number <= 100 else None

    identity = numeric(evidence.get("identity_percent", identity))
    coverage = numeric(evidence.get("reference_coverage_percent"))

    def metric(value, threshold):
        text = f"{value:,.2f}%" if value is not None else "N/A"
        minimum = numeric(threshold)
        if value is not None and minimum is not None and value < minimum:
            return f'<font color="{layout.fail_color.hexval()}">{text}</font>'
        return text

    id_text = metric(identity, evidence.get("minimum_identity_percent"))
    cov_text = metric(coverage, evidence.get("minimum_coverage_percent"))
    style = ParagraphStyle("LambdaMatch", parent=styles["TableBodyPara"], alignment=TA_CENTER)
    return Paragraph(id_text + '<br/><font size="6">Coverage: ' + cov_text + '</font>', style)


def _build_table(df: pd.DataFrame, layout: LayoutConfig, lambda_evidence=None) -> Table:
    # Machine-readable isolate flags travel in TSVs, not as extra report columns.
    df = df.drop(columns=["Isolate QC", "Control QC"], errors="ignore")
    styles = _styles(layout)
    from xml.sax.saxutils import escape
    from bactipipe.qc_status import organism_input
    from bactipipe.lambda_control import is_lambda_control

    if "Expected organism" in df.columns:
        df = df.iloc[sorted(range(len(df)), key=lambda i: not is_lambda_control(df.iloc[i]["Expected organism"]))]

    # Headers: wrap except reserved organism columns
    wrapped_headers = []
    for i, c in enumerate(df.columns):
        txt = str(c)
        if i in layout.reserved_cols:
            wrapped_headers.append(txt)  # no wrapping for these headers
        else:
            wrapped_headers.append(_wrap_header(txt, layout.wrap_threshold))

    # Body rows: wrap only organism columns (reserved_cols)
    data = [wrapped_headers]
    org_cols = set(layout.reserved_cols)  # usually {2, 3}
    for row in df.astype(str).values.tolist():
        control = "Expected organism" in df.columns and is_lambda_control(row[df.columns.get_loc("Expected organism")])
        sample = row[df.columns.get_loc("Sample ID")] if "Sample ID" in df.columns else ""
        r = []
        for i, cell in enumerate(row):
            if control and df.columns[i] == "Expected organism":
                cell = "Lambda (Control)"
            if control and df.columns[i] in {"CheckM completeness", "CheckM contamination"}:
                cell = "N/A"
            if df.columns[i] == "Expected organism" and organism_input(cell)["warning"]:
                if not control:
                    cell = f"{cell or 'Not supplied'} (! analyzed as Unknown)"
            if df.columns[i] in {"Mean quality", "Coverage Depth", "CheckM completeness", "CheckM contamination"} and cell.strip().upper() in {"N/A", "NA", "NAN", ""}:
                cell = "N/A" if control and df.columns[i].startswith("CheckM") else "Not assessed"
            if control and df.columns[i] == "% Match":
                r.append(_lambda_match_cell(cell, (lambda_evidence or {}).get(sample, {}), styles, layout))
            elif i in org_cols:
                r.append(Paragraph(escape(cell), styles["TableBodyPara"]))
            else:
                r.append(cell)
        data.append(r)

    # Column widths
    col_widths = _compute_col_widths(df, layout)
    t = Table(data, repeatRows=1, colWidths=col_widths)

    style = TableStyle([
        # Header row
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#f0f0f0")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.black),
        ("VALIGN", (0,0), (-1,0), "MIDDLE"),
        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
        ("FONTSIZE", (0,0), (-1,0), layout.table_header_font_size),
        ("BOTTOMPADDING", (0,0), (-1,0), 6),

        # Body defaults
        ("ALIGN", (0,1), (-1,-1), "CENTER"),
        ("VALIGN", (0,1), (-1,-1), "MIDDLE"),
        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
        ("FONTSIZE", (0,1), (-1,-1), layout.table_body_font_size),
    ])

    # Alternating row backgrounds
    for i in range(1, len(data)):
        style.add("BACKGROUND", (0,i), (-1,i),
                  colors.HexColor("#fcfcfc") if i % 2 == 1 else colors.white)

    # Keep the compact Lambda coverage annotation on one line without widening
    # the shared Match column or changing the other samples' cell padding.
    if "% Match" in df.columns and "Expected organism" in df.columns:
        match_idx = df.columns.get_loc("% Match")
        for row_idx, expected in enumerate(df["Expected organism"], start=1):
            if is_lambda_control(expected):
                style.add("LEFTPADDING", (match_idx, row_idx), (match_idx, row_idx), 3)
                style.add("RIGHTPADDING", (match_idx, row_idx), (match_idx, row_idx), 3)

    # Match header alignment to body
    n_cols = df.shape[1]
    style.add("ALIGN", (0,0), (-1,0), "CENTER")  # default header center
    for col in layout.left_align_cols:
        if 0 <= col < n_cols:
            style.add("ALIGN", (col, 0), (col, 0), "LEFT")   # header col
            style.add("ALIGN", (col, 1), (col, -1), "LEFT")  # body col

    # Conditional coloring for Overall Quality
    if "Overall Quality" in df.columns:
        oq_idx = df.columns.get_loc("Overall Quality")
        for r_idx, val in enumerate(df["Overall Quality"].astype(str), start=1):
            v = val.strip().lower()
            if v == "pass":
                style.add("TEXTCOLOR", (oq_idx, r_idx), (oq_idx, r_idx), layout.pass_color)
            elif v == "fail":
                style.add("TEXTCOLOR", (oq_idx, r_idx), (oq_idx, r_idx), layout.fail_color)
            elif v == "review":
                style.add("TEXTCOLOR", (oq_idx, r_idx), (oq_idx, r_idx), colors.HexColor("#946200"))

    t.setStyle(style)
    return t


# ──────────────────────────────────────────────────────────────────────────────
# Document assembly
# ──────────────────────────────────────────────────────────────────────────────

def _render_pdf_common(
    df: pd.DataFrame,
    pdf_path: str,
    run_metadata: Optional[List[str]],
    legend: Optional[List[str]],
    title: Optional[str],
    subtitle: Optional[str],
    tools_header: Optional[Sequence[Tuple[str, str]]],
    layout: LayoutConfig,
    run_name: Optional[str],  # used in footer
    lambda_evidence: Optional[dict] = None,
) -> str:
    styles = _styles(layout)

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=layout.page_size,
        rightMargin=layout.right_margin,
        leftMargin=layout.left_margin,
        topMargin=layout.top_margin,
        bottomMargin=layout.bottom_margin,
    )

    story: List[object] = []

    # Title / subtitle
    if title:
        story.append(Paragraph(title, styles["TitleCenter"]))
    if subtitle:
        story.append(Paragraph(subtitle, styles["HeaderBlock"]))
    if title or subtitle:
        story.append(Spacer(1, 0.2 * inch))

    # Dual header band
    band = _build_dual_header_band(run_metadata, tools_header, layout, styles)
    if band:
        story.append(band)
        story.append(Spacer(1, 0.15 * inch))

    # Legend (optional)
    if legend:
        story.append(Paragraph("Legend", styles["Header"]))
        for l in legend:
            story.append(Paragraph("• " + str(l), styles["HeaderBlock"]))
        story.append(Spacer(1, 0.2 * inch))

    # Main QC table
    story.append(Paragraph("Quality Summary Table", styles["Header"]))
    story.append(_build_table(df, layout, lambda_evidence=lambda_evidence))

    # Footer with run name + page number
    footer_func = _page_footer_factory(run_name)

    doc.build(story, onFirstPage=footer_func, onLaterPages=footer_func)
    return pdf_path
