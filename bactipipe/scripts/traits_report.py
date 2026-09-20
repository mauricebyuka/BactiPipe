# bactipipe/scripts/traits_report.py
from __future__ import annotations

import os
import csv
import re
import json
from typing import Dict, List, Tuple, Optional, Set, DefaultDict
from collections import defaultdict
from importlib.resources import files as resource_files
from xml.sax.saxutils import escape as _xml_escape
from datetime import datetime
import socket
# ReportLab
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, TableStyle, Spacer, KeepTogether
from reportlab.platypus.tables import Table as _RLTable, TableStyle as _RLTableStyle

from reportlab.platypus.flowables import Flowable
from reportlab.platypus import Paragraph
from math import radians, sin, cos
from xml.sax.saxutils import escape as _xml_escape

# BactiPipe utilities
from bactipipe.scripts.utils import make_report_header_block


# =============================================================================
# Safe Table wrapper to avoid ReportLab crashes on None rowHeights
# =============================================================================
class SafeTable(_RLTable):
    """A drop-in replacement for ReportLab Table that avoids identity/_culprit crashes."""
    def identity(self, maxLen=None):
        try:
            return super().identity(maxLen)
        except Exception:
            rows = len(getattr(self, "_cellvalues", []) or [])
            cols = len((getattr(self, "_cellvalues", [])[0]) or []) if rows else 0
            return f"SafeTable({rows}x{cols})"

    def _culprit(self, ident=None):
        try:
            return super()._culprit(ident)
        except Exception:
            rows = len(getattr(self, "_cellvalues", []) or [])
            cols = len((getattr(self, "_cellvalues", [])[0]) or []) if rows else 0
            return f"SafeTable({rows}x{cols})"


# =============================================================================
# Public TSV writers (used by run_traits.py)
# =============================================================================

def write_sample_tsvs(sample: str, merged: Dict[str, List[Dict[str, str]]], out_root: str) -> None:
    """
    Writes:
      <out_root>/<sample>.amr.tsv
      <out_root>/<sample>.vf.tsv
      <out_root>/<sample>.summary.tsv
    """
    os.makedirs(out_root, exist_ok=True)

    def _write(path: str, rows: List[Dict[str, str]], schema: List[str]):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=schema, delimiter="\t")
            w.writeheader()
            for r in rows or []:
                w.writerow({k: ("" if r.get(k) is None else r.get(k, "")) for k in schema})

    amr_schema = [
        "sample","determinant","type","phenotype","class","identity","coverage","length",
        "contig","start","end","strand","source","database","accession","tool","subtype"
    ]
    vf_schema = [
        "sample","gene","base_gene","category","note","identity","coverage","length",
        "contig","start","end","strand","source","database","accession","tool"
    ]

    _write(os.path.join(out_root, f"{sample}.amr.tsv"), merged.get("amr", []), amr_schema)
    _write(os.path.join(out_root, f"{sample}.vf.tsv"),  merged.get("vf",  []), vf_schema)

    summary_schema = ["sample","kind","name","phenotype_or_category","identity","coverage","source"]
    summary_rows: List[Dict[str, str]] = []
    for r in merged.get("amr", []) or []:
        summary_rows.append({
            "sample": r.get("sample",""),
            "kind": "AMR",
            "name": r.get("determinant",""),
            "phenotype_or_category": (r.get("phenotype") or r.get("class") or ""),
            "identity": r.get("identity",""),
            "coverage": r.get("coverage",""),
            "source": r.get("source",""),
        })
    for r in merged.get("vf", []) or []:
        summary_rows.append({
            "sample": r.get("sample",""),
            "kind": "Virulence",
            "name": r.get("gene",""),
            "phenotype_or_category": r.get("category",""),
            "identity": r.get("identity",""),
            "coverage": r.get("coverage",""),
            "source": r.get("source",""),
        })
    _write(os.path.join(out_root, f"{sample}.summary.tsv"), summary_rows, summary_schema)


def write_run_summary(all_merged: Dict[str, Dict[str, List[Dict[str, str]]]], out_root: str) -> str:
    os.makedirs(out_root, exist_ok=True)
    path = os.path.join(out_root, "run.summary.tsv")
    schema = ["sample","kind","name","phenotype_or_category","identity","coverage","source"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=schema, delimiter="\t")
        w.writeheader()
        for sample, merged in all_merged.items():
            for r in merged.get("amr", []) or []:
                w.writerow({
                    "sample": sample,
                    "kind": "AMR",
                    "name": r.get("determinant",""),
                    "phenotype_or_category": (r.get("phenotype") or r.get("class") or ""),
                    "identity": r.get("identity",""),
                    "coverage": r.get("coverage",""),
                    "source": r.get("source",""),
                })
            for r in merged.get("vf", []) or []:
                w.writerow({
                    "sample": sample,
                    "kind": "Virulence",
                    "name": r.get("gene",""),
                    "phenotype_or_category": r.get("category",""),
                    "identity": r.get("identity",""),
                    "coverage": r.get("coverage",""),
                    "source": r.get("source",""),
                })
    return path


# =============================================================================
# Catalog loaders (from bactipipe/data)
# =============================================================================
def _load_vf_catalog() -> List[Dict[str, str]]:
    """
    vf_category_map.tsv  columns:
      gene    functional_category   notes
    """
    path = resource_files("bactipipe.data").joinpath("vf_category_map.tsv")
    rows: List[Dict[str, str]] = []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for rec in r:
            g = (rec.get("gene") or "").strip()
            if not g:
                continue
            rows.append({
                "gene": g,
                "category": (rec.get("functional_category") or "Misc").strip() or "Misc",
                "note": (rec.get("notes") or "").strip(),
            })
    rows.sort(key=lambda x: (x["category"].lower(), x["gene"].lower()))
    return rows


def _load_amr_maps() -> Tuple[List[Dict[str,str]], Dict[str,str]]:
    """
    amr_acquired_map.tsv  columns:
      gene    phenotype

    amr_mutation_map.tsv columns:
      mutation(or gene)  phenotype
      (we store phenotype under base gene key, e.g., 'gyrA')
    """
    acq_path = resource_files("bactipipe.data").joinpath("amr_acquired_map.tsv")
    mut_path = resource_files("bactipipe.data").joinpath("amr_mutation_map.tsv")

    acquired: List[Dict[str,str]] = []
    mut_map: Dict[str,str] = {}

    if os.path.exists(acq_path):
        with open(acq_path) as f:
            r = csv.DictReader(f, delimiter="\t")
            for rec in r:
                g = (rec.get("gene") or "").strip()
                if g:
                    acquired.append({"gene": g, "phenotype": (rec.get("phenotype") or "").strip()})
    acquired.sort(key=lambda x: x["gene"].lower())

    if os.path.exists(mut_path):
        with open(mut_path) as f:
            r = csv.DictReader(f, delimiter="\t")
            for rec in r:
                g = (rec.get("mutation") or rec.get("gene") or "").strip()
                if not g:
                    continue
                base = _base_token(g)
                ph  = (rec.get("phenotype") or "").strip()
                if base and ph:
                    mut_map[base] = ph

    return acquired, mut_map


# =============================================================================
# Presence helpers
# =============================================================================
def _base_token(s: str) -> str:
    if not s:
        return ""
    s = s.strip()
    for sep in ("_", "-", "/", " "):
        if sep in s:
            return s.split(sep, 1)[0]
    return s


def _samples_in_run(all_merged: Dict[str, Dict[str, List[Dict[str, str]]]]) -> List[str]:
    return list(all_merged.keys())


def _vf_presence(all_merged: Dict[str, Dict[str, List[Dict[str, str]]]]) -> Dict[str, Set[str]]:
    present: Dict[str, Set[str]] = {}
    for sample, merged in all_merged.items():
        for r in merged.get("vf", []) or []:
            g = _base_token(r.get("gene", ""))
            if not g:
                continue
            present.setdefault(g, set()).add(sample)
    return present

# ===== Helpers to format tool/db info for PDF header =====

def _fmt_tools_for_header(tools: dict) -> list[str]:
    """Return display lines like 'AMRFinderPlus: v4.0.23' with a stable order."""
    if not tools:
        return []
    order = ["AMRFinderPlus", "ABRicate", "VirulenceFinder"]
    lines = []
    for k in order:
        if k in tools:
            lines.append(f"{k}: {tools[k]}")
    for k in sorted(set(tools.keys()) - set(order)):
        lines.append(f"{k}: {tools[k]}")
    return lines

def _fmt_dbs_for_header(dbs: dict) -> list[str]:
    if not dbs:
        return []
    names = {
        "amrfinder_db": "AMRFinder DB",
        "resfinder_db": "ResFinder DB",
        "virulencefinder_db": "VirulenceFinder DB",
        "abricate_resfinder": "ABRicate ResFinder DB",
        "abricate_card": "ABRicate CARD DB",
        "abricate_vfdb": "ABRicate VFDB",
    }
    order = [
        "amrfinder_db",
        "resfinder_db",
        "virulencefinder_db",
        "abricate_resfinder",
        "abricate_card",
        "abricate_vfdb",
    ]
    lines = []
    for k in order:
        if k in dbs:
            lines.append(f"{names.get(k, k)}: {dbs[k]}")
    for k in sorted(set(dbs.keys()) - set(order)):
        lines.append(f"{names.get(k, k)}: {dbs[k]}")
    return lines

# --- AMR collection directly from merged rows (record-driven) ---
_MUT_RX = re.compile(r"^[A-Za-z0-9]+_[A-Z]\d+[A-Z]$")  # e.g., gyrA_S83F

def _first_nonempty(rec: dict, *keys: str) -> str:
    for k in keys:
        v = rec.get(k)
        if v is not None:
            s = str(v).strip()
            if s:
                return s
    return ""

def _is_mutation_amr(rec: dict) -> bool:
    # AMRFinder: Type=AMR always; Subtype=POINT → mutation
    rtype   = _first_nonempty(rec, "type", "Type").lower()
    subtype = _first_nonempty(rec, "subtype", "Subtype").lower()
    det     = _first_nonempty(rec, "determinant", "Element symbol", "name", "gene", "mutation")
    if subtype in {"point", "snp", "snv"}:
        return True
    if rtype in {"mutation","mut","variant","substitution"}:
        return True
    return bool(_MUT_RX.match(det))


def _collect_amr_presence(
    all_merged: dict,
    acquired_catalog: List[Dict[str,str]],
    mut_pheno_map: Dict[str,str],
    *,
    out_root: Optional[str] = None,
    verbose_print: bool = False,
) -> Tuple[Dict[str, Set[str]], Dict[str, str], Dict[str, Set[str]], Dict[str, str]]:
    """
    Returns:
      acq_present:  dict[determinant]       -> set(samples)
      acq_pheno:    dict[determinant]       -> phenotype
      mut_present:  dict[mutation_string]   -> set(samples) (e.g., 'gyrA_S83F')
      mut_pheno:    dict[base_gene]         -> phenotype    (e.g., 'gyrA' -> ...)
    Also writes reports/debug/amr_rows_debug.tsv if out_root is provided.
    """
    acq_present: DefaultDict[str, Set[str]] = defaultdict(set)
    mut_present: DefaultDict[str, Set[str]] = defaultdict(set)
    acq_pheno: Dict[str,str] = {}
    mut_pheno: Dict[str,str] = dict(mut_pheno_map)

    catalog_ph = {d["gene"]: d.get("phenotype","") for d in (acquired_catalog or [])}
    dbg_rows: List[List[str]] = []

    for sample, merged in (all_merged or {}).items():
        rows = (merged.get("amr") or [])
        for r in rows:
            det = _first_nonempty(r, "determinant", "Element symbol", "name", "gene", "mutation")
            if not det:
                continue

            # Prefer explicit phenotype fields; fall back to Class/Subclass if needed
            ph = _first_nonempty(r, "phenotype", "Phenotype")
            if not ph:
                cls = _first_nonempty(r, "Class", "class")
                sub = _first_nonempty(r, "Subclass", "subclass")
                if cls and sub and sub.lower() != cls.lower():
                    ph = f"{cls} / {sub}"
                elif cls:
                    ph = cls

            if _is_mutation_amr(r):
                mut_present[det].add(sample)
                base = det.split("_", 1)[0] if "_" in det else det
                if base and base not in mut_pheno and ph:
                    mut_pheno[base] = ph
                dbg_rows.append([sample,"MUT",det,base,ph or "", _first_nonempty(r,"Subtype","subtype") or "-"])
            else:
                acq_present[det].add(sample)
                if det not in acq_pheno:
                    acq_pheno[det] = ph or catalog_ph.get(det, "")
                dbg_rows.append([sample,"ACQ",det,"-",acq_pheno.get(det,""), _first_nonempty(r,"Subtype","subtype") or "-"])

    if out_root:
        dbg_dir = os.path.join(out_root, "reports", "debug")
        os.makedirs(dbg_dir, exist_ok=True)
        with open(os.path.join(dbg_dir, "amr_rows_debug.tsv"), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["sample","kind","key","base_or_dash","phenotype","subtype"])
            w.writerows(dbg_rows)

        # also dump the merged structure for inspection
        with open(os.path.join(dbg_dir, "all_merged.debug.json"), "w") as f:
            json.dump(all_merged, f, indent=2)

    if verbose_print:
        print(f"[traits_report] AMR collected: acquired={len(acq_present)} determinants, mutations={len(mut_present)}")

    return acq_present, acq_pheno, mut_present, mut_pheno


# =============================================================================
# Safe table helpers
# =============================================================================
class RotatedLabel(Flowable):
    """
    Rotated label anchored at the bottom-left of the table cell.
    Uses a fixed header height (box_h) so the header row stays compact.
    """
    def __init__(self, text_or_para, style, angle: float = 55.0, pad: float = 2.0, box_h: float = 48.0):
        super().__init__()
        self.style = style
        self.angle = angle
        self.pad = pad
        self.box_h = box_h

        if isinstance(text_or_para, Paragraph):
            self._para_in = text_or_para
            self._text = None
        else:
            self._para_in = None
            self._text = str(text_or_para)

    def wrap(self, availWidth, availHeight):
        if self._para_in is not None:
            p = self._para_in
        else:
            p = Paragraph(_xml_escape(self._text) if self._text else "", self.style)

        self._para = p
        self._col_w = float(availWidth)
        # Fixed header height returned to the Table
        return self._col_w, float(self.box_h)

    def draw(self):
        c = self.canv
        c.saveState()
        # Bottom-left anchor inside the cell
        c.translate(self.pad, self.pad)      # <-- anchor at BL, not top
        c.rotate(self.angle)                 # rotate CCW
        # Draw the paragraph with its "lower-left" at (0, 0)
        w0, h0 = self._para.wrapOn(c, 1000, 1000)
        self._para.drawOn(c, 0, 0)
        c.restoreState()

def _normalize_table_data(data: List[List]) -> List[List]:
    """
    Keep Flowables (Paragraph, RotatedLabel, etc.) so they render correctly.
    Stringify only plain scalars; replace None with "".
    """
    out: List[List] = []
    for row in data:
        nr: List = []
        for v in row:
            if isinstance(v, (Paragraph, Flowable)):
                nr.append(v)          # keep as-is so Table can draw it
            elif v is None:
                nr.append("")
            else:
                nr.append(str(v))
        out.append(nr)
    return out

def _safe_table(
    data: list[list],
    *,
    repeat_header: bool = True,
    style=None,
    colWidths=None,
):
    data = _normalize_table_data(data)
    t = SafeTable(data, repeatRows=1 if repeat_header else 0, colWidths=colWidths)
    if isinstance(style, TableStyle):
        t.setStyle(style)
    elif isinstance(style, (list, tuple)):
        t.setStyle(TableStyle(list(style)))
    return t

# =============================================================================
# Run-level PDF generator
# =============================================================================
def _vf_rows_from_merged(all_merged: dict, sample_names: list[str]):
    """
    Build virulence presence/absence grouped by category directly from normalized rows in all_merged.

    Each VF row in all_merged[s]["vf"] must include:
      - gene        (display_name already)
      - category
      - note

    Returns:
      headers: ["Gene", "Note", *sample_names]  (we won't use headers text; we keep your existing "Gene/Description")
      grouped: dict(category -> list of [gene, note, ✓/"" per sample in sample_names order])
    """
    keys = []
    seen = set()
    present = defaultdict(set)  # (cat,gene,note) -> set(sample_name)

    for sname, merged in all_merged.items():
        for r in (merged.get("vf") or []):
            cat = (r.get("category") or "").strip()
            gene = (r.get("gene") or "").strip()
            note = (r.get("note") or "").strip()
            if not cat or not gene:
                continue
            k = (cat, gene, note)
            if k not in seen:
                seen.add(k)
                keys.append(k)
            present[k].add(sname)

    grouped = defaultdict(list)
    for (cat, gene, note) in keys:
        row = [gene, note]
        row.extend("✓" if s in present[(cat, gene, note)] else "" for s in sample_names)
        grouped[cat].append(row)

    headers = ["Gene", "Note", *sample_names]
    return headers, grouped

def render_run_pdf(
    *,
    out_root: str,
    run_name: str,
    accession: str,
    tech_name: str,
    versions: dict,
    thresholds: dict,
    all_merged: dict[str, dict[str, list[dict]]],
    pipeline_title: str = "Virulence and AMR Gene Detection",
    run_amr: bool = True,
    run_vf: bool = True,
):

    reports_dir = os.path.join(out_root, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    pdf_path = os.path.join(reports_dir, f"{accession}_vir_amr_detect.pdf")

    samples = _samples_in_run(all_merged)

    # ---------- Catalogs & presence ----------
    vf_headers, vf_grouped = _vf_rows_from_merged(all_merged, samples)

    acquired_catalog, mut_pheno_map = _load_amr_maps()    # acquired gene list & mutation phenotype map
    amr_acq_present, amr_acq_pheno, amr_mut_present, amr_mut_pheno = _collect_amr_presence(
        all_merged, acquired_catalog, mut_pheno_map, out_root=out_root, verbose_print=False
    )

    # ---------- Styles & document ----------
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("TitleCentered", parent=styles["Title"], alignment=TA_CENTER, fontName="Helvetica-Bold")
    h2_style = ParagraphStyle("H2", parent=styles["Heading2"], spaceBefore=10, spaceAfter=4)

    note_style = ParagraphStyle(
        "NoteWrap",
        parent=styles["Normal"],
        fontSize=9,
        leading=11,
        alignment=TA_LEFT,
        wordWrap="CJK",
    )

    pheno_style = ParagraphStyle(
        "PhenotypeWrap",
        parent=styles["Normal"],
        fontSize=9,
        leading=11,
        alignment=TA_LEFT,
        wordWrap="CJK",
    )

    rot_hdr_style = ParagraphStyle(
        "RotHdr",
        parent=styles["Normal"],
        fontName="Helvetica-Bold",
        fontSize=9,
        leading=11,
    )

    # Titles (left & right boxes)
    box_title_style = ParagraphStyle(
        "BoxTitle",
        parent=styles["Heading3"],
        alignment=TA_LEFT,
        fontName="Helvetica-Bold",
        fontSize=12,
        spaceAfter=4,
        spaceBefore=0,
    )
    left_body_style = ParagraphStyle(
        "LeftHeader",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=9,
        leading=11,
    )

    doc = SimpleDocTemplate(
        pdf_path, pagesize=landscape(letter),
        leftMargin=0.55*inch, rightMargin=0.55*inch, topMargin=0.55*inch, bottomMargin=0.55*inch
    )

    def _on_page(canvas, doc_):
        canvas.setFont("Helvetica", 9)
        canvas.drawRightString(doc_.pagesize[0] - doc_.rightMargin, 0.45*inch, f"BactiPipe - Vir & AMR  | {run_name}/{accession} | Page {doc_.page}")

    story: List = []
    story.append(Paragraph("WGS - Virulence and AMR Gene Detection", title_style))
    story.append(Spacer(1, 6))

    # ---------- Header band: two columns (left: Run Metadata, right: Tool versions & cutoffs)

    # Left box (Run Metadata)
    left_text = make_report_header_block(
        title=pipeline_title,     # <— pipeline title goes here
        run_name=run_name,
        tech_name=tech_name,
        accession=accession,
    )

    # Convert lines to HTML and prepend an explicit “Run Metadata” heading
    left_html = "<br/>".join(_xml_escape(ln) for ln in left_text.splitlines())
    left_cell = [
        Paragraph("Run Metadata", box_title_style),
        Paragraph(left_html, left_body_style),
    ]

    # --- Right block: Tool versions, databases, and cutoffs ---
    tools = (versions or {}).get("tools", {})
    dbs = (versions or {}).get("databases", {})

    right_lines = ["<b>Tool Version Numbers &amp; Cutoffs</b>", "=" * 30]

    # Tools
    if tools:
        for k, v in sorted(tools.items()):
            right_lines.append(f"{_xml_escape(k)}: {_xml_escape(v)}")

    # Databases
    if dbs:
        right_lines.append("<br/><b>Databases</b>")
        for k, v in sorted(dbs.items()):
            right_lines.append(f"{_xml_escape(k)}: {_xml_escape(v)}")

    # Cutoffs
    if thresholds:
        mi = thresholds.get("min_identity")
        mc = thresholds.get("min_coverage")
        right_lines.append("<br/><b>Cutoffs</b>")
        if mi is not None:
            right_lines.append(f"Minimum identity: {float(mi):.0f}%")
        if mc is not None:
            right_lines.append(f"Minimum coverage: {float(mc):.0f}%")

    right_lines.append("=" * 30)
    right_para = Paragraph("<br/>".join(right_lines), styles["Normal"])

    # --- Compose the two-column band ---
    band = Table([[left_cell, right_para]], colWidths=[6.2*inch, 3.8*inch], hAlign="LEFT")
    band.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    story.append(band)
    story.append(Spacer(1, 10))

    # ---------- Column width helpers (fit to page width) ----------
    def _vf_colwidths(doc_width: float, n_samples: int) -> List[float]:
        """
        Width proportions for the virulence factors table.
        """
        fixed_fracs = [0.10, 0.36]
        fixed_w = [doc_width * f for f in fixed_fracs]
        rem = max(doc_width - sum(fixed_w), 1.0)
        samp_w = [rem / n_samples] * n_samples if n_samples > 0 else []
        return fixed_w + samp_w

    def _amr_colwidths(doc_width: float, n_samples: int) -> List[float]:
        """
        Width proportions for the AMR determinants table.
        """
        fixed_fracs = [0.12, 0.34]
        fixed_w = [doc_width * f for f in fixed_fracs]
        rem = max(doc_width - sum(fixed_w), 1.0)
        samp_w = [rem / n_samples] * n_samples if n_samples > 0 else []
        return fixed_w + samp_w
   
    HEADER_ROW_HEIGHT = 42.0  # fixed header row height
    HEADER_TILT_DEGREES = 35  # preferred tilt

    # ------------------------------------------------------------------
    # Virulence section (from merged rows; wrapped Note; category headers)
    # ------------------------------------------------------------------
    if run_vf:
        story.append(Paragraph("Virulence factors (presence/absence per sample)", h2_style))

        has_vf = any(vf_grouped.values())
        if not has_vf:
            story.append(Paragraph("<i>No virulence factors detected across samples</i>", styles["Normal"]))
        else:
            # Rotated sample headers if >5
            if len(samples) > 5:
                sample_headers = [RotatedLabel(str(s), rot_hdr_style, angle=HEADER_TILT_DEGREES, box_h=HEADER_ROW_HEIGHT) for s in samples]
            else:
                sample_headers = [str(s) for s in samples]

            cols = ["Gene", "Description"] + sample_headers
            data: List[List] = [cols]

            for cat in sorted(vf_grouped.keys(), key=str.lower):
                # Category subheader row spanning whole width
                data.append([f"Category: {cat}"] + [""] * (len(cols) - 1))

                # Rows are [gene, note, marks...]; wrap the note into a Paragraph
                for row in vf_grouped[cat]:
                    gene = row[0]
                    note = Paragraph(_xml_escape(row[1]) if row[1] else "", note_style)
                    marks = row[2:]
                    data.append([gene, note] + marks)

            vf_tbl = _safe_table(
                data,
                style=[
                    ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
                    ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                    ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                    ("ALIGN", (0,0), (-1,0), "CENTER"),   # header row
                    ("ALIGN", (2,1), (-1,-1), "CENTER"),  # sample columns start at col=2
                    ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                    ("FONTSIZE", (0,0), (-1,-1), 9),
                    ("TOPPADDING", (0,0), (-1,-1), 3),
                    ("BOTTOMPADDING", (0,0), (-1,-1), 3),
                ],
                colWidths=_vf_colwidths(doc.width, len(samples)),
            )

            total_cols = len(cols)
            for r in range(1, len(data)):
                if isinstance(data[r][0], str) and data[r][0].startswith("Category:"):
                    vf_tbl.setStyle(TableStyle([
                        ("SPAN", (0, r), (total_cols-1, r)),
                        ("BACKGROUND", (0, r), (total_cols-1, r), colors.whitesmoke),
                        ("FONTNAME", (0, r), (total_cols-1, r), "Helvetica-Bold"),
                        ("ALIGN", (0, r), (total_cols-1, r), "LEFT"),
                    ]))

            story.append(vf_tbl)

        story.append(Spacer(1, 14))

    # ------------------------------------------------------------------
    # AMR section (Acquired vs Point mutations; Phenotype wrapped)
    # ------------------------------------------------------------------
    if run_amr:
        story.append(Paragraph("AMR determinants (presence/absence per sample)", h2_style))

        if not (amr_acq_present or amr_mut_present):
            story.append(Paragraph("<i>No AMR determinants detected across samples</i>", styles["Normal"]))
        else:
            if len(samples) > 5:
                sample_headers = [RotatedLabel(str(s), rot_hdr_style, angle=HEADER_TILT_DEGREES, box_h=HEADER_ROW_HEIGHT) for s in samples]
            else:
                sample_headers = [str(s) for s in samples]

            cols = ["Name", "Resistance Phenotype"] + sample_headers
            data: List[List] = [cols]

            # Acquired genes
            kept_acq = sorted([det for det, smp in amr_acq_present.items() if smp], key=str.lower)
            if kept_acq:
                data.append(["Category: Acquired genes"] + [""]*(len(cols)-1))
                for det in kept_acq:
                    ph = amr_acq_pheno.get(det, "")
                    ph_cell = Paragraph(_xml_escape(ph), pheno_style) if ph else ""
                    checks = ["✓" if s in amr_acq_present[det] else "" for s in samples]
                    data.append([det, ph_cell] + checks)

            # Point mutations
            kept_mut = sorted([mut for mut, smp in amr_mut_present.items() if smp], key=str.lower)
            if kept_mut:
                data.append(["Category: Point mutations"] + [""]*(len(cols)-1))
                for mut in kept_mut:
                    base = _base_token(mut)
                    ph = amr_mut_pheno.get(base, "")
                    ph_cell = Paragraph(_xml_escape(ph), pheno_style) if ph else ""
                    checks = ["✓" if s in amr_mut_present[mut] else "" for s in samples]
                    data.append([mut, ph_cell] + checks)

            if len(data) == 1:
                story.append(Paragraph("<i>No AMR determinants detected across samples</i>", styles["Normal"]))
            else:
                amr_tbl = _safe_table(
                    data,
                    style=[
                        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
                        ("BACKGROUND", (0,0), (-1,0), colors.lightgrey),
                        ("FONTNAME", (0,0), (-1,0), "Helvetica-Bold"),
                        ("ALIGN", (0,0), (-1,0), "CENTER"),   # header
                        ("ALIGN", (2,1), (-1,-1), "CENTER"),  # sample cols start at col=2
                        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
                        ("FONTSIZE", (0,0), (-1,-1), 9),
                        ("TOPPADDING", (0,0), (-1,-1), 3),
                        ("BOTTOMPADDING", (0,0), (-1,-1), 3),
                    ],
                    colWidths = _amr_colwidths(doc.width, len(samples)),
                )

                total_cols = len(cols)
                for r in range(1, len(data)):
                    if isinstance(data[r][0], str) and data[r][0].startswith("Category:"):
                        amr_tbl.setStyle(TableStyle([
                            ("SPAN", (0, r), (total_cols-1, r)),
                            ("BACKGROUND", (0, r), (total_cols-1, r), colors.whitesmoke),
                            ("FONTNAME", (0, r), (total_cols-1, r), "Helvetica-Bold"),
                            ("ALIGN", (0, r), (total_cols-1, r), "LEFT"),
                        ]))
                story.append(amr_tbl)

    # ------------------------------------------------------------------
    doc.build(story, onFirstPage=_on_page, onLaterPages=_on_page)
    return pdf_path

def write_consolidated_tsv(
    all_merged: dict,
    out_root: str,
    filename: str = "traits_consolidated.tsv",
    include_header: bool = True,
) -> str:
    """
    Create a run-wide TSV mirroring the PDF layout.

    Blocks (each separated by a blank line):
      1) Virulence:
         Header: Category   Gene   Note   <sample...>     (if include_header)
         Rows:   category   gene   note   1/blank per sample
      2) AMR Acquired:
         Header: Type   Gene   Phenotype   <sample...>    (Type='Acquired')
      3) AMR Mutations:
         Header: Type   Mutation   Phenotype   <sample...> (Type='Mutation')

    all_merged: { sample_display_name: { "vf": [...], "amr": [...] } }
    """
    sample_names = list(all_merged.keys())

    # ---- Collect virulence keys ----
    vf_order, vf_seen = [], set()
    vf_presence = defaultdict(set)

    # ---- AMR acquired ----
    acq_order, acq_seen = [], set()
    acq_presence = defaultdict(set)

    # ---- AMR mutations ----
    mut_order, mut_seen = [], set()
    mut_presence = defaultdict(set)

    for sname, merged in all_merged.items():
        # VF
        for row in (merged.get("vf") or []):
            cat = row.get("category") or row.get("functional_category") or "Misc"
            gene = row.get("gene") or row.get("name") or ""
            note = row.get("note") or row.get("notes") or ""
            key = (cat, gene, note)
            if key not in vf_seen:
                vf_seen.add(key); vf_order.append(key)
            vf_presence[key].add(sname)

        # AMR acquired. Normalized per-sample data stores acquired and mutation
        # records together under "amr" and distinguishes them with "type".
        for row in (merged.get("amr") or []):
            if _is_mutation_amr(row):
                continue
            gene = row.get("determinant") or row.get("name") or row.get("gene") or ""
            phenotype = row.get("phenotype") or row.get("class") or ""
            key = (gene, phenotype)
            if key not in acq_seen:
                acq_seen.add(key); acq_order.append(key)
            acq_presence[key].add(sname)

        # AMR mutations
        for row in (merged.get("amr") or []):
            if not _is_mutation_amr(row):
                continue
            mut = row.get("determinant") or row.get("mutation") or row.get("name") or ""
            phenotype = row.get("phenotype") or ""
            key = (mut, phenotype)
            if key not in mut_seen:
                mut_seen.add(key); mut_order.append(key)
            mut_presence[key].add(sname)

    out_path = os.path.join(out_root, "reports", filename)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")

        # -------- Virulence block --------
        if vf_order:
            if include_header:
                w.writerow(["Category", "Gene", "Note", *sample_names])
            for (cat, gene, note) in vf_order:
                row = [cat, gene, note]
                row.extend("1" if s in vf_presence[(cat, gene, note)] else "" for s in sample_names)
                w.writerow(row)
            w.writerow([])  # spacer

        # -------- AMR acquired --------
        if acq_order:
            if include_header:
                w.writerow(["Type", "Gene", "Phenotype", *sample_names])
            for (gene, phen) in acq_order:
                row = ["Acquired", gene, phen]
                row.extend("1" if s in acq_presence[(gene, phen)] else "" for s in sample_names)
                w.writerow(row)
            w.writerow([])

        # -------- AMR mutations --------
        if mut_order:
            if include_header:
                w.writerow(["Type", "Mutation", "Phenotype", *sample_names])
            for (mut, phen) in mut_order:
                row = ["Mutation", mut, phen]
                row.extend("1" if s in mut_presence[(mut, phen)] else "" for s in sample_names)
                w.writerow(row)

    return out_path
