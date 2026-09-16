import tempfile
from pathlib import Path
import unittest
import json
from unittest.mock import patch
from bactipipe.reproducibility import database_log_record

from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.platypus import Table
from bactipipe.scripts.pdf_renderer import (
    LayoutConfig, _styles, _build_dual_header_band, _report_header_groups, render_pdf_from_rows,
)


ENTRIES=[('BactiPipe','0.1.0'),('BactiPipe','v0.1.0'),('Filtlong','Filtlong v0.2.1'),
         ('Flye','2.9.6-b1802'),('SpeciesFinder','v4.0.4'),('CheckM','1.2.4'),
         ('BactiPipe database','BactiPipe-DB-v1'),
         ('Species identification database','BactiPipe-SpeciesFinder-v4.0.4-20260112-3a72a1e6a9c7'),
         ('Minimum read length','1000 bp'),('Required Mean quality','15'),
         ('BLASTn','2.16.0+'), ('Minimum Lambda identity','99%'),
         ('Minimum Lambda genome coverage','90%'),
         ('Completeness cutoff','90%'),('Contamination cutoff','10%')]


class PDFHeaderTests(unittest.TestCase):
    def test_minimum_depth_immediately_precedes_mean_quality(self):
        for value in ('50X', '30X (Q ≥30); 40X (29 ≤Q <30); 50X (Q <29); identified E. coli: +10X'):
            entries = ENTRIES + [('Minimum depth', value)]
            thresholds = _report_header_groups(entries)[2]
            quality = thresholds.index(('Minimum mean quality', '15'))
            self.assertEqual(thresholds[quality-1], ('Minimum depth', value))
            band = _build_dual_header_band(['Platform: Illumina'], entries, LayoutConfig(), _styles(LayoutConfig()))
            table = next(x for x in band._cellvalues[0][2] if isinstance(x, Table))
            if ';' in value:
                self.assertIn(('SPAN', (0, quality-1), (1, quality-1)), table._spanCmds)

    def test_illumina_policy_and_control_render_together(self):
        headers = ['Sample ID','Mean quality','Expected organism','Identified organism','% Match','Coverage Depth',
                   'Required depth','CheckM completeness','CheckM contamination','Overall Quality']
        rows = [['Control01','31','Lambda','Lambda','99.99','100X','30X','N/A','N/A','Pass'],
                ['S01','31','Escherichia coli','Escherichia coli','99.8','100X','40X','99','0.1','Pass']]
        policy = '30X (Q ≥30); 40X (29 ≤Q <30); 50X (Q <29); identified E. coli: +10X'
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / 'illumina.pdf'
            render_pdf_from_rows(rows, headers, str(output), run_metadata=['Platform: Illumina'],
                tools_header=ENTRIES + [('Minimum depth', policy)],
                lambda_evidence={'Control01': dict(identity_percent=99.99, reference_coverage_percent=99,
                                   minimum_identity_percent=99, minimum_coverage_percent=90)})
            self.assertTrue(output.read_bytes().startswith(b'%PDF'))

    def test_tools_heading_uses_embedded_legible_font(self):
        from reportlab.pdfbase import pdfmetrics
        style = _styles(LayoutConfig())["HeaderRight"]
        self.assertEqual(style.fontName, "BactiPipeVeraBold")
        font = pdfmetrics.getFont(style.fontName)
        self.assertNotEqual(font.face.charToGlyph[ord("i")], font.face.charToGlyph[ord("l")])
        band=_build_dual_header_band(['Run name: Example'],ENTRIES,LayoutConfig(),_styles(LayoutConfig()))
        self.assertEqual({panel[0].style.fontName for panel in band._cellvalues[0]}, {'BactiPipeVeraBold'})

    def test_lambda_identity_and_coverage_are_separately_colored(self):
        import pandas as pd
        from bactipipe.scripts.pdf_renderer import _build_table
        rows=pd.DataFrame([{'Sample ID':'L', 'Expected organism':'Lambda', '% Match':'99.99'},
                           {'Sample ID':'S', 'Expected organism':'Escherichia coli', '% Match':'99.95'}])
        evidence={'L':dict(identity_percent=99.99,reference_coverage_percent=83.17,
                           minimum_identity_percent=99,minimum_coverage_percent=90)}
        table=_build_table(rows,LayoutConfig(),evidence)
        cell=table._cellvalues[1][2]
        self.assertIn('99.99%',cell.getPlainText())
        self.assertIn('Coverage: 83.17%',cell.getPlainText())
        self.assertEqual(cell.text.count('<font color='),1)
        self.assertIn('>83.17%</font>',cell.text)
        self.assertEqual(table._cellvalues[2][2].getPlainText(),'99.95')

    def test_lambda_missing_coverage_is_not_inferred_from_identity(self):
        from bactipipe.scripts.pdf_renderer import _lambda_match_cell
        cell=_lambda_match_cell('100',{},_styles(LayoutConfig()),LayoutConfig())
        self.assertIn('Coverage: N/A',cell.getPlainText())
        self.assertNotIn('<font color=',cell.text)
        cell=_lambda_match_cell('100',dict(identity_percent=98,reference_coverage_percent=100,
                     minimum_identity_percent=99,minimum_coverage_percent=90),_styles(LayoutConfig()),LayoutConfig())
        self.assertIn('>98.00%</font>',cell.text)
        self.assertEqual(cell.text.count('<font color='),1)

    def test_lambda_annotation_stays_on_two_lines_at_production_width(self):
        import pandas as pd
        from reportlab.lib.units import inch
        from bactipipe.scripts.pdf_renderer import _build_table
        headers=['Sample ID','Mean quality','Expected organism','Identified organism','% Match','Coverage Depth',
                 'Required depth','CheckM completeness','CheckM contamination','Overall Quality']
        row=['Lambda','25','Lambda','Lambda','100','300X','50X','N/A','N/A','Pass']
        layout=LayoutConfig(reserved_width=1.7*inch,last3_narrow_width=.8*inch)
        table=_build_table(pd.DataFrame([row],columns=headers),layout,
                           {'Lambda':dict(identity_percent=100,reference_coverage_percent=100)})
        cell=table._cellvalues[1][4]
        self.assertLessEqual(cell.wrap(table._argW[4]-6,1000)[1],20)

    def test_pdf_reads_retained_control_metrics_without_changing_tsv(self):
        from bactipipe.scripts.pdf_renderer import render_pdf_from_tsv
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            tsv=root/'run_qc_summary.tsv'
            original='Sample ID\tExpected organism\t% Match\nL\tLambda\t100\n'
            tsv.write_text(original)
            evidence=root/'qc_out'/'L'/'lambda_control'/'result.json'
            evidence.parent.mkdir(parents=True)
            evidence.write_text(json.dumps(dict(identity_percent=100,reference_coverage_percent=99.5,
                           minimum_identity_percent=99,minimum_coverage_percent=90)))
            with patch('bactipipe.scripts.pdf_renderer._render_pdf_common') as render:
                render_pdf_from_tsv(str(tsv),str(root/'test.pdf'))
            self.assertEqual(render.call_args.kwargs['lambda_evidence']['L']['reference_coverage_percent'],99.5)
            self.assertEqual(tsv.read_text(),original)

    def test_separate_software_and_database_identities(self):
        software,datasets,thresholds=_report_header_groups(ENTRIES)
        self.assertEqual(software.count(('BactiPipe','v0.1.0')),1)
        self.assertIn(('SpeciesFinder','v4.0.4'),software)
        self.assertIn(('Filtlong','v0.2.1'),software)
        self.assertIn(('Flye','v2.9.6-b1802'),software)
        self.assertEqual(datasets,[('SpeciesFinder db','v4.0.4-20260112-3a72a1e6a9c7')])
        self.assertIn(('CheckM','v1.2.4'),software)
        self.assertNotIn(('BactiPipe database','BactiPipe-DB-v1'),software+datasets+thresholds)
        self.assertNotIn(('Minimum read length','1000 bp'),software+datasets+thresholds)
        self.assertIn(('Minimum completeness','90%'),thresholds)
        self.assertIn(('Maximum contamination','10%'),thresholds)
        self.assertIn(('Minimum read length','1000 bp'),ENTRIES) # inputs/provenance unchanged

    def test_explicit_widths_and_no_split_tool_names(self):
        layout=LayoutConfig();styles=_styles(layout)
        band=_build_dual_header_band(['Run name: Example'],ENTRIES,layout,styles)
        self.assertEqual(len(band._cellvalues[0]),3)
        tables=[x for panel in band._cellvalues[0][1:] for x in panel if isinstance(x,Table)]
        self.assertEqual(len(tables),2)
        self.assertEqual(band._cellvalues[0][1][0].getPlainText(),'Tool Versions')
        self.assertEqual(band._cellvalues[0][2][0].getPlainText(),'QC Thresholds')
        for panel, expected in zip(band._cellvalues[0], (2,2,2)):
            separators=[x for x in panel if hasattr(x,'getPlainText') and x.getPlainText().startswith('===')]
            self.assertEqual(len(separators),expected)
        for i,table in enumerate(tables,1):
            self.assertLessEqual(sum(table._argW),band._argW[i] - (12 if i==1 else 0))
            self.assertFalse(any(not k.getPlainText() and not v.getPlainText() for k,v in table._cellvalues))
        for panel in band._cellvalues[0]:
            for index, flow in enumerate(panel):
                if isinstance(flow,Table):
                    line=panel[index-1]
                    width=stringWidth(line.getPlainText(),line.style.fontName,line.style.fontSize)
                    self.assertLess(abs(width-sum(flow._argW)),10)
        for table in tables:
            for label,value in table._cellvalues:
                self.assertGreater(table._argW[0],stringWidth(label.getPlainText(),'Helvetica',9))
                self.assertFalse(label.style.splitLongWords)
                self.assertLessEqual(label.wrap(table._argW[0],1000)[1],11)

    def test_unknown_version_not_invented(self):
        self.assertEqual(_report_header_groups([('CheckM','v?')])[0],[('CheckM','Not recorded')])

    def test_source_database_identifier_takes_precedence(self):
        entries=ENTRIES+[('SpeciesFinder db','bacteria_14_AT_110126')]
        self.assertEqual(_report_header_groups(entries)[1],[('SpeciesFinder db','bacteria_14_AT_110126')])

    def test_other_databases_are_not_reported_as_software_or_speciesfinder(self):
        entries=[('CheckM database','2015_01_16'),('Reference bundle','BactiPipe-DB-v1'),
                 ('Species identification database','BactiPipe-Species-R232-candidate3')]
        self.assertEqual(_report_header_groups(entries),([],[],[]))

    def test_detailed_database_identities_remain_in_log(self):
        bundle={'version':'BactiPipe-DB-v1','components':{'checkm':'2015_01_16'}}
        species={'version':'candidate','components':{'index':{'sha256':'abc'}}}
        with patch('bactipipe.reproducibility.database_identity',side_effect=[bundle,species]):
            record=database_log_record('/references/species')
        self.assertIn(json.dumps(bundle,sort_keys=True),record)
        self.assertIn(json.dumps(species,sort_keys=True),record)

    def test_pdf_build_with_many_rows_and_long_reference_id(self):
        headers=['Sample ID','Mean quality','Expected organism','Identified organism','% Match','Coverage Depth',
                 'Required depth','CheckM completeness','CheckM contamination','Overall Quality']
        row=['Sample01','20','Escherichia coli','Escherichia coli','99.9','120X','50X','99','0.1','Pass']
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'report.pdf'
            render_pdf_from_rows([row]*48,headers,str(path),run_metadata=['Run name: Example','Analysis by: A & B'],
                                 tools_header=ENTRIES,run_name='Example')
            self.assertTrue(path.read_bytes().startswith(b'%PDF'))
            self.assertGreater(path.stat().st_size,3000)


if __name__=='__main__':unittest.main()
