"""Shared controlled-workbook sample-table boundaries (QC and downstream)."""
import re


def is_footer(accession, organism=''):
    """Recognize documented worksheet sections, not missing biological metadata.

    Blank organisms and absent barcodes remain valid sample inputs. Unknown text
    is not silently discarded: the caller validates its accession separately.
    """
    accession=re.sub(r'\s+',' ',str(accession or '')).strip()
    organism=re.sub(r'\s+',' ',str(organism or '')).strip()
    return bool(re.match(r'^\*?\s*Analysis type\s*:',accession,re.I) or
                re.match(r'^(?:Native Barcoding Kit(?: V\d+)?|Sequencing Auxiliary Vials(?: V\d+)?|'
                         r'Flow Cell|BSA|MinKnow Software Version)\s*:',accession,re.I) or
                (not accession and re.fullmatch(r'Lot\s*#\s*\(Kit\)',organism,re.I)))


def valid_accession(value):
    return bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*',value))
