"""Translate optional primary-workbook requests into the downstream contract.

No reference database or Django dependency: reference names are resolved by the
caller at execution time. Blank request columns never request an analysis.
"""
from __future__ import annotations

import csv
from hashlib import sha256
from io import StringIO
import re

from .downstream_manifest import DownstreamManifestError, _parse_rows
from .organism_names import expected_organism_name
from .scripts.caps import CAPS


RELATE_NAMES = {
    'acinetobacter baumannii': 'abaumannii',
    'clostridium perfringens': 'clostridium',
    'campylobacter fetus': 'cfetus', 'campylobacter jejuni': 'cjejuni',
    'escherichia coli': 'ecoli', 'klebsiella pneumoniae': 'klebsiella',
    'listeria monocytogenes': 'listeria', 'pasteurella multocida': 'pmultocida',
    'salmonella enterica': 'salmonella',
    'staphylococcus epidermidis': 'sepidermidis',
    'staphylococcus aureus': 'saureus',
    'staphylococcus pseudintermedius': 'spseudintermedius',
    'streptococcus equi subsp. zooepidemicus': 'szooepidemicus',
    'ornithobacterium rhinotracheale': 'orhinotracheale',
}
HEADERS = ('cohort_id', 'command', 'organism', 'other_organism', 'sample_id',
           'display_name', 'specimen', 'role', 'reference_name', 'run_ani',
           'run_ska', 'run_cgmlst', 'detect_amr', 'detect_virulence',
           'extended_virulence', 'abricate_resfinder', 'abricate_card')
EMPTY = {'', 'none', 'n/a', 'na', '-', '--'}


def text(value):
    return str(value).strip() if value is not None else ''


def report_order(rows, *, name_index=1):
    """Named rows first (case-insensitive); unnamed rows retain worksheet order."""
    return sorted(rows, key=lambda r: (not bool(text(r[name_index])), text(r[name_index]).casefold()))


def organism_token(entered, command, detect_organisms=()):
    name = expected_organism_name(entered)
    if command == 'relate':
        token = RELATE_NAMES.get(name)
        if token and token in CAPS:
            return token, ''
    else:
        lookup = {value.replace('_', ' ').casefold(): value for value in detect_organisms}
        token = lookup.get(name) or lookup.get(name.split(' ')[0])
        if token:
            return token, ''
    # Unsupported species must not borrow a scheme from a different species.
    # Both CLIs support generic analysis for an unlisted genus.
    genus = name.split(' ')[0] if name else 'Unknown'
    return 'Other', genus.capitalize()


def primary_downstream_tsv(path, *, detect_organisms=()):
    """Return a validated manifest TSV, or '' when no analysis was requested."""
    from openpyxl import load_workbook
    with_sheet = load_workbook(path, read_only=True, data_only=True)
    try:
        sheet = with_sheet.worksheets[0]
        # Keep the sample-table boundary; never ingest protocol/kit notes below it.
        rows = iter(sheet.iter_rows(values_only=True))
        prefix = []
        columns = None
        for row_number, row in enumerate(rows, 1):
            prefix.append(row)
            names = [re.sub(r'\s+', ' ', text(v)).casefold() for v in row]
            if 'accession number' in names and 'bacteria species' in names:
                if 'r/d*' not in names and 'r/d' not in names:
                    return ''
                columns = {v: i for i, v in enumerate(names) if v}
                break
            if row_number >= 40:
                return ''
        if columns is None:
            return ''
        refs = {}
        for i in range(4):
            row = prefix[4+i] if len(prefix) > 4+i else ()
            label = text(row[4] if len(row) > 4 else '').upper().rstrip(':')
            if label not in {'', f'R{i+1}'}:
                raise DownstreamManifestError(f'E{i+5}: expected R{i+1}, found {label!r}.')
            value = text(row[5] if len(row) > 5 else '')
            refs[f'R{i+1}'] = list(dict.fromkeys(
                v.strip() for v in re.split(r'[,;\n]+', value)
                if v.strip().casefold() not in EMPTY
            ))
        grouped = {}
        started = False
        seen = set()
        for row in rows:
            row_number += 1
            def cell(name):
                i = columns.get(name)
                return text(row[i]) if i is not None and i < len(row) else ''
            sid, entered = cell('accession number'), cell('bacteria species')
            if started and not sid and not entered:
                break
            if not sid and not entered:
                continue
            started = True
            request = cell('r/d*') or cell('r/d')
            if request.casefold() in EMPTY:
                continue
            if not sid or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', sid):
                raise DownstreamManifestError(f'Row {row_number}: invalid Accession Number {sid!r}.')
            if sid in seen:
                raise DownstreamManifestError(f'Row {row_number}: duplicate sample {sid!r}.')
            seen.add(sid)
            requests = list(dict.fromkeys(re.split(r'[,;/+\s]+', request.upper())))
            if any(v not in {'R', 'R1', 'R2', 'R3', 'R4', 'D'} for v in requests):
                raise DownstreamManifestError(f'Row {row_number} ({sid}): use R1–R4, R or D; combine with commas.')
            if entered.casefold() == 'lambda':
                raise DownstreamManifestError(f'Row {row_number}: Lambda cannot request bacterial downstream analysis.')
            if entered.casefold() == 'other':
                entered = cell('other bacteria species') or entered
            accession = sid.rsplit('-', 1)[0] if '-' in sid else sid
            for command in dict.fromkeys('detect' if req == 'D' else 'relate' for req in requests):
                organism, other = organism_token(entered, command, detect_organisms)
                references = {}
                if command == 'relate':
                    for req in requests:
                        for reference in refs.get(req, []):
                            references.setdefault(reference.casefold(), reference)
                reference_keys = tuple(sorted(references))
                key = (accession, command, organism, other, reference_keys)
                group = grouped.setdefault(key, {'references': [references[k] for k in reference_keys], 'rows': []})
                record = dict(command=command, organism=organism, other_organism=other,
                              sample_id=sid, display_name=cell('isolate id'),
                              specimen=cell('host and source'), role='sample')
                if command == 'relate':
                    record.update(run_ani='yes', run_ska='yes',
                                  run_cgmlst='yes' if CAPS.get(organism, {}).get('cgmlst') else 'no')
                else:
                    record.update(detect_amr='yes', detect_virulence='yes',
                                  extended_virulence='no', abricate_resfinder='no', abricate_card='no')
                group['rows'].append(record)
        if not grouped:
            return ''
        output = StringIO()
        writer = csv.DictWriter(output, fieldnames=HEADERS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key, group in grouped.items():
            accession, command, organism, other, _ = key
            # Stable identity per organism/reference, independent of row order or other cohorts.
            suffix = sha256(repr(key).encode()).hexdigest()[:10]
            cohort_id = f'{accession[:48]}-{command}-{suffix}'
            records = group['rows']
            if command == 'relate' and len(records) == 1 and not group['references']:
                records[0]['run_ska'] = 'no'  # Pairwise analysis needs two genomes.
            used_ids = {r['sample_id'] for r in records}
            for ref in group['references']:
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', ref):
                    raise DownstreamManifestError(f'Invalid reference name {ref!r}; separate names with commas, semicolons or newlines.')
                ref_sid = ref
                while ref_sid in used_ids:
                    ref_sid = 'REF-' + ref_sid
                used_ids.add(ref_sid)
                writer.writerow(dict(cohort_id=cohort_id, command=command, organism=organism,
                                     other_organism=other, sample_id=ref_sid, display_name='reference',
                                     role='reference', reference_name=ref))
            for record in records:
                writer.writerow(dict(cohort_id=cohort_id, **record))
        value = output.getvalue()
        _parse_rows(csv.reader(StringIO(value), delimiter='\t'))
        return value
    finally:
        with_sheet.close()
