"""Apply verified species confirmation without changing measured QC values."""
import csv
import io
from .qc_status import interpret, threshold_verdict


def update_row(row, identification, *, platform='', quality=15, completeness=90,
               contamination=10, selected_mean_quality=None):
    row = dict(row)
    if identification.get('status') not in {'confirmed_species', 'no_confident_identification'}:
        return row
    confirmed = identification['status'] == 'confirmed_species'
    row['Identified organism'] = identification['reported_organism'] if confirmed else 'No confident identification'
    row['% Match'] = str(identification['match_percent']) if confirmed else 'N/A'
    if confirmed and platform.lower() == 'illumina':
        from .scripts.qc_plan import illumina_minimum_coverage
        if selected_mean_quality is not None:
            row['Required depth'] = str(illumina_minimum_coverage(
                selected_mean_quality, identification['reported_organism'])) + 'X'
        elif identification['reported_organism'].casefold() == 'escherichia coli':
            # Do not infer an unrounded quality tier from the displayed mean.
            row['Required depth'] = 'N/A'
    def value(key):
        return str(row.get(key, '')).rstrip('Xx')
    grades = [threshold_verdict(value('Mean quality'), quality),
              threshold_verdict(value('Coverage Depth'), value('Required depth')),
              threshold_verdict(value('CheckM completeness'), completeness),
              threshold_verdict(value('CheckM contamination'), contamination, maximum=True)]
    verdict = interpret(row.get('Expected organism', ''), row['Identified organism'], grades,
        isolate_qc_failure=str(row.get('Isolate QC', '')).lower() == 'fail',
        isolate_qc_review=str(row.get('Isolate QC', '')).lower() == 'review')
    row['Overall Quality'] = verdict['overall']
    return row


def update_summary(text, confirmations, **kwargs):
    lines = text.splitlines(True)
    index = next(i for i, line in enumerate(lines) if line.startswith('Sample ID\t'))
    reader = csv.DictReader(lines[index:], delimiter='\t')
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=reader.fieldnames, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in reader:
        entry = confirmations.get(row.get('Sample ID'))
        if entry:
            row = update_row(row, entry['identification'],
                selected_mean_quality=entry.get('selected_mean_quality'), **kwargs)
        writer.writerow(row)
    return ''.join(lines[:index]) + output.getvalue()
