"""Shared interpretation of reference-qualified relatedness table columns."""


def reference_column_groups(headers):
    """Return (first comparison column, reference IDs) for paired v2 columns."""
    for start, header in enumerate(headers):
        if not header.endswith(': SNP distance'):
            continue
        refs = []
        for i in range(start, len(headers), 2):
            ref = headers[i].removesuffix(': SNP distance')
            if i + 1 >= len(headers) or headers[i:i+2] != [ref + ': SNP distance', ref + ': ANI (%)']:
                return len(headers), []
            refs.append(ref)
        return start, refs
    return len(headers), []


def report_sections(sections, references_per_page=2):
    """Keep wide PDF tables readable; repeat common columns for later refs."""
    for title, rows in sections:
        start, refs = reference_column_groups(rows[0]) if rows else (0, [])
        if len(refs) <= references_per_page:
            yield title, rows, False
            continue
        for offset in range(0, len(refs), references_per_page):
            indices = list(range(start)) + list(range(start + 2*offset, min(len(rows[0]), start + 2*(offset+references_per_page))))
            selected = [[row[i] if i < len(row) else '' for i in indices] for row in rows]
            yield title + ' — ' + ', '.join(refs[offset:offset+references_per_page]), selected, offset > 0
