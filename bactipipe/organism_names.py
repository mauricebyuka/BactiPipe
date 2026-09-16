"""Shared expected-organism names, without fuzzy taxonomy or strain confirmation.

Normalization changes comparison resolution, not the original sample-sheet text.
It is deliberately independent of the QC genome-size catalog: identification may
resolve organisms outside that list. Catalog membership is checked separately.
"""
import re
import unicodedata


def normalize_organism_name(value: object) -> str:
    text = unicodedata.normalize("NFKC", str(value or ""))
    return re.sub(r"\s+", " ", text).strip().casefold()


# Historical name reviewed for species-level expectation only, not serotyping.
ALIASES = {"salmonella typhimurium": "salmonella enterica"}
# Common unlabelled serovars in laboratory sheets. Other serovars can use an
# explicit "serovar"/"serotype" label; never guess an arbitrary trailing word.
SALMONELLA_SEROVARS = frozenset({
    "typhimurium", "enteritidis", "dublin", "newport", "heidelberg", "kentucky",
    "infantis", "montevideo",
})


def expected_organism_name(value: object) -> str:
    """Strip recognizable below-species qualifiers; retain explicit subspecies.

Keep unknown suffixes, alternative/mixed names, and uncertainty qualifiers
unchanged. No abbreviation expansion, spelling repair, or two-word truncation.
"""
    name = normalize_organism_name(value)
    if name in ALIASES:
        return ALIASES[name]
    if re.search(r"[/+&;]|\b(?:and|or|complex|group|mixed|cf|aff)\b", name):
        return name
    match = re.fullmatch(
        r"(?P<species>[a-z][a-z-]* [a-z][a-z-]*)"
        r"(?: (?:(?:subsp|ssp)\.?|subspecies) (?P<subspecies>[a-z][a-z-]*))?"
        r"(?P<tail> .+)?", name,
    )
    if not match:
        return name
    species = match['species']
    subspecies = match['subspecies']
    base = species + (f" subsp. {subspecies}" if subspecies else "")
    tail = (match['tail'] or '').strip()
    if not tail:
        return base
    if tail.startswith('(') and tail.endswith(')'):
        tail = tail[1:-1].strip()
    # An unparsed subspecies designation must never be erased as a strain.
    if re.search(r"\b(?:subsp|ssp|subspecies)\b", tail):
        return name
    if species == "salmonella enterica":
        if tail in SALMONELLA_SEROVARS:
            return base
        for serovar in SALMONELLA_SEROVARS:
            if tail.startswith(serovar + ' '):
                tail = tail[len(serovar):].strip()
                break
    if re.fullmatch(r"(?:strain|str\.?|isolate|serovar|serotype|ser\.?|sv\.?):?\s+\S.*", tail):
        return base
    if re.fullmatch(r"(?:atcc|nctc|dsm|dsmz|ccug|cip)\s+\d+[a-z0-9.-]*", tail):
        return base
    # Recognizable unlabelled laboratory strain/serotype codes, e.g. K-12,
    # ATCC25922 or O157:H7. Pure numbers/unknown words alone stay unrecognized.
    if re.fullmatch(r"(?=[a-z0-9:_.-]*[a-z])(?=[a-z0-9:_.-]*\d)[a-z0-9][a-z0-9:_.-]*", tail):
        return base
    return name
