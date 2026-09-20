"""Exact AMR organism matching and reviewed VirulenceFinder database selection."""
import re


VF_DATABASES = {
    "Escherichia": "virulence_ecoli,stx",
    "Enterococcus": "virulence_ent",
    "Enterococcus_faecalis": "virulence_ent",
    "Enterococcus_faecium": "virulence_ent,virulence_entfm_entls",
    "Enterococcus_lactis": "virulence_ent,virulence_entfm_entls",
    "Listeria": "listeria",
    "Staphylococcus_aureus": "s.aureus_exoenzyme,s.aureus_hostimm,s.aureus_toxin",
}


def organism_key(value):
    return "_".join(str(value or "").split()).casefold()


SIMPLE_NAMES = {
    'salmonella_enterica': 'Salmonella',
    'escherichia_coli': 'Escherichia',
}


def display_organism_name(value):
    """Readable input names; preserve species distinctions where tools need them."""
    return SIMPLE_NAMES.get(organism_key(value), str(value).replace('_', ' '))


def amrfinder_organism(value, supported):
    # No genus guessing: Enterococcus must not become E. faecalis/faecium.
    lookup = {organism_key(item): item for item in supported}
    exact = lookup.get(organism_key(value))
    if exact:
        return exact
    # Only the reviewed legacy group aliases may bridge a species/group name.
    wanted = organism_key(display_organism_name(value))
    matches = [item for item in supported if organism_key(display_organism_name(item)) == wanted]
    return matches[0] if len(matches) == 1 else None


def virulence_databases(value):
    key = organism_key(value)
    exact = {organism_key(name): db for name, db in VF_DATABASES.items()}
    if key in exact:
        return exact[key]
    # These databases are genus-wide; Staphylococcus is intentionally not.
    for genus in ("enterococcus", "escherichia", "listeria"):
        if key.startswith(genus + "_"):
            return exact[genus]
    return None


def parse_amrfinder_organisms(stdout, stderr=""):
    """Only parse an explicit list or a plain one-token-per-line list, never logs."""
    text = stdout + "\n" + stderr
    for prefix in ("Available --organism options:", "Possible organisms:"):
        if prefix in text:
            payload = text.split(prefix, 1)[1].splitlines()
            lines = []
            for line in payload:
                line = line.strip()
                if not line:
                    if lines:
                        break
                    continue
                if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_, ]*", line):
                    break
                tokens = [item.strip() for item in line.split(",") if item.strip()]
                if not all(re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", item) for item in tokens):
                    break
                lines.extend(tokens)
            if lines:
                return sorted(set(lines))
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if lines and all(re.fullmatch(r"[A-Z][A-Za-z0-9_]*", line) for line in lines):
        return sorted(set(lines))
    raise RuntimeError("Cannot read AMRFinderPlus organism options; check its installation/database.")
