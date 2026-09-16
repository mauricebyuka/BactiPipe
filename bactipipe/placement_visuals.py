"""Presentation-only tree and ANI figures; never executes identification tools.

Design adapted from the user's P2631472 diagnostic R figures: highlighted query,
nearest-reference context, and ANI/AF threshold guides. Unlike those case-specific
scripts, all labels and values here come from the supplied evidence. No bootstrap,
type-strain, outgroup, or diagnostic conclusion is inferred by the renderer.
"""
from copy import deepcopy
import math
from pathlib import Path
import re
import textwrap

STYLE_VERSION = 'placement-visuals-v6-compact-tree-label'
INK = '#243746'
MUTED = '#64748b'
RED = '#bd343f'
BLUE = '#2463a5'
GREEN = '#287749'


def _style():
    import matplotlib
    matplotlib.use('Agg')
    return matplotlib.rc_context({'font.family':'DejaVu Sans', 'font.size':10,
        'text.parse_math':False, 'svg.fonttype':'none', 'svg.hashsalt':STYLE_VERSION,
        'axes.titleweight':'bold', 'figure.facecolor':'white', 'savefig.facecolor':'white'})


def _save(figure, destination):
    import matplotlib.pyplot as plt
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=170, bbox_inches='tight', pad_inches=.20,
                   metadata={'Creator':'BactiPipe'})
    plt.close(figure)


def _title(figure, title, subtitle):
    figure.text(.035,.965,'BACTIPIPE  /  SPECIES INVESTIGATION', color=MUTED, fontsize=9, weight='bold')
    figure.text(.035,.916,title, color=INK, fontsize=19, weight='bold')
    figure.text(.035,.873,subtitle, color=MUTED, fontsize=10)


def tip_label(label):
    match = re.search(r'\[(GC[AF]_\d+\.\d+)\]$', label)
    acc = match[1] if match else ''
    name = label[:match.start()].strip() if match else label
    cluster = ''
    if ' (GTDB: ' in name and name.endswith(')'):
        name, cluster = name[:-1].split(' (GTDB: ', 1)
    return name, acc, cluster


def named_reference(name):
    """Figure-only exclusion; never removes competitors from identification."""
    name = str(name or '').strip()
    return bool(name and not re.search(
        r'^Unresolved reference\b|\b(?:uncultured|unclassified|bacterium|archaeon)\b|\bsp(?:\.|\d|\s|$)',
        name, re.I))


def prepare_tree(tree, query_name, maximum_references=10, reference_names=None):
    """Prune only when needed; retain branch lengths/topology without rerooting."""
    from .phylogenetic_placement import _distances_from_query, _induced_tree
    tree = deepcopy(tree)
    leaves = tree.get_terminals()
    if len({tip.name for tip in leaves}) != len(leaves):
        raise ValueError('Tree tips must be unique')
    query = next((tip for tip in leaves if tip.name == query_name), None)
    if query is None:
        raise ValueError('The exact sample leaf is absent from the tree')
    if any(clade.branch_length is not None and (not math.isfinite(clade.branch_length) or clade.branch_length < 0)
           for clade in tree.find_clades()):
        raise ValueError('Invalid branch length in saved tree')
    distance = _distances_from_query(tree, query)
    neighbors = sorted((tip for tip in leaves if tip is not query), key=lambda tip:(distance[tip],tip.name or ''))
    if not neighbors:
        raise ValueError('Tree has no reference tips')
    limit = max(1,min(10,int(maximum_references)))
    nearest_name = neighbors[0].name
    if reference_names is not None:
        names = set(reference_names) & {tip.name for tip in neighbors}
        if not names:
            raise ValueError('The shared panel has no references in this tree')
        nearest_name = next(tip.name for tip in neighbors if tip.name in names)
        tree = _induced_tree(tree, {query_name, *names})
    elif len(neighbors) > limit:
        tree = _induced_tree(tree, {query_name, *[tip.name for tip in neighbors[:limit]]})
    tree.ladderize()
    return tree, nearest_name


def shared_reference_panel(evidence, newick=None, *, query_name='query', outgroup=None, maximum_references=10):
    """Bounded union joined ONLY by exact versioned assembly accession.

    Reserve up to five places for ANI evidence, four for nearest marker neighbors,
    and one for the supplied outgroup. Reference IDs are identifiers, not ranks.
    Historical trees cannot acquire new tips without biological placement.
    """
    from Bio import Phylo
    from .phylogenetic_placement import _distances_from_query
    limit = max(3, min(11, int(maximum_references)))
    measured, _, _ = ani_rows(evidence, maximum_references=10)
    measured = [row for row in measured if named_reference(row.get('species_name'))]
    # Include exact measurements for tree references even if a different strain
    # of their species was selected by the ANI figure's representative ranking.
    all_measured = {}
    for value in evidence.get('reference_audit') or evidence.get('closest_species', []):
        if not named_reference(value.get('species_name')):
            continue
        rows, _, _ = ani_rows({'reference_audit':[value],
                              'identification_thresholds': evidence.get('identification_thresholds', {})})
        if rows and re.fullmatch(r'GC[AF]_\d+\.\d+', str(value.get('accession',''))):
            all_measured[value['accession']] = rows[0]
    tree_refs = []
    if newick:
        tree = Phylo.read(str(newick), 'newick')
        query = next((tip for tip in tree.get_terminals() if tip.name == query_name), None)
        if query is None:
            raise ValueError('The exact sample leaf is absent from the tree')
        distances = _distances_from_query(tree, query)
        for tip in sorted((t for t in tree.get_terminals() if t is not query),
                          key=lambda t:(distances[t],t.name or '')):
            name, acc, cluster = tip_label(tip.name or '')
            if acc and (named_reference(name) or acc == (outgroup or {}).get('reference_accession')):
                tree_refs.append({'accession':acc, 'species_name':name,
                    'tree_label':tip.name, 'tree_distance':distances[tip], 'cluster':cluster})
    trees = {row['accession']:row for row in tree_refs}
    selected = []
    def add(acc):
        if acc and acc not in selected and len(selected)<limit:
            selected.append(acc)
    out_acc = (outgroup or {}).get('reference_accession')
    if out_acc in trees:
        add(out_acc)
    # Keep the nearest marker neighbor even if ANI candidates fill their quota.
    if tree_refs:
        add(tree_refs[0]['accession'])
    for row in measured[:5]:
        if row.get('accession') in all_measured:
            add(row['accession'])
    for row in tree_refs:
        add(row['accession'])
    for row in measured:
        if row.get('accession') in all_measured:
            add(row['accession'])
    panel = []
    for index, acc in enumerate(sorted(selected), 1):
        numeric = all_measured.get(acc, {})
        row = {**numeric, **trees.get(acc, {}), 'accession':acc,
               'panel_id':f'R{index:02d}', 'has_ani':bool(numeric), 'in_tree':acc in trees,
               'is_outgroup':acc == out_acc}
        row.setdefault('species_name', acc)
        row['availability'] = ('Tree and ANI' if numeric and acc in trees else
                               'ANI comparison only' if numeric else 'ANI not measurable')
        panel.append(row)
    return panel


def render_tree(newick, destination, *, query_name='query', sample_name=None, maximum_references=10, outgroup=None, panel=None):
    from Bio import Phylo
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    if panel is not None and not any(row.get('in_tree') for row in panel):
        with _style():
            figure=plt.figure(figsize=(12.5,5))
            _title(figure,'Phylogenetic neighborhood',sample_name or query_name)
            figure.text(.08,.55,'No named-species references are available in this saved tree.',color=MUTED)
            figure.text(.08,.45,'Detailed placement evidence remains available in the downloads.',color=MUTED)
            _save(figure,destination)
        return
    tree, nearest = prepare_tree(Phylo.read(str(newick),'newick'), query_name, maximum_references,
        reference_names=[row['tree_label'] for row in panel if row.get('in_tree')] if panel is not None else None)
    by_accession = {row['accession']:row for row in panel or []}
    leaves = tree.get_terminals()
    count = len(leaves)
    sample_name = sample_name or query_name
    depth = tree.depths()
    offset = depth[tree.root]
    depth = {clade:value-offset for clade,value in depth.items()}
    span = max(depth.values())
    scale = span if span > 0 else 1
    y = {tip:count-1-i for i,tip in enumerate(leaves)}
    for clade in tree.find_clades(order='postorder'):
        if clade.clades:
            y[clade] = (y[clade.clades[0]]+y[clade.clades[-1]])/2
    query = next(tip for tip in leaves if tip.name==query_name)
    reference = next(tip for tip in leaves if tip.name==nearest)
    ancestor = tree.common_ancestor(query,reference)
    with _style():
        figure = plt.figure(figsize=(12.5,max(6.5, .55*count+2.0)))
        axes = figure.add_axes([.04,.17,.40,.63])
        labels = figure.add_axes([.46,.17,.51,.63], sharey=axes)
        _title(figure,'Phylogenetic neighborhood',f'{sample_name}  •  {count-1} reference genomes  •  Conserved-marker placement')
        tips = ancestor.get_terminals()
        bottom,top = min(y[t] for t in tips)-.43,max(y[t] for t in tips)+.43
        axes.add_patch(Rectangle((depth[ancestor]-.008*scale,bottom),
            scale-depth[ancestor]+.04*scale,top-bottom,facecolor='#fff0f1',edgecolor='none',zorder=0))
        for tip in (query,reference):
            labels.axhspan(y[tip]-.43,y[tip]+.43,color='#fff0f1' if tip is query else '#edf4fc',zorder=0)
        for clade in tree.find_clades():
            if clade.clades:
                axes.plot([depth[clade]]*2,[y[clade.clades[0]],y[clade.clades[-1]]],color=INK,lw=1.35,zorder=2)
                for child in clade.clades:
                    axes.plot([depth[clade],depth[child]],[y[child]]*2,
                        color=RED if child is query else INK,lw=2.4 if child is query else 1.35,zorder=2)
        for tip in leaves:
            color = RED if tip is query else BLUE if tip is reference else INK
            axes.plot(depth[tip],y[tip],'o',color=color,markersize=6.4,zorder=3)
            axes.plot([depth[tip]+.015*scale,scale*1.06],[y[tip]]*2,
                      color='#d8dfe6',lw=.7,linestyle=(0,(2,3)),zorder=1)
            name,acc,cluster = tip_label(tip.name or '')
            shared = by_accession.get(acc)
            if shared:
                name = f"{shared['panel_id']}  {shared['species_name']}"
            if tip is query:
                name,acc,cluster = sample_name,'DIAGNOSTIC ISOLATE',''
            detail = acc + ('  ·  GTDB: '+cluster if cluster else '')
            if outgroup and acc == outgroup.get('reference_accession'):
                detail += '  ·  OUTGROUP'
            title_lines = textwrap.wrap(name, 52) or ['Unnamed reference']
            labels.text(.015,y[tip]+.10,'\n'.join(title_lines),color=color,
                        fontsize=11,fontstyle='normal' if tip is query else 'italic',
                        weight='bold' if tip in (query,reference) else 'normal',va='center')
            labels.text(.015,y[tip]-.23,'\n'.join(textwrap.wrap(detail,78)),color=MUTED,fontsize=7.9,va='center')
        axes.set_xlim(-.025*scale,1.08*scale)
        axes.set_ylim(-.6,count-.4)
        labels.set_xlim(0,1)
        axes.axis('off'); labels.axis('off')
        if span > 0:
            base = 10**math.floor(math.log10(span/5))
            bar = max(v*base for v in (1,2,5,10) if v*base<=span/3)
            length = .4*bar/(1.105*scale)
            figure.add_artist(Line2D([.055,.055+length],[.12,.12],transform=figure.transFigure,color=INK,lw=1.8))
            figure.text(.055,.087,f'{bar:g} substitutions per site',fontsize=8,color=MUTED)
        figure.text(.47,.115,'Isolate',color=RED,fontsize=9,weight='bold')
        figure.text(.56,.115,'Nearest reference',color=BLUE,fontsize=9,weight='bold')
        figure.text(.72,.115,'Other references',color=INK,fontsize=9)
        unmatched = [row for row in panel or [] if not row.get('in_tree')]
        note = 'Reference IDs are shared with the ANI figure, not similarity ranks. Evolutionary distance is not ANI.'
        if unmatched:
            note += '\nANI comparison only (not present in this placement): ' + ', '.join(row['panel_id'] for row in unmatched)
        else:
            note += '\nShading marks the nearest displayed reference branch. Sample bootstrap support is not inferred.'
        figure.text(.035,.033,note,fontsize=8.5,color=MUTED)
        _save(figure,destination)


def ani_rows(evidence, maximum_references=10):
    rows = []
    for value in evidence.get('reference_audit') or evidence.get('closest_species', []):
        try:
            values = [float(value[k]) for k in ('ani','aligned_fraction_query','aligned_fraction_reference')]
        except (KeyError,TypeError,ValueError):
            continue
        if any(not math.isfinite(v) or not 0<=v<=100 for v in values):
            continue
        rows.append({**value, 'ani':values[0], 'aligned_fraction_query':values[1],
            'aligned_fraction_reference':values[2], 'minimum_af':min(values[1:])})
    policy = evidence.get('identification_thresholds',{})
    ani_gate = float(policy.get('minimum_ani',95)); af_gate = float(policy.get('minimum_aligned_fraction',65))
    mode = policy.get('aligned_fraction_mode', 'both')  # Unstamped historical evidence.
    if mode not in {'either', 'both'}:
        raise ValueError('Invalid aligned-fraction mode')
    for row in rows:
        row['decision_af'] = (max(row['aligned_fraction_query'], row['aligned_fraction_reference'])
                              if mode == 'either' else row['minimum_af'])
    if not 0<=ani_gate<=100 or not 0<=af_gate<=100:
        raise ValueError('Invalid identification thresholds')
    rows.sort(key=lambda row:(row['ani']>=ani_gate and row['decision_af']>=af_gate,
                             row['decision_af'],row['ani'],row.get('accession','')),reverse=True)
    selected,seen = [],set()
    for row in rows:
        # Show separate GTDB clusters if available; never confuse another strain
        # with an additional species. No expected-organism bias in selection.
        group = row.get('gtdb_taxonomy') or row.get('species_taxid') or row.get('species_name') or row.get('accession')
        if group in seen: continue
        seen.add(group); selected.append(row)
        if len(selected)>=max(1,min(10,int(maximum_references))): break
    return selected,ani_gate,af_gate


def render_ani(evidence, destination, *, sample_name='Sample', maximum_references=10, panel=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    rows,ani_gate,af_gate = ani_rows(evidence,maximum_references)
    af_label = ('Larger' if evidence.get('identification_thresholds', {}).get('aligned_fraction_mode') == 'either'
                else 'Smaller')
    rows = [row for row in rows if named_reference(row.get('species_name'))]
    if panel is not None:
        rows = panel
    numeric = [row for row in rows if row.get('has_ani',True)]
    with _style():
        figure = plt.figure(figsize=(13.2,max(7.0,.70*len(rows)+2.5)))
        _title(figure,'Whole-genome similarity',f'{sample_name}  •  ANI and alignment evidence  •  Up to {maximum_references} representative references')
        axes = figure.add_axes([.075,.21,.44,.57])
        table = figure.add_axes([.55,.18,.43,.64]); table.axis('off')
        table.set_xlim(0,1); table.set_ylim(0,1)
        low = min(75, math.floor(min((row['ani'] for row in numeric), default=75)/5)*5)
        low = max(0,low)
        axes.add_patch(Rectangle((af_gate,ani_gate),100-af_gate,100-ani_gate,
            facecolor='#e9f4ec',edgecolor='none',zorder=0))
        axes.axvline(af_gate,color=GREEN,lw=1,ls=(0,(4,3)))
        axes.axhline(ani_gate,color=GREEN,lw=1,ls=(0,(4,3)))
        axes.text((af_gate+100)/2,(ani_gate+100)/2,'Numeric thresholds met',
            ha='center',va='center',color=GREEN,fontsize=9,weight='bold')
        axes.set_xlim(-2,102); axes.set_ylim(low-1,101)
        axes.set_xlabel(f'{af_label} of query / reference aligned fractions (%)',fontsize=9,color=INK)
        axes.set_ylabel('Average nucleotide identity (ANI, %)',fontsize=10,color=INK)
        axes.grid(color='#edf0f3',lw=.8); axes.set_axisbelow(True)
        for spine in axes.spines.values(): spine.set_color('#cbd5df')
        table.text(0,1,'REFERENCE',fontsize=8,color=MUTED,weight='bold')
        for x,label in ((.66,'ANI'),(.80,'QUERY AF'),(.96,'REF AF')):
            table.text(x,1,label,fontsize=8,color=MUTED,weight='bold',ha='right')
        annotations=[]
        highest_af = max(numeric,key=lambda row:row['decision_af']) if numeric else None
        for index,row in enumerate(rows):
            rank=row.get('panel_id',str(index+1))
            step=.85/max(1,len(rows))
            y=.91-index*step
            table.axhspan(y-step+.043,y+.043,color='#f0f3f6' if index%2==0 else 'white',zorder=0)
            table.axhline(y-step+.043,color='#dce2e8',lw=.5)
            name=row.get('species_name') or row.get('accession','Unnamed reference')
            lines=textwrap.wrap(name,30)
            table.text(.09,y+.025,'\n'.join(lines),fontsize=8.6,color=INK,va='top',linespacing=1.15)
            detail = row.get('accession','')
            if row.get('in_tree') is False:
                detail += ' · ANI only'
            if row.get('is_outgroup'):
                detail += ' · Outgroup'
            # Offset in points from the same name anchor, so each accession
            # remains attached to its own name even when that name wraps.
            table.annotate(detail,(.09,y+.025),xycoords='axes fraction',
                xytext=(0,-(10*len(lines)+3)),textcoords='offset points',
                fontsize=7.2,color=MUTED,va='top')
            if not row.get('has_ani',True):
                table.text(.01,y+.016,rank,fontsize=8,color=MUTED,weight='bold')
                table.text(.96,y+.005,'Not measurable',fontsize=8,ha='right',color=MUTED)
                continue
            meets=row['ani']>=ani_gate and row['decision_af']>=af_gate
            color=GREEN if meets else RED
            axes.scatter(row['decision_af'],row['ani'],s=110,color=color,edgecolors='white',lw=1,zorder=3)
            # Numbers link to the table without long overlapping organism labels.
            annotation=axes.annotate(str(rank),(row['decision_af'],row['ani']),
                xytext=(8,7 if index%2==0 else -12),textcoords='offset points',fontsize=8,color=INK,
                bbox=dict(facecolor='white',edgecolor='#dce2e8',boxstyle='round,pad=.15'),
                arrowprops=dict(arrowstyle='-',color=MUTED,lw=.5),zorder=4)
            annotations.append(annotation)
            table.text(.01,y+.016,str(rank),fontsize=9,color=color,weight='bold',va='center')
            for x,key,gate in ((.66,'ani',ani_gate),(.80,'aligned_fraction_query',af_gate),(.96,'aligned_fraction_reference',af_gate)):
                table.text(x,y+.005,f'{row[key]:.2f}',fontsize=9,ha='right',va='center',
                    color=GREEN if row[key]>=gate else RED,weight='bold' if row[key]<gate else 'normal')
        if not numeric:
            axes.text(.42,.35,'No reportable ANI comparisons',transform=axes.transAxes,ha='center',color=MUTED)
        if not rows:
            table.text(.01,.82,'No numeric ANI evidence was retained.\nMissing values are not zero.',color=MUTED,fontsize=10)
        # Label placement is visual only: never jitter the measured points.
        # Greedily fit numbered labels around nearby dots in display coordinates.
        figure.canvas.draw()
        renderer=figure.canvas.get_renderer()
        occupied=[]
        from matplotlib.transforms import Bbox
        point_boxes=[]
        for row in numeric:
            px,py=axes.transData.transform((row['decision_af'],row['ani']))
            point_boxes.append(Bbox.from_bounds(px-7,py-7,14,14))
        for annotation in annotations:
            for dx,dy in [(x,y) for x in (10,-20,25,-40,42,-60,60,-80,80) for y in (10,-12,24,-25,38,-38,54,-54)]:
                annotation.set_position((dx,dy))
                annotation.update_positions(renderer)
                annotation.update_bbox_position_size(renderer)
                box=annotation.get_bbox_patch().get_window_extent(renderer).expanded(1.25,1.25)
                if (axes.bbox.contains(box.x0,box.y0) and axes.bbox.contains(box.x1,box.y1)
                    and not any(box.overlaps(other) for other in occupied+point_boxes)):
                    occupied.append(box)
                    break
        figure.text(.055,.11,f'Guides: ANI ≥ {ani_gate:g}%  ·  Query AF ≥ {af_gate:g}%  ·  Reference AF ≥ {af_gate:g}%',fontsize=10,color=INK)
        note = 'Reference IDs are not ranks. Low-AF comparisons represent small shared regions, not species identification.'
        if highest_af:
            af_summary = 'maximum' if af_label == 'Larger' else 'minimum'
            note += f"\nHighest {af_summary} aligned fraction: {highest_af.get('panel_id',highest_af.get('accession',''))} ({highest_af['decision_af']:.2f}%). Missing estimates are not zero."
        figure.text(.055,.052,note,fontsize=8.5,color=MUTED)
        _save(figure,destination)
