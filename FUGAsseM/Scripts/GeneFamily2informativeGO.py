import sys
import re
import argparse
import pandas as pd
from collections import Counter
from fugassem.common.utils import warn, qw
from fugassem.common.dictation import polymap, col2dict

sys.setrecursionlimit(50000)

# Constants
ALLOW_CROSSTALK = False
c_prop_pattern = r"^([^\s]*?): (.*)"
c_goid_pattern = r"(GO:\d{7})"

c_namespace_convert = qw("""
biological_process BP
molecular_function MF
cellular_component CC
""", as_dict=True)

c_relationship_subtypes = qw("""
part_of
""")

c_singular_props = qw("""
id
name
namespace
""")

parentage_types = Counter()

# Stanza class
class Stanza:
    def __init__(self):
        self.lines = []
        self.props = {}

    def __getitem__(self, key):
        if key in c_singular_props:
            assert len(self.props[key]) == 1, f"not a singular relationship: {key} {self.props[key]}"
            return self.props[key][0]
        return self.props[key]

    def __contains__(self, key):
        return key in self.props

    def get(self, key, default=None):
        return self[key] if key in self else default

    def populate(self):
        for line in self.lines:
            match = re.search(c_prop_pattern, line)
            if match:
                propname, line = match.groups()
                self.props.setdefault(propname, []).append(line)

    def get_parent_goids(self):
        parent_goids = []
        for line in self.props.get("is_a", []):
            parent_goids.append(line.split()[0])
            parentage_types["is_a"] += 1
        for line in self.props.get("relationship", []):
            items = line.split()
            if items[0] in c_relationship_subtypes:
                parent_goids.append(items[1])
                parentage_types["relationship:" + items[0]] += 1
        for line in self.props.get("intersection_of", []):
            items = line.split()
            for k in [0, 1]:
                if "GO:" in items[k]:
                    parent_goids.append(items[k])
                    parentage_types["intersection_of"] += 1
        return [k for k in parent_goids if re.search(c_goid_pattern, k)]

# OBO Parser
class OBOParser:
    def __init__(self, p_obo):
        self.stanzas = []
        IN_TERMS = False
        with open(p_obo) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    self.stanzas.append(Stanza())
                    IN_TERMS = True
                elif IN_TERMS and line == "[Typedef]":
                    IN_TERMS = False
                elif IN_TERMS and line != "":
                    self.stanzas[-1].lines.append(line)
        for s in self.stanzas:
            s.populate()

    def extract_terms(self):
        for stanza in self.stanzas:
            yield Term(stanza)

# Term class
class Term:
    def __init__(self, stanza):
        self.goid = stanza["id"]
        self.name = stanza["name"]
        self.namespace = stanza["namespace"]
        self.namespace_short = c_namespace_convert[self.namespace]
        self.is_obsolete = "is_obsolete" in stanza
        self.replaced_by = stanza.get("replaced_by", None)
        self.alt_ids = stanza.get("alt_id", [])
        self.parent_ids = stanza.get_parent_goids()
        self.parents = set()
        self.children = set()
        self.genes = set()
        self.depth = None
        self.is_root = False
        self.is_leaf = False
        self.progeny = None
        self.progeny_genes = None
        self.is_informative = False
        self.is_pruned = False
        self.is_acceptable = True

    def __repr__(self):
        return f"{self.goid}: [{self.namespace_short}] {self.name}"

    def add_gene(self, gene):
        self.genes.add(gene)

    def get_progeny(self):
        if self.progeny is None:
            self.progeny = {self}
            for child in self.children:
                self.progeny.update(child.get_progeny())
        return self.progeny

    def get_progeny_genes(self):
        if self.progeny_genes is None:
            self.progeny_genes = set(self.genes)
            for child in self.children:
                self.progeny_genes.update(child.get_progeny_genes())
        return self.progeny_genes

# Ontology class
class Ontology:
    def __init__(self, p_obo):
        self.terms = {}
        self.idmap = {}
        self.roots = []
        self.leaves = []
        self.attached_genes = set()

        for term in OBOParser(p_obo).extract_terms():
            if not re.search(c_goid_pattern, term.goid):
                continue
            if term.is_obsolete:
                if term.replaced_by:
                    self.idmap[term.goid] = term.replaced_by
            else:
                self.terms[term.goid] = term
                self.idmap[term.goid] = term.goid
                for alt_id in term.alt_ids:
                    self.idmap[alt_id] = term.goid

        for cterm in self.terms.values():
            for parent_id in cterm.parent_ids:
                pterm = self.terms.get(parent_id)
                if not pterm:
                    continue
                if ALLOW_CROSSTALK or pterm.namespace == cterm.namespace:
                    cterm.parents.add(pterm)
                    pterm.children.add(cterm)
                else:
                    warn("Cross-ontology relationship ignored:", pterm, cterm)
                    parentage_types["Cross-ontology [ignored]"] += 1

        for term in self.terms.values():
            if not term.parents:
                term.is_root = True
                self.roots.append(term)
            if not term.children:
                term.is_leaf = True
                self.leaves.append(term)

        def recurse_depth(term):
            for child in term.children:
                d = term.depth + 1
                if child.depth is None or child.depth > d:
                    child.depth = d
                    recurse_depth(child)

        for root in self.roots:
            root.depth = 0
            recurse_depth(root)

    def iter_terms(self):
        return self.terms.values()

    def attach_genes(self, mapping):
        for gene, goids in mapping.items():
            self.attached_genes.add(gene)

            # Fix nested list problem
            flat_goids = []
            for g in goids:
                if isinstance(g, list):
                    flat_goids.extend(g)
                else:
                    flat_goids.append(g)

            for goid in flat_goids:
                goid = self.idmap.get(goid, goid)
                if isinstance(goid, list):
                    goid = goid[0]
                if goid in self.terms:
                    self.terms[goid].add_gene(gene)

    def prune(self, goid):
        def recurse(term):
            term.is_pruned = True
            for c in term.children:
                if not c.is_pruned:
                    recurse(c)
        recurse(self.terms[goid])

    def set_informative(self, threshold):
        def recurse(term):
            if len(term.get_progeny_genes()) >= threshold:
                term.is_informative = True
                for child in term.children:
                    if len(child.get_progeny_genes()) >= threshold:
                        term.is_informative = False
                        break
            else:
                for parent in term.parents:
                    recurse(parent)
        for leaf in self.leaves:
            recurse(leaf)

# Argument parsing
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("obo", help="GO-provided obo file")
    parser.add_argument("--mapping", help="gene↔GO mapping file")
    parser.add_argument("--allowed-genes", help="optional list of allowed genes")
    parser.add_argument("--flip", action="store_true", help="flip gene/GO mapping direction")
    parser.add_argument("--depth", type=int, help="only include terms at this depth")
    parser.add_argument("--grep", help="filter terms by name regex")
    parser.add_argument("--prune", help="only include descendants of this GO term")
    parser.add_argument("--namespace", nargs="+", choices=["BP", "MF", "CC"], help="only include terms from these namespaces")
    parser.add_argument("--informative", help="min # or fraction of genes to be informative")
    parser.add_argument("--ignore-progeny", action="store_true", help="do not inherit annotations through DAG")
    parser.add_argument("--terms-only", action="store_true", help="only output GO terms")
    parser.add_argument("--outfile", default=None, help="output file")
    return parser.parse_args()

# Main
def main():
    args = get_args()
    obo = Ontology(args.obo)

    warn("Summary of relationship types:")
    for k in sorted(parentage_types):
        warn(k, parentage_types[k])

    mapping = {}
    if args.mapping:
        mapping = polymap(args.mapping, reverse=args.flip)

        # Flatten nested lists
        for k, v in mapping.items():
            flat = []
            for item in v:
                if isinstance(item, list):
                    flat.extend(item)
                else:
                    flat.append(item)
            mapping[k] = flat

        if args.allowed_genes:
            allowed = col2dict(args.allowed_genes)
            mapping = {k: v for k, v in mapping.items() if k in allowed}

        obo.attach_genes(mapping)
        warn("# of attached genes:", len(obo.attached_genes))

    if args.informative:
        threshold = float(args.informative)
        if threshold < 1:
            threshold *= len(obo.attached_genes)
        obo.set_informative(int(threshold))
        for term in obo.iter_terms():
            if not term.is_informative:
                term.is_acceptable = False

    if args.prune:
        obo.prune(args.prune)
        for term in obo.iter_terms():
            if not term.is_pruned:
                term.is_acceptable = False

    if args.depth is not None:
        for term in obo.iter_terms():
            if term.depth != args.depth:
                term.is_acceptable = False

    if args.grep:
        for term in obo.iter_terms():
            if not re.search(args.grep, term.name):
                term.is_acceptable = False

    if args.namespace:
        for term in obo.iter_terms():
            if term.namespace_short not in args.namespace:
                term.is_acceptable = False

    # Output
    with open(args.outfile, "w") if args.outfile else sys.stdout as fh:
        for term in obo.iter_terms():
            if term.is_acceptable:
                outline = [str(term)]
                if not args.terms_only:
                    genes = term.genes if args.ignore_progeny else term.get_progeny_genes()
                    outline += sorted(genes)
                fh.write("\t".join(outline) + "\n")

if __name__ == "__main__":
    main()
