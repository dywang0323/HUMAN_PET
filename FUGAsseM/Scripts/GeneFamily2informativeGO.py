#!/usr/bin/env python3
import sys
import re
import argparse
import pandas as pd
from collections import Counter

sys.setrecursionlimit(50000)

# Patterns and constants
c_prop_pattern = r"^([^\s]*?): (.*)"
c_goid_pattern = r"(GO:\d{7})"
ALLOW_CROSSTALK = False

c_namespace_convert = {
    "biological_process": "BP",
    "molecular_function": "MF",
    "cellular_component": "CC"
}

c_relationship_subtypes = {"part_of"}
c_singular_props = {"id", "name", "namespace"}
parentage_types = Counter()


def warn(*args):
    print("WARNING (GENEFAMILY2INFORMATIVEGO.PY):", *args, file=sys.stderr)


# =========================
# GO OBO Parsing Classes
# =========================

class Stanza:
    def __init__(self):
        self.lines = []
        self.props = {}

    def __getitem__(self, key):
        if key in c_singular_props:
            assert len(self.props[key]) == 1
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
        return [k for k in parent_goids if re.search(c_goid_pattern, k)]


class OBOParser:
    def __init__(self, path):
        self.stanzas = []
        IN_TERMS = False
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    self.stanzas.append(Stanza())
                    IN_TERMS = True
                elif IN_TERMS and line == "[Typedef]":
                    IN_TERMS = False
                elif IN_TERMS and line:
                    self.stanzas[-1].lines.append(line)
        for s in self.stanzas:
            s.populate()

    def extract_terms(self):
        for stanza in self.stanzas:
            yield Term(stanza)


class Term:
    def __init__(self, stanza):
        self.goid = stanza["id"]
        self.name = stanza["name"]
        self.namespace = stanza["namespace"]
        self.namespace_short = c_namespace_convert.get(self.namespace, "NA")
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
            for cterm in self.children:
                self.progeny.update(cterm.get_progeny())
        return self.progeny

    def get_progeny_genes(self):
        if self.progeny_genes is None:
            self.progeny_genes = set(self.genes)
            for cterm in self.children:
                self.progeny_genes.update(cterm.get_progeny_genes())
        return self.progeny_genes


class Ontology:
    def __init__(self, obo_path):
        self.terms = {}
        self.idmap = {}
        self.roots = []
        self.leaves = []
        self.attached_genes = set()

        for term in OBOParser(obo_path).extract_terms():
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

        for term in self.terms.values():
            for parent_id in term.parent_ids:
                parent = self.terms.get(parent_id)
                if not parent:
                    continue
                if ALLOW_CROSSTALK or parent.namespace == term.namespace:
                    term.parents.add(parent)
                    parent.children.add(term)
                else:
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
            for goid in (goids if isinstance(goids, list) else [goids]):
                if isinstance(goid, dict):
                    continue
                goid = str(goid)
                goid = self.idmap.get(goid, goid)
                if goid in self.terms:
                    self.terms[goid].add_gene(gene)

    def set_informative(self, threshold):
        def recurse(term):
            if len(term.get_progeny_genes()) >= threshold:
                term.is_informative = True
                for c in term.children:
                    if len(c.get_progeny_genes()) >= threshold:
                        term.is_informative = False
                        break
            else:
                for p in term.parents:
                    recurse(p)

        for leaf in self.leaves:
            recurse(leaf)


# ============ MAIN ============

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("obo", help="GO .obo file")
    parser.add_argument("--mapping", required=True, help="Gene to GO mapping file (tab-delimited)")
    parser.add_argument("--informative", type=float, required=True, help="Minimum number of genes per term")
    parser.add_argument("--namespace", choices=["BP", "MF", "CC"], help="Namespace filter")
    parser.add_argument("--terms-only", action="store_true", help="Only output informative terms")
    parser.add_argument("--outfile", required=True, help="Output file path")
    return parser.parse_args()


def main():
    args = parse_args()

    # Load gene-to-GO mapping
    mapping = {}
    df = pd.read_csv(args.mapping, sep="\t", header=None, names=["Gene", "GO"])
    for _, row in df.iterrows():
        mapping.setdefault(row["Gene"], []).append(row["GO"])

    # Load GO Ontology
    obo = Ontology(args.obo)
    obo.attach_genes(mapping)
    warn("# of attached genes:", len(obo.attached_genes))

    # Set informative
    threshold = int(args.informative)
    obo.set_informative(threshold)

    # Apply namespace filter
    terms = []
    for term in obo.iter_terms():
        if not term.is_informative:
            continue
        if args.namespace and term.namespace_short != args.namespace:
            continue
        terms.append(term)

    warn(f"# of informative terms: {len(terms)}")

    with open(args.outfile, "w") as out:
        for term in sorted(terms, key=lambda x: x.goid):
            out.write(str(term) + "\n")


if __name__ == "__main__":
    main()
