#!/usr/bin/env python3
"""
get_taxa_by_rank.py — Get all taxIDs of a given rank under a clade.

Usage:
    python get_taxa_by_rank.py 2759 phylum
    python get_taxa_by_rank.py 2759 phylum --out eukaryote_phyla.txt
    python get_taxa_by_rank.py 2759 class --out eukaryote_classes.txt
"""

import argparse
import sys
from ete3 import NCBITaxa

NCBI = NCBITaxa()


def get_taxa_at_rank(root_taxid: int, rank: str) -> list[tuple[int, str]]:
    """Return all (taxid, name) pairs at the given rank under root_taxid."""
    descendants = NCBI.get_descendant_taxa(root_taxid, intermediate_nodes=True)
    ranks = NCBI.get_rank(descendants)
    hits = [taxid for taxid, r in ranks.items() if r == rank]
    names = NCBI.get_taxid_translator(hits)
    return sorted(names.items(), key=lambda x: x[1])  # sort by name


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("taxid", type=int, help="Root clade taxID")
    parser.add_argument("rank", type=str, help="Rank to filter (e.g. phylum, class, order)")
    parser.add_argument("--out", "-o", metavar="FILE",
                        help="Output file (default: print to stdout)")
    args = parser.parse_args()

    results = get_taxa_at_rank(args.taxid, args.rank)

    if not results:
        print(f"No taxa found at rank '{args.rank}' under taxID {args.taxid}.", file=sys.stderr)
        sys.exit(1)

    lines = [f"{taxid}\t{name}" for taxid, name in results]

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"Wrote {len(results)} taxa to {args.out}")
    else:
        print("\n".join(lines))


if __name__ == "__main__":
    main()
