"""
Fetches descendant taxonomic IDs for a given parent taxonomic ID.
Utilizes the local NCBI taxonomy database via the ete3 library.
"""

from ete3 import NCBITaxa

NCBI = NCBITaxa()

def get_descendant_taxids(taxid: int) -> list[int]:
    """Get all descendant taxonomic IDs (including the input taxid and intermediate nodes)."""
    tree_full = NCBI.get_descendant_taxa(parent = taxid, intermediate_nodes=True)
    return tree_full + [taxid]

def get_descendant_organisms_taxids(taxid: int) -> list[int]:
    """Get all descendant organism taxonomic IDs (not including the input taxid and intermediate nodes)."""
    tree = NCBI.get_descendant_taxa(parent = taxid)
    return tree

def get_name_from_taxid(taxid: int) -> str:
    """Get the scientific name for a given taxonomic ID."""
    names = NCBI.get_taxid_translator([taxid])
    return names.get(taxid, "Unknown")

def get_rank_from_taxid(taxid: int) -> str:
    """Get the taxonomic rank for a given taxonomic ID."""
    ranks = NCBI.get_rank([taxid])
    return ranks.get(taxid, "Unknown")