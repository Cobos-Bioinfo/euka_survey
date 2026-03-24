# Eukaryote Survey

A set of tools to track and build a database of eukaryotic taxonomic IDs and their associated genomic features (assemblies, annotations, and RNA-seq reads). 

## Project Structure

- `build_db.py` - The main orchestrator script that runs the pipeline to fetch data and construct the database.
- `scripts/build_db/` - Subdirectory containing the modules used to fetch data from different sources:
  - `get_taxids.py`: Retrieves descendant taxonomic IDs from NCBI using `ete3`.
  - `get_assemblies.py`: Fetches assembly information via the NCBI datasets CLI.
  - `get_annotations.py`: Retrieves annotation data from Annotrieve.
  - `get_reads.py`: Fetches short and long RNA-seq reads from the ENA portal API.
  - `build_database.py`: Handles the SQLite database creation and data insertion.

## Prerequisites

1. **Conda Environment**: Create and activate the Conda environment using the provided `environment.yml` file:
   ```bash
   conda env create -f environment.yml
   conda activate euk_survey
   ```

2. **NCBI Datasets CLI**: The `get_assemblies.py` script requires the `datasets` command-line tool.
   Install instructions: [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)

## Usage

### 1. Build the Database
To fetch current data and build the local SQLite database, run:
```bash
python pipeline_build_db.py
```
This will create a dated SQLite database file (e.g., `eukaryote_taxid_features_YYYY_MM_DD.db`) in the root directory. Expected runtime depends on network delays, but usually completes in ~2-3 minutes.

## Query the Database
Use `get_taxa_by_rank.py` to obtain a tsv file of taxonomic IDs - names for a specified taxonomic rank (e.g., "family", "genus", "species"):
```bash
# Print to terminal
python get_taxa_by_rank.py 2759 phylum

# Save to file (tab-separated: taxid + name)
python get_taxa_by_rank.py 2759 phylum --out eukaryote_phyla.txt

# Feed the taxIDs directly into the CLI (cuts the first column)
python get_taxa_by_rank.py 2759 phylum --out eukaryote_phyla.txt
cut -f1 eukaryote_phyla.txt | xargs python query_clade.py --db eukaryote_taxid_features_2026_03_19.db
```
The output file is tab-separated (taxid + name) rather than taxid-only so it's human-readable, but `query_clade.py --file` only needs the taxid column — which is why the `cut -f1` above is useful if you want to pipe them directly.
A quick note on valid rank strings for ete3: `superkingdom`, `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species` are the standard ones. The script will just return nothing if you typo the rank, so it prints a clear error in that case.

# Get all species AND subspecies
## Retrieving All Species and Subspecies from NCBI Taxonomy

### Objective
Retrieve all taxonomic IDs (taxIDs) belonging to species **and** subspecies under a given parent taxon using the `ete3` library and its NCBI taxonomy database.

### Problem
`ete3`'s `get_descendant_taxa()` behaves as follows depending on the `collapse_subspecies` flag:
- `collapse_subspecies=True` — returns **species only**, subspecies are collapsed into their parent species.
- `collapse_subspecies=False` — returns **subspecies only**, parent species that have subspecies under them are excluded.

The issue is that collapse_subspecies=True triggers an extremely expensive operation when called on a large clade like Eukaryotes (taxid 2759).


What's happening: When collapse_subspecies=True, ete3's get_descendant_taxa() doesn't just return taxids — it:

Traverses the entire tree to find all descendants
For each taxid, queries the database to check its rank
Filters out subspecies/varietas/formas by comparing ranks
Validates parent-child relationships in the lineage

For Eukaryotes, you're dealing with:

~2.7 million taxids with collapse=False
Potentially millions of database lookups with collapse=True

### Solution
Bypass `get_descendant_taxa()` entirely and query the ete3 SQLite database directly using a **single recursive CTE**. This walks the full subtree in one pass inside the database engine, streaming rows out one at a time and keeping memory usage constant regardless of tree size. The query filters for all target ranks — `species`, `subspecies`, `varietas`, `forma`, and `strain` — in a single call, replicating the union result without the memory overhead.
