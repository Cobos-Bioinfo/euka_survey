# Eukaryote Survey

## Description
This project provides an automated data aggregation pipeline and exploration toolkit that catalogues genomic sequencing data across the entire Eukaryotic tree of life. It tracks whole genome assemblies, functional annotations, and RNA-seq reads (divided into short- and long-reads) to build a fast, queryable local SQLite database.

**Main Use Case:** Identifying clades, families, or species that lack specific types of sequencing data, or discovering clades rich in genomic resources for comparative studies. 

**Target Users:** Bioinformaticians, evolutionary biologists, and comparative genomicists conducting broad taxonomic surveys or meta-analyses who need to evaluate available molecular resources without manually navigating NCBI or ENA portals.

## Features
- **Data Aggregation**: Automatically fetches data from NCBI Datasets, Annotrieve API, and EBI ENA API.
- **High-Performance Taxonomy Traversal**: Uses a highly optimized SQLite CTE query to bypass standard Python-level tree traversals, ensuring massive trees (like all Eukaryotes) are expanded rapidly.
- **Local Relational Database**: Consolidates gathered metrics into a portable SQLite database for rapid downstream querying and summation.
- **Detailed TSV Exports**: Dumps organism-by-organism resource metrics into TSVs ready for graphing or further downstream analysis.

## Project Structure
- `pipeline_build_db.py`: The main execution script. Orchestrates data collection across modules and triggers database construction.
- `query_clade.py`: The primary CLI reporting tool to explore and summarize genomic data availability for specific taxonomic clades.
- `get_taxa_by_rank.py`: A helper CLI tool to extract all descendent taxa at a specified taxonomic rank (e.g., phylum, class) under a given root.
- `ete_utils.py`: The optimized routing and taxonomy traversal backbone using recursive SQL Common Table Expressions (CTE).
- `build_db`: Directory containing specific modular fetchers:
  - `get_assemblies.py`: Retrieves sequenced genome assembly stats using the NCBI Datasets CLI.
  - `get_annotations.py`: Retrieves functional annotation frequencies from the Annotrieve API.
  - `get_reads.py`: Fetches long and short RNA-seq read runs from the EBI ENA portal.
  - `build_database.py`: Handles creation, population, and `INSERT OR REPLACE` logic for the SQLite database.

## Installation
1. **Clone the repository:**
   ```bash
   git clone https://github.com/Cobos-Bioinfo/euka_survey.git
   cd euka_survey
   ```
2. **Install Python dependencies:**
   Use the provided environment.yml file to create a conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate euka_tracker
   ```
3. **Install External Requirements:**
   Installing the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) tool in your PATH is required for gathering assembly data.

### Docker Option:
A Docker image option is in development. For users preferring containerization, a Dockerfile will be provided to encapsulate all dependencies and tools, allowing execution without local environment setup.


## Usage

### 1. Build the Database
Construct the local SQLite database containing organism metrics. (Note: Currently hardcoded for all Eukaryotes, taxID 2759):
```bash
python pipeline_build_db.py
```
*Output: `eukaryote_taxid_features_YYYY_MM_DD.db`*

### 2. Discover Clades by Rank
Find taxonomic IDs for specific lineages (e.g., finding the taxIDs for all eukaryotic phyla):
```bash
python get_taxa_by_rank.py <root_taxid> <rank>
```

### 3. Query Target Clades
Summarize absolute feature counts and unique organism breakdowns for lineages of interest. Can be done natively, via piped input, or from a tracking file:
```bash
# Query an individual clade
python query_clade.py <taxid> --db eukaryote_taxid_features_YYYY_MM_DD.db

# Output detailed per-species statistics to TSV
python query_clade.py <taxid> --db eukaryote_taxid_features_YYYY_MM_DD.db --tsv ./results_directory/
```

## Workflow Overview
1. **Initialization:** User ensures dependencies and NCBI Datasets CLI tool are correctly configured in their environment.
2. **Data Aggregation:** `pipeline_build_db.py` identifies relevant descendants via `ete_utils.py`. The pipeline successively queries external databases (NCBI, Annotrieve, ENA) via the specialized submodules to retrieve count dictionaries in memory. Finally, `build_database.py` generates the cohesive, localized SQLite file.
3. **Exploration & Discovery:** User explores specific ranks using `get_taxa_by_rank.py` to isolate taxids of interest.
4. **Summary & Evaluation:** Users execute `query_clade.py` against the constructed database, yielding high-level statistical summaries or detailed per-species TSV dumps for downstream evaluation.

## Notes / Limitations
- **Exclusion of Human/Mouse data**: RNA-seq runs for humans (taxID 9606) and mice (taxID 10090) are explicitly hardcoded to be excluded from ENA queries. This is an intentional project design to avoid significant API bloat and delays for highly sequenced model organisms.
- **Hardcoded Root**: The root database creation script (pipeline_build_db.py) relies on a hardcoded top-level taxonomic target (Eukaryota; 2759).

## Future Improvements
- Implement a new CLI flag to include human and mice RNA-seq data for users interested in those taxa, while maintaining the default exclusion for general surveys.
- Add a Dockerfile and pre-built image for users who prefer containerized environments.
- Expand the database schema to include additional metadata fields (e.g., assembly quality metrics)
