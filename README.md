
# AMR Metagenomics Pipeline

A Snakemake-based pipeline for identifying antimicrobial resistance (AMR) genes in metagenomic sequencing data and investigating their potential bacterial taxonomic origins.

The pipeline combines multiple AMR detection approaches:

- **StarAMR/ResFinder:** Detection of resistance genes in assembled sequences.
- **CARD/RGI (contig):** Detection of resistance-associated genes in assembled sequences.
- **CARD/RGI (reads):** Read-based resistance analysis using KMA.
- **argNorm:** Standardization of resistance annotations using CARD ARO identifiers and ontology information.
- **hAMRonization:** Conversion of AMR detection results into standardized formats.
- **Kraken2 and BLAST:** Taxonomic investigation of resistance-associated scaffolds selected using StarAMR/ResFinder.

The workflow also generates per-sample Excel reports comparing AMR detections across the different analysis methods and preserving their supporting evidence.


## Requirements

The pipeline runs on Linux (also tested under WSL2) and requires Conda and Snakemake.

Four Conda environments are used:

- `ar-pipeline` — Main workflow and analysis tools
- `rgi` — CARD/RGI analysis
- `argnorm` — AMR annotation normalization
- `hamronization` — Standardized AMR reporting

Environment definitions are provided in `envs/`.

Kraken2, BLAST and CARD reference databases must be installed separately.

## Installation

Clone the repository and enter the project directory.

Create the required Conda environments:

```bash
conda env create -f envs/ar-pipeline.yaml
conda env create -f envs/rgi.yaml
conda env create -f envs/argnorm.yaml
conda env create -f envs/hamronization.yaml
```

Activate the main environment before running the pipeline:

```bash
conda activate ar-pipeline
```

The workflow also requires separately downloaded Kraken2, BLAST and CARD databases. Their locations must be specified in `config.yaml`.


## Input data

The pipeline requires paired-end Illumina FASTQ files in the `inputs/` directory.

Files must follow this naming convention:

- `SampleName_1.fastq.gz`
- `SampleName_2.fastq.gz`

Create a `samples.csv` file using `samples.example.csv` as a template. List one sample identifier per line, without the FASTQ suffixes.

For example, `ExampleSample_S1_L001` corresponds to:

- `inputs/ExampleSample_S1_L001_1.fastq.gz`
- `inputs/ExampleSample_S1_L001_2.fastq.gz`

Only samples listed in `samples.csv` will be processed.

## Configuration

Pipeline settings are defined in `config.yaml`.

Before running the workflow, update the following settings:

- `kraken2.db` — Path to the Kraken2 database.
- `blast.db` — Path to the BLAST nucleotide database.
- `card.database_path` — Path to the CARD database.
- `threads` settings — Number of CPU threads assigned to individual tools (configured separately under `metaspades`, `staramr`, `kraken2`, `blast` and `card`).

The configuration also defines input/output directories, the sample list and software/database version metadata.

Database paths are machine-specific and must be adjusted for each installation.

## Usage

Run the pipeline from the project directory after preparing the input files, sample list, configuration and Conda environments.

Activate the main environment:

```bash
conda activate ar-pipeline
```

Preview the workflow without executing any analyses:

```bash
snakemake --cores 6 --dry-run
```

Run the complete pipeline:

```bash
snakemake --cores 6
```

The `--cores` parameter specifies the total number of CPU cores available to Snakemake. Adjust it according to the available hardware.

The workflow processes the samples listed in `samples.csv` and writes results to the output directory specified in `config.yaml` (default: `amr/`).

## Output files

Results are saved in the directory specified by `output_folder` in `config.yaml` (default: `amr/`).

The main outputs include:

- **metaSPAdes:** Assembled scaffolds (`{sample}_metaspades/scaffolds.fasta`).
- **StarAMR:** Resistance gene detections (`{sample}_staramr/`).
- **CARD/RGI:** Contig- and read-based AMR detection results (`card/`).
- **Kraken2 and BLAST:** Taxonomic analysis of resistance-associated scaffolds.
- **hAMRonization:** Standardized AMR detection tables and summary reports (`hamronized/`).
- **argNorm:** Normalized resistance annotations (`argnorm/`).
- **Cross-tool comparison:** Per-sample Excel reports (`comparison/reports/`).

Each sample's comparison folder contains:

- `AMR_cross_tool_comparison.xlsx` — Comparison of ARO detections across StarAMR, RGI contigs and RGI reads.
- `AMR_detection_evidence.xlsx` — Detailed detection evidence, including original tool-specific information.

Sample report folder names omit Illumina sample and lane suffixes (e.g., ExampleSample_S1_L001 becomes ExampleSample).
