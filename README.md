# TEA - Transposable Element Analysis Pipeline

A comprehensive Snakemake pipeline for transposable element (TE) differential abundance analysis using STAR, RepeatMasker, and TEtranscripts.

## Overview

TEA (Transposable Element Analysis) processes RNA-seq data to identify differentially expressed transposable elements between experimental conditions. The pipeline uses a sequential workflow that combines multiple industry-standard tools:

1. Automatic download of reference genomes and annotations from Gencode
2. STAR index creation and RNA-seq alignment with optimized TE detection parameters
3. RepeatMasker annotation of transposable elements
4. TEtranscripts for differential expression analysis between treatment and control groups
5. Quality control with FastQC and MultiQC

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
  - [Sample Configuration](#sample-configuration)
  - [Reference Genome Configuration](#reference-genome-configuration)
  - [Computational Resources](#computational-resources)
- [Running TEA](#running-tea)
  - [Basic Usage](#basic-usage)
  - [Dry Run](#dry-run)
  - [Advanced Options](#advanced-options)
- [Pipeline Steps](#pipeline-steps)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [License](#license)

## Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://github.com/mamba-org/mamba)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- Tools (automatically installed with the environment.yaml file):
  - STAR aligner
  - RepeatMasker
  - TEtranscripts
  - FastQC
  - MultiQC
  - Samtools
  - Python 3.9+

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Leran10/TEA.git
   cd TEA
   ```

2. Create a minimal conda environment (to avoid disk quota issues):
   ```bash
   # Option 1: Create the full environment (may require significant disk space)
   conda env create -f environment.yaml
   
   # Option 2: Create a minimal environment on systems with limited space
   conda create -n te_minimal python=3.9 snakemake star repeatmasker fastqc samtools pip
   conda activate te_minimal
   pip install TEtranscripts==2.2.3
   ```

3. Activate the environment:
   ```bash
   # If using the full environment
   conda activate te_analysis
   
   # Or if using the minimal environment
   conda activate te_minimal
   ```

4. Verify the installation:
   ```bash
   # Check that TEtranscripts was installed
   TEtranscripts --version
   
   # Check that all other tools are available
   STAR --version
   RepeatMasker -h | head -n 5
   fastqc --version
   ```

## Quick Start

1. Edit `config.yaml` to specify your samples and reference genome (see [Configuration](#configuration))
2. Place your paired-end FASTQ files in `data/raw_fastq/` with naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
3. Run the pipeline:
   ```bash
   ./TEA -c 16
   ```

## Configuration

All pipeline settings are managed in the `config.yaml` file. Here are the key parameters you need to configure:

### Sample Configuration

```yaml
# Samples and conditions
samples:
  - "sample1"
  - "sample2"
  - "sample3"
  - "sample4"
  - "sample5"
  - "sample6"

treatment_samples:
  - "sample2"
  - "sample3"
  - "sample5"

control_samples:
  - "sample1"
  - "sample4"
  - "sample6"
```

- `samples`: List only the sample base names (without the _R1.fastq.gz or _R2.fastq.gz suffixes)
- `treatment_samples`: List the sample IDs in your treatment group
- `control_samples`: List the sample IDs in your control group

For each sample in your list, you must have two paired-end FASTQ files in the `data/raw_fastq/` directory following this naming convention:
- `{sample}_R1.fastq.gz` - First read in the pair
- `{sample}_R2.fastq.gz` - Second read in the pair

For example, if you list `tumor1` in your samples, you should have these files:
- `data/raw_fastq/tumor1_R1.fastq.gz`
- `data/raw_fastq/tumor1_R2.fastq.gz`

Additionally, create a metadata file at `metadata/sample_info.txt` with the following format:
```
SampleID    Condition    Batch
sample1     control      1
sample2     treatment    1
sample3     treatment    1
sample4     control      1
sample5     treatment    1
sample6     control      1
```

### Reference Genome Configuration

```yaml
# Reference files
# Specify paths to your reference genome and annotation files
genome_fasta: "/path/to/your/genome.fa"
gtf_file: "/path/to/your/annotation.gtf"
te_gtf_file: "/path/to/your/te_annotation.gtf"  # Path to your pre-existing TE GTF file
```

- `genome_fasta`: Full path to your genome FASTA file
- `gtf_file`: Full path to your gene annotation GTF file
- `te_gtf_file`: Full path to your TE annotation GTF file (pre-existing)

> **Note:** You must download and prepare these reference files yourself. The files can be located anywhere on your system - just provide the full absolute paths in the config file.
> 
> Good sources for reference files include:
> - [Gencode](https://www.gencodegenes.org/): High-quality annotations for human and mouse
> - [Ensembl](https://www.ensembl.org/): Reference genomes and annotations for multiple species
> - [UCSC Genome Browser](https://genome.ucsc.edu/): Various genome builds and annotations

### Computational Resources

```yaml
# Computational resources
threads:
  star_index: 12
  star_align: 8
  repeatmasker: 8
  samtools: 8
  tetx: 8
```

Adjust threads based on your available computing resources.

## Running TEA

### Basic Usage

After configuring `config.yaml` and placing your FASTQ files in the correct location:

```bash
./TEA -c <number_of_cores>
```

### Command-line Options

The `TEA` command provides a convenient wrapper around Snakemake with the following options:

```
Usage:
  TEA [options] [snakemake_args]

Options:
  -h, --help       Show help message and exit
  -n, --dry-run    Perform a dry run (no actual execution)
  -c, --cores N    Specify number of cores to use (default: 1)
  --profile NAME   Use the specified Snakemake profile
  --configfile FILE Use the specified config file instead of default
```

### Dry Run

To check the workflow without executing any commands:

```bash
./TEA -n
```

### Advanced Usage

Run a specific part of the pipeline:

```bash
# Run only up to the STAR alignment
./TEA -c 16 results/star/{sample}/Aligned.sortedByCoord.out.bam
```

Run with a cluster/HPC system:

```bash
./TEA --profile slurm
```

You can pass any additional Snakemake arguments after the TEA options. For more advanced Snakemake options, see the [Snakemake documentation](https://snakemake.readthedocs.io/).

## Pipeline Steps

The pipeline executes in the following sequence:

1. **Reference File Requirements**
   - You must provide your own genome FASTA, gene annotation GTF, and TE annotation GTF files
   - Specify their paths in the config.yaml file

2. **STAR Indexing**
   - Create STAR index with your provided reference genome and GTF files (using the paths from config.yaml)

3. **RNA-seq Alignment with STAR**
   - Align reads with optimized parameters for TE detection:
     - `--outFilterMultimapNmax 100`
     - `--winAnchorMultimapNmax 100`
     - `--outMultimapperOrder Random`
     - `--runRNGseed 777`
     - `--outSAMmultNmax 1`
     - `--outSAMtype BAM Unsorted`
     - `--outFilterType BySJout`
     - `--alignSJoverhangMin 8`
     - `--alignSJDBoverhangMin 1`
     - `--outFilterMismatchNmax 999`
     - `--alignIntronMin 20`
     - `--alignIntronMax 1000000`
     - `--alignMatesGapMax 1000000`

4. **BAM Processing**
   - Sort and index BAM files using samtools

5. **TEtranscripts Differential Expression Analysis**
   - Compare treatment vs control samples directly using the DESeq method
   - Uses the following inputs:
     - Sorted BAM files (treatment and control)
     - Gene annotation GTF 
     - TE annotation GTF (pre-existing file specified in config.yaml)


## Output Files

The pipeline generates the following outputs:

- **STAR Alignment**:
  - `results/star_index/`: STAR index files
  - `results/star/{sample}/Aligned.out.bam`: Unsorted BAM files
  - `results/star/{sample}/Aligned.sortedByCoord.out.bam`: Sorted BAM files
  - `results/star/{sample}/Aligned.sortedByCoord.out.bam.bai`: BAM index files
  - `results/star/{sample}/Log.final.out`: Alignment statistics

- **TEtranscripts**:
  - `results/tetranscripts/differential_expression.tsv`: Differential expression results
  - `results/tetranscripts/`: Various additional output files from TEtranscripts

## Troubleshooting

- **STAR indexing fails**: Check memory requirements; STAR indexing can require substantial memory.
- **Reference file errors**: Verify that your genome_fasta, gtf_file, and te_gtf_file paths in config.yaml are correct and the files exist.
- **TEtranscripts fails**: Check that your FASTQ files are properly formatted and named.
- **Installation issues**: If you encounter "disk quota exceeded" or memory errors:
  - Use the minimal installation option in the Installation section
  - Avoid R/Bioconductor packages which can be large
  - Install packages one by one to identify problematic dependencies
  - Use mamba instead of conda for more efficient package resolution: `mamba create -n te_minimal ...`

For more detailed troubleshooting, check the log files in the relevant output directories.

## License

This pipeline is available under the MIT License.