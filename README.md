# TEA - Transposable Element Analysis Pipeline
# TE Family Analysis Pipeline

A Snakemake pipeline for transposable element (TE) differential abundance analysis using STAR, RepeatMasker, and TEtranscripts.

## Overview

This pipeline processes RNA-seq data to identify differentially expressed transposable elements between experimental conditions using a sequential workflow:

1. Automatic download of reference genomes and annotations from Gencode (or use custom references)
2. STAR index creation (if index does not exist)
3. RNA-seq alignment with STAR using specific parameters for TE analysis
4. BAM file sorting and indexing
5. RepeatMasker to identify repetitive elements in the reference genome
6. Convert RepeatMasker output to GTF format for TEtranscripts
7. TEtranscripts analysis comparing treatment vs control groups directly
8. Quality control with FastQC and MultiQC throughout the process

## Prerequisites

- Snakemake
- STAR aligner
- RepeatMasker
- TEtranscripts (including TEcount and TEDiff)
- FastQC
- MultiQC
- Python with relevant packages

## Usage

1. Clone this repository
2. Edit `config.yaml` to specify your samples, reference files, and computational resources
3. Create the necessary directory structure:
   ```
   mkdir -p data/raw_fastq metadata results
   ```
4. Place your paired-end FASTQ files in `data/raw_fastq/` with naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
5. Create a sample metadata file at `metadata/sample_info.txt` following the TEtranscripts format
6. Run the pipeline:
   ```
   snakemake --cores 16
   ```

## Configuration

Edit `config.yaml` to customize the pipeline:

- `samples`: List of all sample IDs
- `treatment_samples`: List of treatment sample IDs
- `control_samples`: List of control sample IDs
- `genome_build`: Choose from predefined Gencode genome builds: "GRCh38", "GRCh37", "GRCm39", "GRCm38", or "custom"
- For custom genomes:
  - `custom_genome`: Path to custom genome FASTA file
  - `custom_gtf`: Path to custom gene annotation GTF file
- `reference_dir`: Directory to store downloaded reference files
- `species`: Species for RepeatMasker (e.g., "human", "mouse")
- `threads`: Computational resources for each rule

## Pipeline Steps and Workflow

The pipeline executes in the following sequence:

1. Download reference genome and GTF annotation if not already available
2. Create STAR index if not already available
3. Align RNA-seq data with STAR (outputs unsorted BAM)
4. Sort and index BAM files
5. Run RepeatMasker on the reference genome
6. Convert RepeatMasker output to GTF format
7. Run TEtranscripts using:
   - Sorted BAM files from treatment and control samples
   - Gene annotation GTF
   - TE annotation GTF derived from RepeatMasker
8. Generate FastQC and MultiQC reports

## Output

Results will be organized in the `results/` directory:
- `results/star_index/`: STAR index files
- `results/star/`: STAR alignment files (both unsorted and sorted BAM)
- `results/repeatmasker/`: RepeatMasker output files
- `results/te_annotation/`: TE annotation in GTF format (converted from RepeatMasker)
- `results/tetranscripts/`: TE differential expression analysis results
- `results/qc/`: Quality control reports (FastQC and MultiQC)

## License

This pipeline is available under the MIT License.