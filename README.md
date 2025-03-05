# Transposable Element (TE) Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This Snakemake pipeline performs comprehensive analysis of transposable elements in RNA-seq data using STAR, RepeatMasker, and TEtranscripts. It automates the full workflow from RepeatMasker annotation through differential expression analysis of TEs.

## Pipeline Steps

1. **RepeatMasker TE Annotation**: Runs RepeatMasker on the reference genome to identify transposable elements and create a TE.gtf file:
   ```
   RepeatMasker --species murinae -pa {threads} -gff -dir ./repeatmasker/ genome.fa
   ```

2. **STAR Alignment**: Aligns RNA-seq reads to the reference genome with optimized parameters for TE analysis, including settings for multi-mapped reads

3. **BAM Sorting and Indexing**: Prepares BAM files for downstream analysis

4. **TEtranscripts Differential Expression**: Performs differential expression analysis of transposable elements between treatment and control groups using DESeq:
   ```
   TEtranscripts -t treatment1.bam treatment2.bam treatment3.bam \
                 -c control1.bam control2.bam control3.bam \
                 --TE resources/repeatmasker/TE.gtf \
                 --GTF /home/reference/gene_annotation.gtf \
                 --sortByPos \
                 --outdir ./results \
                 --mode multi \
                 --DESeq
   ```

5. **MultiQC**: Generates quality control reports for all samples

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [STAR](https://github.com/alexdobin/STAR)
- [RepeatMasker](https://www.repeatmasker.org/)
- [TEtranscripts](https://hammelllab.labsites.cshl.edu/software/#TEtranscripts)
- [Samtools](http://www.htslib.org/)
- [MultiQC](https://multiqc.info/)

## Installation

```bash
# Clone this repository
git clone https://github.com/yourusername/TE_family_analysis.git
cd TE_family_analysis

# Install dependencies (if using conda)
conda env create -f environment.yaml
conda activate te_analysis
```

## Usage

1. Edit the `config.yaml` file to specify your samples and paths to reference files
2. Place your FASTQ files in the `data/fastq/` directory with naming pattern `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
3. Run the pipeline:

```bash
# Dry run to check workflow
snakemake -n

# Run the pipeline with specified cores
snakemake --cores 16

# Run on a cluster with SLURM
snakemake --cluster "sbatch -p {partition} -c {threads}" --jobs 10
```

## Configuration

Edit the `config.yaml` file to customize:

- Sample names and grouping (treatment vs. control samples)
- Reference genome path for RepeatMasker and STAR alignment
- Gene annotation GTF path (e.g., Gencode)
- RepeatMasker species selection (murinae)
- TEtranscripts parameters (multi/unique mapping mode, DESeq analysis)
- Computational resources (threads and memory)

## Directory Structure

```
TE_family_analysis/
├── Snakefile            # Main workflow file
├── config.yaml          # Configuration parameters
├── environment.yaml     # Conda environment file
├── scripts/
│   └── prepare_te_gtf.py # Utility script for TE annotation conversion
├── data/
│   └── fastq/           # Input FASTQ files (sample1_R1.fastq.gz, etc.)
├── resources/
│   └── repeatmasker/    # RepeatMasker outputs
│       ├── GRCm39.genome.fa.out  # RepeatMasker raw output
│       └── TE.gtf       # Converted GTF for TEtranscripts
└── results/
    ├── star/            # STAR alignment results
    │   └── sample1/     # Per-sample directories
    │       ├── Aligned.out.bam
    │       └── Aligned.sortedByCoord.out.bam
    ├── tetranscripts/   # TEtranscripts differential expression results
    │   └── differential/
    │       ├── DESeq.gene.txt
    │       └── DESeq.TE.txt
    └── multiqc/         # MultiQC reports
        └── multiqc_report.html
```

## Customization

- Modify the Snakefile to add or remove analysis steps
- Adjust parameters for each tool in the respective rule

## Quick Start

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/te-analysis-pipeline.git
   cd te-analysis-pipeline
   ```

2. Set up your environment with the required dependencies:
   ```bash
   conda env create -f environment.yaml
   conda activate te_analysis
   ```

3. Copy and modify the example config:
   ```bash
   cp example_config.yaml config.yaml
   # Edit config.yaml with your settings
   ```

4. Run the pipeline:
   ```bash
   snakemake --cores 8
   ```

## Test Data

For testing purposes, you can generate minimal test data:
```bash
./prepare_test_data.sh
```

This creates small FASTQ files and reference files that can be used to test the pipeline's workflow.

## Citation

If you use this pipeline in your research, please cite:

```
Wang, L. (2025). TE-Analysis-Pipeline: A Snakemake workflow for transposable element analysis. 
https://github.com/yourusername/te-analysis-pipeline
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.