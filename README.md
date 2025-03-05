# TEA (Transposable Element Analysis)

A Snakemake-based pipeline for transposable element family differential abundance analysis.

## Overview

TEA is a comprehensive and user-friendly pipeline for analyzing transposable element (TE) expression in RNA-seq data. The pipeline integrates several powerful tools including STAR for alignment, RepeatMasker for TE annotation, and TEtranscripts for differential abundance analysis.

## Features

- **RNA-seq alignment**: Optimized STAR alignment for transposable element analysis
- **Transposable element annotation**: Identification and classification of TEs using RepeatMasker
- **Differential abundance analysis**: Statistical analysis of TE family expression differences between conditions
- **Quality control**: Integrated MultiQC reporting
- **Reproducibility**: Containerized workflow with Conda environment

## Requirements

- Conda/Mamba
- Snakemake
- 16+ GB RAM (recommended for RepeatMasker and STAR alignment)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/TE_family_analysis.git
cd TE_family_analysis
```

2. Create the Conda environment:
```bash
conda env create -f environment.yaml
```

3. Activate the environment:
```bash
conda activate te_analysis
```

## Configuration

1. Modify the `config.yaml` file with your specific parameters:
   - Sample names
   - Treatment and control groups
   - Reference genome and annotation paths
   - Computational resources allocation

2. Example configuration:
```yaml
# Sample information
samples:
  - sample1
  - sample2
  - sample3

# Sample groups for differential expression analysis
treatment_samples:
  - sample2
  - sample3

control_samples:
  - sample1

# Reference and annotation files
reference_genome: "/path/to/genome.fa"
gene_annotation: "/path/to/genes.gtf"
star_index: "/path/to/star_index" 

# RepeatMasker settings
repeatmasker_species: "human"  # CRITICAL: Must match your organism (human, mouse, rat, etc.)

# Computational resources
threads:
  star: 8
  repeatmasker: 8  
  tetranscripts: 8

memory:
  star: 16000
  repeatmasker: 8000
  tetranscripts: 8000
```

## Input Data Preparation

### Required Reference Files

Before running the pipeline, you need to prepare the following reference files:

#### Recommended Directory Structure

We recommend organizing your reference files as follows:
```
TEA/
├── resources/
│   ├── genome/               # Place reference genomes here
│   │   └── GRCh38.fa
│   ├── annotations/          # Place gene annotations here
│   │   └── gencode.v38.gtf
│   └── star_index/           # Place STAR indices here
│       └── GRCh38/
└── data/
    └── fastq/               # Place input FASTQ files here
```

1. **Reference Genome**: 
   - Download species-specific genome (e.g., human GRCh38, mouse GRCm39)
   - Example source: [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) or [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html)
   - Recommended file format: FASTA (.fa)
   - Recommended location: `resources/genome/`
   - Create directory: `mkdir -p resources/genome`
   - Update `config.yaml` with the absolute path: `reference_genome: "/absolute/path/to/TEA/resources/genome/GRCh38.fa"`

2. **Gene Annotation**:
   - Download GTF file for your reference genome
   - Example source: [GENCODE](https://www.gencodegenes.org/human/) or [Ensembl](https://www.ensembl.org/info/data/ftp/index.html)
   - Recommended file format: GTF (.gtf)
   - Recommended location: `resources/annotations/`
   - Create directory: `mkdir -p resources/annotations`
   - Update `config.yaml` with the absolute path: `gene_annotation: "/absolute/path/to/TEA/resources/annotations/gencode.v38.gtf"`

3. **STAR Index**:
   - If you don't have a pre-built STAR index, you'll need to create one
   - Recommended location: `resources/star_index/`
   - Create directory: `mkdir -p resources/star_index/GRCh38`
   - Update `config.yaml` with the absolute path: `star_index: "/absolute/path/to/TEA/resources/star_index/GRCh38"`
   - To create a STAR index:
     ```bash
     # Create STAR index directory
     mkdir -p resources/star_index/GRCh38
     
     # Generate STAR index
     STAR --runMode genomeGenerate \
          --genomeDir resources/star_index/GRCh38 \
          --genomeFastaFiles resources/genome/GRCh38.fa \
          --sjdbGTFfile resources/annotations/gencode.v38.gtf \
          --sjdbOverhang 100 \
          --runThreadN 8
     ```
   
> **Important**: Always use absolute paths in config.yaml, not relative paths

4. **RepeatMasker Database**:
   - The pipeline uses RepeatMasker to identify transposable elements
   - RepeatMasker requires specific species libraries that should be installed with the tool
   - Verify your RepeatMasker installation includes the species library you need:
     ```bash
     # List available species libraries
     RepeatMasker -species help | grep -i "your_species"
     ```
   - If your species is not available, you may need to:
     1. Install Dfam/RepBase libraries (requires license)
     2. Or use a closely related species instead
     3. Or build custom libraries (advanced users)

### RNA-seq Input Files

Place your paired-end RNA-seq FASTQ files in the `data/fastq/` directory with the naming convention:
- `{sample}_R1.fastq.gz` - Forward reads
- `{sample}_R2.fastq.gz` - Reverse reads

Where `{sample}` corresponds exactly to the sample names defined in your `config.yaml` file.

## Complete Workflow

Here's the step-by-step process to run the entire pipeline from start to finish:

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/TEA.git
   cd TEA
   ```

2. **Create and activate the conda environment**:
   ```bash
   conda env create -f environment.yaml
   conda activate te_analysis
   ```

3. **Prepare reference files**:
   - Download reference genome FASTA file
   - Download gene annotation GTF file
   - Create STAR index if needed (see instructions in Input Data Preparation section)

4. **Edit the configuration file**:
   ```bash
   cp example_config.yaml config.yaml
   # Edit config.yaml with your text editor to add:
   # - Sample names
   # - Treatment/control assignments
   # - Paths to reference files
   # - RepeatMasker species parameter
   ```
   
   The RepeatMasker `species` parameter is critical and must match your reference genome:
   ```yaml
   # Set the correct species for RepeatMasker
   repeatmasker_species: "human"  # Options include: human, mouse, rat, etc.
   ```
   
   Common species values:
   - For human genomes: `human`
   - For mouse genomes: `mouse` or `murinae`
   - For rat genomes: `rat`
   - For other species, check available libraries: `RepeatMasker -species help | less`

5. **Prepare input data**:
   ```bash
   # Create directories
   mkdir -p data/fastq
   
   # Copy or link your FASTQ files to the data directory
   # Make sure they follow the naming convention: {sample}_R1.fastq.gz and {sample}_R2.fastq.gz
   # Where {sample} matches the sample names in your config.yaml
   ```

6. **Run the complete pipeline**:
   ```bash
   # Run everything with a single command
   snakemake --cores <number_of_cores> --use-conda
   ```

7. **Check results**:
   ```bash
   # TE differential analysis results
   less results/tetranscripts/differential/TEcount.DE_results
   
   # QC report
   firefox results/multiqc/multiqc_report.html
   ```

## Usage Options

1. **Perform a dry run** to check the workflow without executing:
   ```bash
   snakemake -n
   ```

2. **Run the pipeline** with specific resource allocation:
   ```bash
   snakemake --cores <number_of_cores> --use-conda
   ```

3. **Run specific parts** of the pipeline:
   ```bash
   # Run only the alignment step
   snakemake --cores 8 results/star/sample1/Aligned.sortedByCoord.out.bam

   # Run only RepeatMasker annotation
   snakemake --cores 8 resources/repeatmasker/TE.gtf

   # Run only differential expression analysis
   snakemake --cores 8 results/tetranscripts/differential
   ```

4. **Resume a failed run**:
   ```bash
   snakemake --cores <number_of_cores> --rerun-incomplete
   ```

5. **Generate a workflow graph**:
   ```bash
   snakemake --dag | dot -Tpdf > workflow.pdf
   ```

6. **Run on a computing cluster**:
   ```bash
   # PBS/Torque
   snakemake --cluster "qsub -l nodes=1:ppn={threads}" --jobs 100

   # SLURM
   snakemake --cluster "sbatch --ntasks={threads} --mem={resources.mem_mb}M" --jobs 100

   # Using profiles (recommended)
   snakemake --profile <your_cluster_profile>
   ```

## Output

The pipeline generates the following outputs:

- `results/star/`: STAR alignment files (BAM and indices)
- `results/repeatmasker/`: RepeatMasker TE annotations
- `results/tetranscripts/differential/`: Differential abundance analysis results
- `results/multiqc/multiqc_report.html`: Quality control report

### Analyzing Results

#### Examining alignment statistics
```bash
# View alignment summary for a sample
cat results/star/sample1/Log.final.out

# Compare mapping rates across samples
multiqc results/star/ -o multiqc_star_only
```

#### Working with TEtranscripts output
```bash
# View differentially expressed TEs with significant p-values
cat results/tetranscripts/differential/TEcount.DE_results | awk '$8 < 0.05'

# Extract top 20 upregulated TE families (sorted by fold change)
cat results/tetranscripts/differential/TEcount.DE_results | sort -k4,4nr | head -n 20 > top_upregulated_TEs.txt

# Plot results using R (example command)
Rscript -e "library(ggplot2); data <- read.table('results/tetranscripts/differential/TEcount.DE_results', header=TRUE); ggplot(data[data\$padj < 0.05,], aes(x=log2FoldChange, y=-log10(padj))) + geom_point() + theme_minimal(); ggsave('volcano_plot.pdf')"
```

#### Visualizing alignments
```bash
# View alignments in a specific genomic region using IGV
igv.sh -g /path/to/reference_genome.fa results/star/sample1/Aligned.sortedByCoord.out.bam results/star/sample2/Aligned.sortedByCoord.out.bam

# Extract reads mapping to a specific TE family
samtools view results/star/sample1/Aligned.sortedByCoord.out.bam | grep "L1HS" > L1HS_reads.txt
```

## Citation

If you use this pipeline, please cite:

- STAR: Dobin et al. (2013)
- RepeatMasker: Smit et al.
- TEtranscripts: Jin et al. (2015)

## Troubleshooting

### Common Issues

#### Errors with RepeatMasker
```
Error: RepeatMasker: command not found
```
- Make sure you've activated the conda environment: `conda activate te_analysis`
- Ensure RepeatMasker is properly installed and configured with repeat libraries
- Check if the RepeatMasker Species is valid with: `RepeatMasker -species help`

#### STAR alignment issues
```
Error: STAR: command not found
```
- Make sure you've activated the conda environment: `conda activate te_analysis`
- Check that the STAR index path is correct in your config file

#### TEtranscripts failures
```
Error: GTF file parsing error
```
- Ensure the GTF files (gene and TE) have the correct format
- Verify that the RepeatMasker output was properly converted to GTF format
- Make sure the script `prepare_te_gtf.py` ran successfully

#### Snakemake workflow issues
```
Error: Missing input files for rule
```
- Run with `--debug` flag to get more detailed error information
- Check that all input files exist in the expected locations
- Verify that your config.yaml has all required fields

### Performance Tips

1. Increase threads for compute-intensive steps in config.yaml
2. For large genomes, consider using a computer with 32+ GB RAM
3. For cluster environments, adjust the resources section in the Snakefile
4. Use the `--resources` flag to limit resources: `snakemake --resources mem_mb=32000`

## License

This project is licensed under the terms of the MIT license.

## Contact

For questions or issues, please open an issue on GitHub or contact the maintainers directly.