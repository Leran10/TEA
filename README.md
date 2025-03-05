# CRISPR Screen Analysis Pipeline

A Snakemake-based pipeline for CRISPR screen analysis, integrating BBDuk.sh for trimming and MAGeCK for statistical analysis.

## Overview

This pipeline processes CRISPR screening data through the following steps:
1. Read trimming and filtering using BBDuk.sh (from BBTools)
2. sgRNA counting and analysis using MAGeCK

## Requirements

- Conda or Mamba
- Snakemake
- 16+ GB RAM recommended

## Installation

1. Clone this repository:
```bash
git clone https://github.com/username/crispr_analysis.git
cd crispr_analysis
```

2. Create and activate the conda environment:
```bash
conda env create -f environment.yaml
conda activate crispr_analysis
```

## Input Data

### FASTQ Files
Place your FASTQ files in the main directory. The sample names in `config.yaml` should match the prefixes of your FASTQ files.

### CRISPR Library Files
The pipeline provides two ways to specify the sgRNA library:

1. **Use a predefined library by name**:
   - Predefined libraries are stored in the `libraries/` directory
   - Each library includes optimal trimming parameters for that specific library design
   - In `config.yaml`, specify the library by name:
     ```yaml
     # Use a predefined library
     library_name: "GeCKOv2_Human"
     ```

2. **Provide a custom library file**:
   - Create a library file in MAGeCK format (tab-delimited):
     - Column 1: sgRNA ID
     - Column 2: sgRNA sequence
     - Column 3: Gene ID
   - In `config.yaml`, specify the path to your file:
     ```yaml
     # Use a custom library file
     library_file: "path/to/your/library.txt"
     ```

Available predefined libraries include GeCKOv2 (human/mouse) and Brunello. See `libraries/README.md` for a complete list.

## Configuration

Edit the `config.yaml` file to customize the pipeline:

1. **Samples**: List all sample file names (without the `.fastq` extension)
2. **Sample Labels**: Map full sample names to shorter labels for reporting
3. **Comparisons**: Define test vs. control sample comparisons for MAGeCK analysis:
   ```yaml
   comparisons:
     SampleAvsB:                # Name for this comparison (will be used in output filenames)
       test: "SampleB"          # Treatment sample(s) - passed to MAGeCK with -t flag
       control: "SampleA"       # Control sample(s) - passed to MAGeCK with -c flag
   ```
   You can define multiple comparisons, and the pipeline will run separate MAGeCK analyses for each one. For each comparison, specify which sample is the test (experimental) group and which is the control (reference) group.
4. **Library Selection** (choose one option):
   - **Predefined library**:
     ```yaml
     library_name: "GeCKOv2_Human"  # Uses optimal parameters for this library design
     ```
   - **Custom library file**:
     ```yaml
     library_file: "path/to/your/library.txt"
     ```

5. **BBDuk Trimming Parameters** (critical to customize for your CRISPR construct - automatically set if using a predefined library):
   - `vector_sequence`: The 5' vector sequence preceding your sgRNA (for left trimming)
   - `backbone_sequence`: The 3' backbone sequence following your sgRNA (for right trimming)
   - Additional parameters for fine-tuning:
     - `bbduk_kmer_length`: Length of k-mers for matching (default: 8)
     - `bbduk_reverse_complement`: Whether to consider reverse complement (f=false, t=true)
     - `bbduk_mismatch`: Allow mismatches? (f=false, t=true)
     - `bbduk_hamming_distance`: Number of mismatches allowed (0=exact match)
     - `bbduk_min_length`/`bbduk_max_length`: Expected length of sgRNA after trimming

6. **Resources**: Configure computational resources (threads, memory)

## Running the Pipeline

1. **Dry run** to check workflow without execution:
```bash
snakemake -n
```

2. **Run the pipeline**:
```bash
snakemake --cores <number_of_cores>
```

3. **Run with conda environments** (recommended):
```bash
snakemake --cores <number_of_cores> --use-conda
```

4. **For cluster execution** (SLURM example):
```bash
snakemake --cluster "sbatch --ntasks={threads} --mem={resources.mem_mb}M" --jobs 20
```

## Pipeline Steps

### BBDuk Trimming
The pipeline trims reads in 4 steps:
1. **Left filter**: Filter reads containing the vector sequence
2. **Left trim**: Trim off the vector sequence
3. **Right filter**: Filter reads containing the backbone sequence
4. **Right trim**: Trim off the backbone sequence and keep only reads of the expected length

### MAGeCK Analysis
1. **Count**: Count sgRNAs in each sample using the processed FASTQ files
2. **Test**: Perform statistical analysis of sgRNA enrichment/depletion

For the MAGeCK test step, the pipeline determines which samples are test (treatment) and which are control samples from your `config.yaml` file. In the configuration file, you define comparisons like this:

```yaml
comparisons:
  VirusVsControl:       # Name of the comparison
    test: "Virus"       # Treatment sample (will be passed to MAGeCK with -t flag)
    control: "Mock"     # Control sample (will be passed to MAGeCK with -c flag)
  DrugVsVirus:
    test: "Drug"
    control: "Virus"
```

You can also specify multiple samples for each group by using a list:

```yaml
comparisons:
  MultipleReplicates:
    test: ["VirusRep1", "VirusRep2", "VirusRep3"]    # Multiple treatment replicates
    control: ["MockRep1", "MockRep2", "MockRep3"]    # Multiple control replicates
```

The pipeline will automatically run separate MAGeCK test analyses for each comparison you define. This allows you to perform multiple comparisons in a single pipeline run.

Additional MAGeCK parameters can be configured in the config.yaml file:

```yaml
# MAGeCK parameters
mageck_fdr: 0.05                     # False discovery rate threshold
mageck_normalization_method: "median"  # Options: median, total, control
```

## Output

Results are saved in the `results` directory (or as configured in `config.yaml`):

- **BBDuk Results**: `results/bbduk_output/`
  - `*_rtrim.fastq`: Final trimmed FASTQ files
  - `logs/`: BBDuk log files
  - `bbduk_summary_report.html`: Interactive HTML report with statistics and visualizations
  - `bbduk_summary_stats.tsv`: Tab-separated file with detailed statistics
  - `figures/`: Directory containing plots from the summary report

- **MAGeCK Results**: `results/mageck_output/`
  - `mageck_counts.count.txt`: sgRNA count table
  - `*.gene_summary.txt`: Gene-level statistics for each comparison
  - `*.sgrna_summary.txt`: sgRNA-level statistics for each comparison

### BBDuk Summary Report

The pipeline automatically generates a comprehensive summary report after the BBDuk processing steps. This HTML report includes:

1. **Processing Statistics**:
   - Initial read counts for each sample
   - Reads retained at each processing step
   - Final retention rates
   - Detailed statistics from BBDuk logs

2. **Visualizations**:
   - Bar plots showing read counts at each stage
   - Retention percentage for each sample
   - Reference line at 50% to identify potential issues

3. **Troubleshooting Guidance**:
   - Suggestions for improving retention rates
   - Links to diagnostic information

This report helps you diagnose any potential issues with the trimming process and verify that the library-specific trimming parameters are correct for your data.

## Customizing for Different CRISPR Libraries

The most important parameters to customize are the vector and backbone sequences that surround your sgRNA sequence. Different CRISPR libraries use different vector designs. Here are examples for common CRISPR libraries:

### Example: Customizing for GeCKO v2 Library
```yaml
vector_sequence: "CACCG"          # 5' vector sequence
backbone_sequence: "GTTTAAGAGC"   # 3' backbone sequence
bbduk_min_length: 20              # sgRNA length
bbduk_max_length: 20
```

### Example: Customizing for Brunello Library
```yaml
vector_sequence: "AAACACCG"
backbone_sequence: "GTTTAAGAGC"
bbduk_min_length: 20
bbduk_max_length: 20
```

### Example: Allowing 1 Mismatch
If your data has sequencing errors, you might want to allow 1 mismatch:
```yaml
bbduk_mismatch: "t"
bbduk_hamming_distance: 1
```

## Troubleshooting

- **BBDuk errors**: Check logs in `results/bbduk_output/logs/`
- **Missing library file**: Ensure the CRISPR library file exists and is formatted correctly
- **Memory issues**: Increase memory in `config.yaml` for bbduk or mageck as needed
- **Low sgRNA counts**: Verify that the vector and backbone sequences match your library design

## Citation

If you use this pipeline, please cite:
- BBTools/BBDuk: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/
- MAGeCK: Li W, et al. "MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens." Genome Biology 15.12 (2014): 554.
- Snakemake: Köster, J. and Rahmann, S. "Snakemake - A scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.