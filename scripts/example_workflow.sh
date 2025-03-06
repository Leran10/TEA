#!/bin/bash
# Example workflow script to demonstrate how to run the pipeline

# Create conda environment
echo "Creating conda environment from environment.yaml..."
conda env create -f environment.yaml

# Activate environment
echo "Activating conda environment..."
conda activate te_analysis

# Create directory structure
echo "Creating directory structure..."
mkdir -p data/raw_fastq metadata results

# Download test data (this is just a placeholder - replace with actual data)
echo "In a real workflow, you would place your FASTQ files in data/raw_fastq/"
echo "Example file naming: data/raw_fastq/sample1_R1.fastq.gz, data/raw_fastq/sample1_R2.fastq.gz"

# Edit the config file
echo "Editing config.yaml to match your setup..."
echo "Key parameters to check:"
echo "  1. genome_build: Choose from GRCh38, GRCh37, GRCm39, GRCm38 (or custom)"
echo "  2. treatment_samples and control_samples: Ensure these match your experiment"
echo "  3. species: For RepeatMasker (must match your genome)"

# Create sample metadata file
echo "Creating sample metadata file..."
cat > metadata/sample_info.txt << EOF
SampleID	Condition	Batch
sample1	control	1
sample2	treatment	1
sample3	treatment	1
sample4	control	1
sample5	treatment	1
sample6	control	1
EOF

# Run the pipeline
echo "Running the Snakemake pipeline..."
echo "In a real workflow, you would run: snakemake --cores <number_of_cores>"
echo "For a dry run: snakemake -n"

# When everything is ready, you can run:
# snakemake --cores 16