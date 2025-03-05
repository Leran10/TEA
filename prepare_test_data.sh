#!/bin/bash
# Script to generate minimal test data for pipeline verification
# This creates small FASTQ files that can be used to test the pipeline workflow

set -e  # Exit on any error

# Create directories
mkdir -p test_data/fastq
mkdir -p test_data/reference

echo "Preparing minimal test files..."

# Create a mini reference genome (just for testing)
cat > test_data/reference/mini_genome.fa << EOF
>chr1
AGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGTACGTACGTACGTACGTAGCTGATCGATCGATGTACGATCG
ATCGATCGTACGTACGTACGTACGTAGCTGATCGATCGATGTACGATCG
ATCGATCGTACGTACGTACGTACGTAGCTGATCGATCGATGTACGATCG
>chr2
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ATCGATCGATCGTACGTACGTACGTAGCTGATCGATCGAGTACGTACGT
ATCGATCGATCGTACGTACGTACGTAGCTGATCGATCGAGTACGTACGT
ATCGATCGATCGTACGTACGTACGTAGCTGATCGATCGAGTACGTACGT
EOF

# Create a mini GTF file
cat > test_data/reference/mini_annotation.gtf << EOF
chr1	ENSEMBL	gene	1	50	.	+	.	gene_id "ENSG00000001"; gene_name "GENE1";
chr1	ENSEMBL	exon	1	50	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST00000001"; gene_name "GENE1";
chr2	ENSEMBL	gene	1	50	.	-	.	gene_id "ENSG00000002"; gene_name "GENE2";
chr2	ENSEMBL	exon	1	50	.	-	.	gene_id "ENSG00000002"; transcript_id "ENST00000002"; gene_name "GENE2";
EOF

# Create small FASTQ files for each sample
for sample in sample1 sample2 sample3 sample4 sample5 sample6; do
  echo "Creating test FASTQ for ${sample}..."
  # R1 file
  cat > test_data/fastq/${sample}_R1.fastq << EOF
@${sample}_read1
AGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@${sample}_read2
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

  # R2 file
  cat > test_data/fastq/${sample}_R2.fastq << EOF
@${sample}_read1
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGACT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@${sample}_read2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

  # Compress files
  gzip test_data/fastq/${sample}_R1.fastq
  gzip test_data/fastq/${sample}_R2.fastq
done

echo "Test data creation complete!"
echo "To use the test data, update config.yaml with:"
echo "reference_genome: \"$(pwd)/test_data/reference/mini_genome.fa\""
echo "gene_annotation: \"$(pwd)/test_data/reference/mini_annotation.gtf\""
echo ""
echo "NOTE: This is minimal test data for workflow validation only."
echo "These files are too small for meaningful biological analysis."