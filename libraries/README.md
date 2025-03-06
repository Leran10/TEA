# CRISPR Library Files

This directory contains sgRNA library files for use with the CRISPR analysis pipeline.

## Library Format

Each library file should be in MAGeCK-compatible tab-delimited format with the following columns:
- sgRNA ID
- sgRNA sequence
- Gene ID

Example:
```
sgRNA_1	ACGTACGTACGTACGTACGT	GENE1
sgRNA_2	GCTAGCTAGCTAGCTAGCTA	GENE1
sgRNA_3	ATCGATCGATCGATCGATCG	GENE2
```

## Available Libraries

The following libraries are included:

1. **GeCKOv2_Human.txt**
   - Human genome-wide CRISPR knockout library
   - 6 sgRNAs per gene
   - Source: [Addgene #1000000048](https://www.addgene.org/pooled-library/zhang-human-gecko-v2/)

2. **GeCKOv2_Mouse.txt**
   - Mouse genome-wide CRISPR knockout library
   - 6 sgRNAs per gene
   - Source: [Addgene #1000000052](https://www.addgene.org/pooled-library/zhang-mouse-gecko-v2/)

3. **Brunello_Human.txt**
   - Human genome-wide CRISPR knockout library
   - 4 sgRNAs per gene
   - Source: [Addgene #73178](https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/)

## Using Libraries

To use a library in your analysis, specify its path in the config.yaml file:

```yaml
# Single library
library_file: "libraries/GeCKOv2_Human.txt"

# OR use the library_name parameter to point to a predefined library
library_name: "GeCKOv2_Human"
```

## Adding New Libraries

To add a new library:

1. Place the library file in this directory
2. Update this README.md to document the new library
3. Add the library details to the `libraries.yaml` file