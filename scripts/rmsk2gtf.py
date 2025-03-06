#!/usr/bin/env python3
"""
Convert RepeatMasker output to GTF format for use with TEtranscripts.
Usage: python rmsk2gtf.py -i <repeatmasker.out> -o <output.gtf>
"""

import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Convert RepeatMasker out to GTF format")
    parser.add_argument('-i', '--input', required=True, help="RepeatMasker output file")
    parser.add_argument('-o', '--output', required=True, help="Output GTF file")
    return parser.parse_args()

def convert_rmsk_to_gtf(input_file, output_file):
    """Convert RepeatMasker .out file to GTF format for TEtranscripts"""
    
    print(f"Converting {input_file} to GTF format...")
    
    # Track unique TE identifiers
    seen_ids = {}
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write('# GTF file converted from RepeatMasker output\n')
        
        # Skip the header lines in the RepeatMasker output
        for i in range(3):
            next(infile)
            
        # Process each line
        for line in infile:
            fields = line.strip().split()
            
            if len(fields) < 15:
                continue
                
            # Extract relevant information
            score = fields[0]
            chromosome = fields[4]
            start = int(fields[5])
            end = int(fields[6])
            strand = fields[8]
            if strand == '+':
                strand = '+'
            elif strand == 'C':
                strand = '-'
            else:
                strand = '.'
                
            te_class = fields[10]
            te_family = fields[11]
            te_name = fields[10] + ":" + fields[11]
            
            # Create unique ID for this element
            base_id = f"{te_class}:{te_family}:{chromosome}:{start}:{end}"
            if base_id in seen_ids:
                seen_ids[base_id] += 1
                te_id = f"{base_id}.{seen_ids[base_id]}"
            else:
                seen_ids[base_id] = 1
                te_id = f"{base_id}.1"
            
            # Write GTF line
            gtf_line = [
                chromosome,
                "RepeatMasker",
                "exon",
                str(start),
                str(end),
                score,
                strand,
                ".",
                f'gene_id "{te_id}"; transcript_id "{te_id}"; family_id "{te_family}"; class_id "{te_class}";'
            ]
            
            outfile.write('\t'.join(gtf_line) + '\n')
    
    print(f"Conversion complete. Output written to {output_file}")

def main():
    args = parse_args()
    convert_rmsk_to_gtf(args.input, args.output)

if __name__ == "__main__":
    main()