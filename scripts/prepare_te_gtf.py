#!/usr/bin/env python3
"""
Script to convert RepeatMasker output to a GTF file suitable for TEtranscripts.
Usage: python prepare_te_gtf.py <input_rm_out> <output_gtf>
"""

import sys
import os
import pandas as pd

def convert_repeatmasker_to_gtf(input_file, output_file):
    """
    Convert RepeatMasker .out file to GTF format.
    """
    if not os.path.exists(input_file):
        print(f"Input file {input_file} does not exist.", file=sys.stderr)
        sys.exit(1)
    
    # Skip the header lines in RepeatMasker .out file
    with open(input_file, 'r') as f:
        for i in range(3):  # Skip the first 3 lines
            next(f)
        
        # Prepare output GTF file
        with open(output_file, 'w') as out:
            for line in f:
                fields = line.strip().split()
                if len(fields) < 15:
                    continue
                
                # Extract relevant information from RepeatMasker output
                seqname = fields[4]
                start = int(fields[5])
                end = int(fields[6])
                strand = fields[8] if fields[8] in ['+', '-'] else '+'
                repeat_class = fields[10]
                repeat_family = fields[9]
                
                # Write in GTF format
                gtf_line = (
                    f"{seqname}\tRepeatMasker\texon\t{start}\t{end}\t.\t{strand}\t.\t"
                    f"gene_id \"{repeat_family}\"; transcript_id \"{repeat_family}:{seqname}:{start}-{end}\"; "
                    f"family_id \"{repeat_family}\"; class_id \"{repeat_class}\";\n"
                )
                out.write(gtf_line)
    
    print(f"Converted {input_file} to {output_file} in GTF format.")

def create_te_gtf_from_repeatmasker_database(input_file, output_file):
    """
    Alternative method: Create a TE GTF directly from RepeatMasker database.
    """
    # This is a placeholder - in practice you would parse RepeatMasker's library
    # and generate a GTF file with TE annotations.
    pass

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prepare_te_gtf.py <input_rm_out> <output_gtf>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    convert_repeatmasker_to_gtf(input_file, output_file)