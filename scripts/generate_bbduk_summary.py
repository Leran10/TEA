#!/usr/bin/env python3
"""
Generate summary report for BBDuk processing steps.
This script parses BBDuk log files and generates HTML and TSV summaries.
"""

import os
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def count_reads_in_fastq(fastq_file):
    """Count the number of reads in a FASTQ file."""
    if not os.path.exists(fastq_file):
        return 0
        
    # Use wc -l to count lines and divide by 4 (each read has 4 lines in FASTQ)
    result = os.popen(f"wc -l {fastq_file}").read().strip().split()[0]
    try:
        return int(result) // 4
    except:
        return 0

def parse_bbduk_log(log_file):
    """Parse BBDuk log file to extract stats."""
    if not os.path.exists(log_file):
        return {}
        
    stats = {}
    with open(log_file, 'r') as f:
        for line in f:
            # Look for pattern like "Input:          1000 reads"
            input_match = re.search(r'Input:\s+(\d+)\s+reads', line)
            if input_match:
                stats['input_reads'] = int(input_match.group(1))
                
            # Look for pattern like "Matched:        800 reads (80.00%)"
            matched_match = re.search(r'Matched:\s+(\d+)\s+reads\s+\(([0-9.]+)%\)', line)
            if matched_match:
                stats['matched_reads'] = int(matched_match.group(1))
                stats['matched_percent'] = float(matched_match.group(2))
                
            # Look for pattern like "Result:         800 reads (80.00%)"
            result_match = re.search(r'Result:\s+(\d+)\s+reads\s+\(([0-9.]+)%\)', line)
            if result_match:
                stats['result_reads'] = int(result_match.group(1))
                stats['result_percent'] = float(result_match.group(2))
    
    return stats

def main():
    # Get input parameters from Snakemake
    vector_seq = snakemake.params.vector
    backbone_seq = snakemake.params.backbone
    samples = snakemake.params.samples
    
    bbduk_dir = os.path.dirname(snakemake.output.html)
    logs_dir = os.path.join(bbduk_dir, "logs")
    
    # Initialize summary data
    summary_data = []
    
    # Process each sample
    for sample in samples:
        sample_data = {
            'Sample': sample,
            'Initial Reads': count_reads_in_fastq(f"{sample}.fastq")
        }
        
        # Count reads at each stage
        sample_data['Left Filter Reads'] = count_reads_in_fastq(os.path.join(bbduk_dir, f"{sample}_lfilter.fastq"))
        sample_data['Left Trim Reads'] = count_reads_in_fastq(os.path.join(bbduk_dir, f"{sample}_ltrim.fastq"))
        sample_data['Right Filter Reads'] = count_reads_in_fastq(os.path.join(bbduk_dir, f"{sample}_rfilter.fastq"))
        sample_data['Right Trim Reads'] = count_reads_in_fastq(os.path.join(bbduk_dir, f"{sample}_rtrim.fastq"))
        sample_data['Discarded Reads'] = count_reads_in_fastq(os.path.join(bbduk_dir, f"{sample}_rtrim_miss.fastq"))
        
        # Calculate percentages
        if sample_data['Initial Reads'] > 0:
            sample_data['Left Filter %'] = round(sample_data['Left Filter Reads'] / sample_data['Initial Reads'] * 100, 2)
            sample_data['Left Trim %'] = round(sample_data['Left Trim Reads'] / sample_data['Initial Reads'] * 100, 2)
            sample_data['Right Filter %'] = round(sample_data['Right Filter Reads'] / sample_data['Initial Reads'] * 100, 2)
            sample_data['Right Trim %'] = round(sample_data['Right Trim Reads'] / sample_data['Initial Reads'] * 100, 2)
            sample_data['Discarded %'] = round(sample_data['Discarded Reads'] / sample_data['Initial Reads'] * 100, 2)
            sample_data['Final %'] = round(sample_data['Right Trim Reads'] / sample_data['Initial Reads'] * 100, 2)
        
        # Parse log files for additional stats
        for step in ['lfilter', 'ltrim', 'rfilter', 'rtrim']:
            log_file = os.path.join(logs_dir, f"{sample}_{step}_log.txt")
            stats = parse_bbduk_log(log_file)
            for key, value in stats.items():
                sample_data[f"{step}_{key}"] = value
                
        summary_data.append(sample_data)
    
    # Create a DataFrame
    df = pd.DataFrame(summary_data)
    
    # Save TSV file
    df.to_csv(snakemake.output.stats, sep='\t', index=False)
    
    # Generate plots for HTML report
    fig_dir = os.path.join(os.path.dirname(snakemake.output.html), "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    # Plot read counts at each stage
    plt.figure(figsize=(12, 6))
    plot_df = df.melt(id_vars=['Sample'], 
                     value_vars=['Initial Reads', 'Left Filter Reads', 'Left Trim Reads', 
                                 'Right Filter Reads', 'Right Trim Reads'],
                     var_name='Stage', value_name='Read Count')
    sns.barplot(x='Sample', y='Read Count', hue='Stage', data=plot_df)
    plt.xticks(rotation=45)
    plt.title('Read Counts at Each Processing Stage')
    plt.tight_layout()
    read_count_plot = os.path.join(fig_dir, "read_counts.png")
    plt.savefig(read_count_plot)
    
    # Plot percentages
    plt.figure(figsize=(12, 6))
    if 'Final %' in df.columns:
        sns.barplot(x='Sample', y='Final %', data=df)
        plt.axhline(y=50, color='r', linestyle='--')
        plt.ylim(0, 100)
        plt.title('Percentage of Reads Retained After Processing')
        plt.tight_layout()
        pct_plot = os.path.join(fig_dir, "retention_percentage.png")
        plt.savefig(pct_plot)
    
    # Generate HTML report
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>BBDuk Processing Summary</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1, h2, h3 {{ color: #2c3e50; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            tr:nth-child(even) {{ background-color: #f9f9f9; }}
            .overview {{ background-color: #eef7fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
            .plot {{ text-align: center; margin: 20px 0; }}
            .plot img {{ max-width: 100%; height: auto; }}
        </style>
    </head>
    <body>
        <h1>BBDuk Processing Summary</h1>
        
        <div class="overview">
            <h2>Processing Overview</h2>
            <p><strong>Vector Sequence:</strong> {vector_seq}</p>
            <p><strong>Backbone Sequence:</strong> {backbone_seq}</p>
            <p><strong>Total Samples:</strong> {len(samples)}</p>
            <p><strong>Processing Steps:</strong></p>
            <ol>
                <li><strong>Left Filter:</strong> Filter reads containing the vector sequence</li>
                <li><strong>Left Trim:</strong> Trim off the vector sequence</li>
                <li><strong>Right Filter:</strong> Filter reads containing the backbone sequence</li>
                <li><strong>Right Trim:</strong> Trim off the backbone sequence and keep only reads of the expected length</li>
            </ol>
        </div>
        
        <h2>Processing Results</h2>
        
        <table>
            <tr>
                <th>Sample</th>
                <th>Initial Reads</th>
                <th>Final Reads</th>
                <th>Retention Rate</th>
            </tr>
    """
    
    # Add rows for each sample
    for _, row in df.iterrows():
        html_content += f"""
            <tr>
                <td>{row['Sample']}</td>
                <td>{row['Initial Reads']:,}</td>
                <td>{row['Right Trim Reads']:,}</td>
                <td>{row['Final %']}%</td>
            </tr>
        """
    
    html_content += """
        </table>
        
        <h2>Detailed Statistics</h2>
        
        <div class="plot">
            <h3>Read Counts at Each Processing Stage</h3>
            <img src="figures/read_counts.png" alt="Read Counts Plot">
        </div>
        
        <div class="plot">
            <h3>Percentage of Reads Retained</h3>
            <img src="figures/retention_percentage.png" alt="Retention Percentage Plot">
        </div>
        
        <h2>Complete Data Table</h2>
    """
    
    # Add full table
    html_content += "<table><tr>"
    for col in df.columns:
        html_content += f"<th>{col}</th>"
    html_content += "</tr>"
    
    for _, row in df.iterrows():
        html_content += "<tr>"
        for col in df.columns:
            value = row[col]
            if isinstance(value, (int, float)) and not pd.isna(value):
                if isinstance(value, int):
                    html_content += f"<td>{value:,}</td>"
                else:
                    html_content += f"<td>{value:.2f}</td>"
            else:
                html_content += f"<td>{value}</td>"
        html_content += "</tr>"
    
    html_content += """
        </table>
        
        <h2>Summary</h2>
        <p>This report summarizes the BBDuk processing steps for CRISPR sgRNA read trimming and filtering.</p>
        <p>The final reads represent sgRNAs that have been properly trimmed and are ready for MAGeCK analysis.</p>
        <p>If retention rates are low (&lt; 50%), consider checking:</p>
        <ul>
            <li>Vector and backbone sequences are correct for the library design</li>
            <li>Input FASTQ files are from the expected CRISPR library</li>
            <li>Trimming parameters match the expected sgRNA length</li>
        </ul>
        
    </body>
    </html>
    """
    
    with open(snakemake.output.html, 'w') as f:
        f.write(html_content)
    
    print(f"BBDuk summary report generated: {snakemake.output.html}")
    print(f"BBDuk statistics saved: {snakemake.output.stats}")

if __name__ == "__main__":
    main()