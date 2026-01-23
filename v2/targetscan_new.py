#!/usr/bin/env python3
"""
Parallel TargetScan miRNA Target Prediction Tool
Optimized for multi-core execution on pre-split UTR datasets
"""

import os
import argparse
import multiprocessing
import subprocess
import json
from typing import List, Dict

# Configuration
TARGETSCAN_DIR = "/opt/TargetScan"
TARGETSCAN_SCRIPT1 = os.path.join(TARGETSCAN_DIR, "TargetScan_70/targetscan_70.pl")
TARGETSCAN_SCRIPT2 = os.path.join(TARGETSCAN_DIR, "TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl")
UTR_PARTS_DIR = os.path.join(TARGETSCAN_DIR, "Datasets/3utr/")
BLN_BINS_DIR = os.path.join(TARGETSCAN_DIR, "Datasets/bln_bins/")
MIR_FAMILY_INFO = os.path.join(TARGETSCAN_DIR, "Datasets/miR_Family_Info.json")
NUM_UTR_PARTS = 64  # Default number of UTR parts in TargetScan dataset

def parse_fasta(fasta_file: str) -> List[Dict[str, str]]:
    """Parse miRNA FASTA file and filter valid sequences (17-30nt)."""
    sequences = []
    current_header = ""
    current_sequence = ""
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and 17 <= len(current_sequence) <= 30:
                    sequences.append({
                        "header": current_header,
                        "sequence": current_sequence
                    })
                current_header = line[1:]  # Remove '>'
                current_sequence = ""
            else:
                current_sequence += line.upper().replace("T", "U")  # Convert to RNA
        
        # Add last sequence
        if current_header and 17 <= len(current_sequence) <= 30:
            sequences.append({
                "header": current_header,
                "sequence": current_sequence
            })
    
    return sequences

def prepare_targetscan_input(sequence: str, header: str, output_dir: str) -> str:
    """
    Create TargetScan input file for a single miRNA.
    Returns path to created input file.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Get seed region (positions 2-8)
    seed = sequence[1:8]
    
    # Load miRNA family info
    with open(MIR_FAMILY_INFO, 'r') as f:
        mir_family_info = json.load(f)
    
    # Determine species ID
    species_id = ';'.join(mir_family_info[seed]) if seed in mir_family_info else '9606'
    
    # Format identifier (TargetScan-specific)
    identifier = header.split(",")[0].replace('hsa-', '')
    
    # Write input file
    input_path = os.path.join(output_dir, f"{header}_targetscan.txt")
    with open(input_path, 'w') as f:
        f.write(f"{identifier}\t{seed}\t{species_id}\n")
    
    return input_path

def process_utr_part(args: tuple):
    """
    Process a single UTR part with TargetScan (worker function for parallel processing).
    Args: (mirna_input_path, part_id, output_dir)
    """
    mirna_input, part_id, output_dir = args
    
    # Define file paths
    utr_file = os.path.join(UTR_PARTS_DIR, f"targetscan_utr_part_{part_id}.txt")
    bln_bins_file = os.path.join(BLN_BINS_DIR, f"targetscan_median_bls_bins_part_{part_id}.txt")
    out1 = os.path.join(output_dir, f"part_{part_id}_out1.txt")
    out2 = os.path.join(output_dir, f"part_{part_id}_out2.txt")
    
    # Run TargetScan Step 1
    subprocess.run([
        "perl", TARGETSCAN_SCRIPT1,
        mirna_input, utr_file, out1
    ], check=True)
    
    # Run TargetScan Step 2
    with open(out2, 'w') as f_out2, open(os.devnull, 'w') as devnull:
        subprocess.run([
            "perl", TARGETSCAN_SCRIPT2,
            mirna_input, out1, bln_bins_file
        ], stdout=f_out2, stderr=devnull, check=True)

def merge_partial_results(output_dir: str, header: str):
    """
    Merge partial results from parallel processing into final output files.
    """
    # File paths for final outputs
    final_out1 = os.path.join(output_dir, f"{header}_Targetscan_output_sort.txt")
    final_out2 = os.path.join(output_dir, f"{header}_Targetscan_output.txt")
    
    # Merge Step 1 results
    with open(final_out1, 'w') as f_out1:
        for part_id in range(NUM_UTR_PARTS):
            part_file = os.path.join(output_dir, f"part_{part_id}_out1.sort.txt")
            if os.path.exists(part_file):
                with open(part_file, 'r') as f_part:
                    if part_id == 0:
                        f_out1.write(f_part.readline())  # Keep header from first file
                    else:
                        f_part.readline()  # Skip headers in other files
                    f_out1.write(f_part.read())
                os.remove(part_file)
            part_file = os.path.join(output_dir, f"part_{part_id}_out1.txt")
            if os.path.exists(part_file):
                os.remove(part_file)
    
    # Merge Step 2 results
    with open(final_out2, 'w') as f_out2:
        for part_id in range(NUM_UTR_PARTS):
            part_file = os.path.join(output_dir, f"part_{part_id}_out2.txt")
            if os.path.exists(part_file):
                with open(part_file, 'r') as f_part:
                    if part_id == 0:
                        f_out2.write(f_part.readline())  # Keep header
                    else:
                        f_part.readline()  # Skip header
                    f_out2.write(f_part.read())
                os.remove(part_file)

def run_targetscan_parallel(sequences: List[Dict[str, str]], num_cores: int, output_dir: str):
    """
    Run TargetScan prediction for all miRNAs using parallel processing.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for seq in sequences:
        print(f"Processing miRNA: {seq['header']}")
        
        # Step 1: Prepare TargetScan input file
        mirna_input = prepare_targetscan_input(
            seq['sequence'],
            seq['header'],
            output_dir
        )
        
        # Step 2: Process all UTR parts in parallel
        with multiprocessing.Pool(num_cores) as pool:
            pool.map(process_utr_part, [
                (mirna_input, part_id, output_dir)
                for part_id in range(NUM_UTR_PARTS)
            ])
        
        # Step 3: Merge partial results
        merge_partial_results(output_dir, seq['header'])
        
        # Cleanup
        os.remove(mirna_input)
        
        print(f"Completed processing: {seq['header']}")

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Parallel TargetScan miRNA Target Prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file containing miRNA sequences")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for results")
    parser.add_argument("-c", "--cores", type=int, default=multiprocessing.cpu_count(),
                        help="Number of CPU cores to use")
    
    args = parser.parse_args()
    
    # Validate paths
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    # Process sequences
    print("Parsing miRNA sequences...")
    sequences = parse_fasta(args.input)
    print(f"Found {len(sequences)} valid miRNA sequences")
    
    # Run TargetScan
    print("Starting TargetScan predictions...")
    run_targetscan_parallel(sequences, args.cores, args.output)
    
    print("\nTargetScan processing completed successfully!")
    print(f"Results saved to: {args.output}")

if __name__ == "__main__":
    main()