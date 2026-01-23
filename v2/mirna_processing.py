#!/usr/bin/env python3
"""
MicroRNA Sequence Manipulation Tool
===================================

A comprehensive tool for retrieving and modifying microRNA sequences with support for:
- Precursor sequence extraction
- Nucleotide modifications
- Sequence shifting
- Combined modifications and shifts
"""

import argparse
import json
import os
import sys
from typing import Tuple, List, Optional, Dict, Union

# Version information
__version__ = "0.0.2"

# Constants
MARUTE_PRE_MIRNA = '/opt/resources/mature_pre_mirna_ext.json'
MIN_LENGTH = 17
MAX_LENGTH = 30
VALID_NUCLEOTIDES = {'A', 'T', 'C', 'G', 'U'}

def load_json(file_path: str) -> Dict:
    """Load and validate JSON data from file."""
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        if not isinstance(data, dict):
            raise ValueError("JSON data should be a dictionary")
        return data
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{file_path}' was not found.")
    except json.JSONDecodeError:
        raise ValueError(f"The file '{file_path}' is not a valid JSON file.")

def find_mirna_sequence(data: Dict, mirna_id: str, pre_id: Optional[str] = None) -> Tuple[Optional[str], Optional[str], Optional[int], Optional[int], List[str]]:
    """Retrieve mature and precursor sequences with validation.
    Returns: (mature_seq, pre_seq, start, end, available_pre_ids)
    """
    if mirna_id not in data:
        return None, None, None, None, []
    
    entries = data[mirna_id]
    available_pre_ids = [entry.get('pre_id', f"precursor_{i+1}") for i, entry in enumerate(entries)]
    
    # If no pre_id specified and multiple options exist, return available pre_ids
    if not pre_id and len(entries) > 1:
        return None, None, None, None, available_pre_ids
    
    # Select the appropriate entry
    if pre_id:
        selected_entry = None
        for entry in entries:
            if entry.get('pre_id') == pre_id:
                selected_entry = entry
                break
        if not selected_entry:
            raise ValueError(f"Precursor ID '{pre_id}' not found for {mirna_id}")
    else:
        selected_entry = entries[0]
    
    # Extract sequence information
    required_fields = ['mature_seq', 'ext_pre_seq', 'ext_mature_loc_start', 'ext_mature_loc_end']
    
    if not all(field in selected_entry for field in required_fields):
        raise ValueError(f"Missing required fields in data for {mirna_id}")
    
    mature_seq = selected_entry['mature_seq'].upper()
    pre_seq = selected_entry['ext_pre_seq'].upper()
    start = selected_entry['ext_mature_loc_start']
    end = selected_entry['ext_mature_loc_end']
    
    # Validate sequences
    for seq, name in [(mature_seq, 'mature'), (pre_seq, 'precursor')]:
        if not all(nuc in VALID_NUCLEOTIDES for nuc in seq):
            raise ValueError(f"Invalid nucleotides found in {name} sequence for {mirna_id}")
    
    return mature_seq, pre_seq, start, end, available_pre_ids

def prompt_pre_id_selection(available_pre_ids: List[str]) -> str:
    """Prompt user to select a precursor ID from available options."""
    print(f"Multiple precursor sequences available for this miRNA:")
    for i, pre_id in enumerate(available_pre_ids, 1):
        print(f"{i}. {pre_id}")
    
    while True:
        try:
            selection = input("Enter the number of your choice: ")
            idx = int(selection) - 1
            if 0 <= idx < len(available_pre_ids):
                return available_pre_ids[idx]
            print(f"Please enter a number between 1 and {len(available_pre_ids)}")
        except ValueError:
            print("Please enter a valid number")

def validate_modification(mod: str, max_pos: int) -> Tuple[int, str, str]:
    """Validate and parse a modification string."""
    try:
        position_str, change = mod.split(':')
        original, new = change.split('|')
        position = int(position_str)
        
        if position < 1 or position > max_pos:
            raise ValueError(f"Position {position} is out of range (1-{max_pos})")
        if original.upper() not in VALID_NUCLEOTIDES or new.upper() not in VALID_NUCLEOTIDES:
            raise ValueError("Nucleotides must be A, T, C, G, or U")
            
        return position - 1, original.upper(), new.upper()  # Convert to 0-based index
    except ValueError as e:
        raise ValueError(f"Invalid modification '{mod}': {str(e)}")

def apply_modifications(sequence: str, modifications: List[str]) -> Tuple[str, List[str]]:
    """Apply multiple nucleotide modifications to a sequence."""
    modified_seq = sequence
    successful_mods = []
    
    for mod in modifications:
        try:
            pos, original, new = validate_modification(mod, len(sequence))
            if modified_seq[pos] != original:
                    raise ValueError(f"Expected '{original}' at position {pos + 1}, found '{modified_seq[pos]}'")
                
            modified_seq = modified_seq[:pos] + new + modified_seq[pos + 1:]
            successful_mods.append(mod)
        except ValueError as e:
            print(f"Warning: {str(e)}. Skipping modification.")
    
    return modified_seq, successful_mods

def apply_shift(pre_seq: str, mature_seq: str, start: int, end: int, shift: str) -> Tuple[str, str]:
    """Apply sequence shift based on precursor coordinates."""
    try:
        left_shift, right_shift = map(int, shift.split('|'))
        new_start = start + left_shift
        new_end = end + right_shift
        
        if new_start < 0 or new_end > len(pre_seq) or new_start >= new_end:
            raise ValueError("Shift would result in invalid sequence coordinates")
            
        shifted_seq = pre_seq[new_start:new_end]
        return shifted_seq, f"{left_shift}|{right_shift}"
    except ValueError as e:
        raise ValueError(f"Invalid shift '{shift}': {str(e)}")

def handle_combined_operation(pre_seq: str, mature_seq: str, start: int, end: int,
                            modifications: List[str], shift: str) -> Tuple[str, str]:
    """Handle combined shift and modification with position adjustment."""
    # First apply shift
    shifted_seq, shift_info = apply_shift(pre_seq, mature_seq, start, end, shift)
    # left_shift = int(shift.split('|')[0])
    # Then apply modifications with adjusted positions
    modified_seq = shifted_seq
    successful_mods = []
    
    for mod in modifications:
        try:
            pos_str, change = mod.split(':')
            original, new = change.split('|')
            original_pos = int(pos_str)
            adjusted_pos = original_pos  - 1 #- left_shift  # Adjust for shift
            
            if adjusted_pos < 0 or adjusted_pos >= len(shifted_seq):
                raise ValueError(f"Adjusted position {adjusted_pos + 1} out of bounds")
            # Handle T/U conversion case
            current_base = modified_seq[adjusted_pos]
            expected_base = original.upper()
            if current_base == "U" and expected_base == "T":
                pass
            elif current_base != expected_base:          
                raise ValueError(f"Expected '{expected_base}' at position {original_pos}, found '{current_base}'")
            # Apply modification    
            modified_seq = modified_seq[:adjusted_pos] + new.upper() + modified_seq[adjusted_pos + 1:]
            successful_mods.append(mod)
        except ValueError as e:
            print(f"Warning: {str(e)}. Skipping modification: {mod}")
    
    if successful_mods:
        return modified_seq, f"{'&'.join(successful_mods)},{shift_info}"
    return shifted_seq, shift_info

def validate_sequence_length(sequence: str, name: str) -> bool:
    """Check if sequence length is within recommended range."""
    length = len(sequence)
    if length < MIN_LENGTH or length > MAX_LENGTH:
        print(f"Warning: {name} sequence length {length} is outside recommended range ({MIN_LENGTH}-{MAX_LENGTH})")
        return False
    return True

def write_fasta(output_path: str, sequences: List[Tuple[str, str]]) -> None:
    """Write sequences to FASTA file with validation."""
    try:
        with open(output_path, 'w') as f:
            for header, seq in sequences:
                if not header.startswith('>'):
                    header = f">{header}"
                f.write(f"{header}\n{seq}\n")
    except IOError as e:
        raise IOError(f"Failed to write output file: {str(e)}")

def main():
    parser = argparse.ArgumentParser(
        description="MicroRNA Sequence Manipulation Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('mirna_id', help="microRNA ID (e.g., hsa-miR-495-3p)")
    parser.add_argument('-o', '--output', default='mirna.fa',
                      help="Output FASTA file path")
    parser.add_argument('-m', '--modification', action='append',
                      help="Modification in format 'position:original|new' (e.g., '5:A|G')")
    parser.add_argument('-s', '--shift',
                      help="Shift in format 'left|right' (e.g., '-4|-6')")
    parser.add_argument('-b', '--both', action='store_true',
                      help="Apply both modifications and shift (requires both -m and -s)")
    parser.add_argument('--pre-id',
                      help="Specify precursor ID when multiple options exist")
    parser.add_argument('--strict', action='store_true',
                      help="Exit if any modification or shift fails")
    
    args = parser.parse_args()

    try:
        # Validate arguments
        if args.both and not (args.modification and args.shift):
            raise ValueError("--both requires both --modification and --shift")
        
        # Load and validate data
        data = load_json(MARUTE_PRE_MIRNA)
        mature_seq, pre_seq, start, end, available_pre_ids = find_mirna_sequence(data, args.mirna_id, args.pre_id)
        
        # Handle case where multiple precursors exist
        if mature_seq is None and available_pre_ids:
            if args.pre_id:
                raise ValueError(f"Specified precursor ID '{args.pre_id}' not found. Available options: {', '.join(available_pre_ids)}")
            
            if not sys.stdin.isatty():
                raise ValueError(f"Multiple precursor sequences available for {args.mirna_id}. Please specify one with --pre-id: {', '.join(available_pre_ids)}")
            
            selected_pre_id = prompt_pre_id_selection(available_pre_ids)
            mature_seq, pre_seq, start, end, _ = find_mirna_sequence(data, args.mirna_id, selected_pre_id)
        
        if not mature_seq:
            raise ValueError(f"No sequence found for {args.mirna_id}")
        
        validate_sequence_length(mature_seq, "Mature")
        
        # Prepare sequences
        sequences = [(f"{args.mirna_id},WT", mature_seq)]
        
        # Handle operations based on arguments
        if args.both:
            mod_shift_seq, info = handle_combined_operation(
                pre_seq, mature_seq, start, end, args.modification, args.shift
            )
            if validate_sequence_length(mod_shift_seq, "Modified-shifted"):
                sequences.append((f"{args.mirna_id},{info},modified_shifted", mod_shift_seq))
        else:
            if args.modification:
                mod_seq, mods = apply_modifications(mature_seq, args.modification)
                if mods and validate_sequence_length(mod_seq, "Modified"):
                    sequences.append((f"{args.mirna_id},{'&'.join(mods)},modified", mod_seq))
            
            if args.shift:
                shifted_seq, shift_info = apply_shift(pre_seq, mature_seq, start, end, args.shift)
                if validate_sequence_length(shifted_seq, "Shifted"):
                    sequences.append((f"{args.mirna_id},{shift_info},shifted", shifted_seq))
        
        # Write output
        write_fasta(args.output, sequences)
        print(f"Successfully wrote {len(sequences)} sequences to {args.output}")
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()