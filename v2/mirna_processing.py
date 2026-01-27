#!/usr/bin/env python3
from __future__ import print_function
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


# Version information
__version__ = "0.0.2"

# Constants
MARUTE_PRE_MIRNA = '/opt/resources/mature_pre_mirna_ext.json'
MIN_LENGTH = 17
MAX_LENGTH = 30
VALID_NUCLEOTIDES = {'A', 'T', 'C', 'G', 'U'}

def load_json(file_path):
    """Load and validate JSON data from file."""
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        if not isinstance(data, dict):
            raise ValueError("JSON data should be a dictionary")
        return data
    except FileNotFoundError:
        raise FileNotFoundError("The file '{}' was not found.".format(file_path))
    except json.JSONDecodeError:
        raise ValueError("The file '{}' is not a valid JSON file.".format(file_path))

def find_mirna_sequence(data, mirna_id, pre_id=None):
    """Retrieve mature and precursor sequences with validation.
    Returns: (mature_seq, pre_seq, start, end, available_pre_ids)
    """
    if mirna_id not in data:
        return None, None, None, None, []
    
    entries = data[mirna_id]
    available_pre_ids = [entry.get('pre_id', "precursor_{}".format(i + 1)) for i, entry in enumerate(entries)]
    
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
            raise ValueError("Precursor ID '{}' not found for {}".format(pre_id, mirna_id))
    else:
        selected_entry = entries[0]
    
    # Extract sequence information
    required_fields = ['mature_seq', 'ext_pre_seq', 'ext_mature_loc_start', 'ext_mature_loc_end']
    
    if not all(field in selected_entry for field in required_fields):
        raise ValueError("Missing required fields in data for {}".format(mirna_id))
    
    mature_seq = selected_entry['mature_seq'].upper()
    pre_seq = selected_entry['ext_pre_seq'].upper()
    start = selected_entry['ext_mature_loc_start']
    end = selected_entry['ext_mature_loc_end']
    
    # Validate sequences
    for seq, name in [(mature_seq, 'mature'), (pre_seq, 'precursor')]:
        if not all(nuc in VALID_NUCLEOTIDES for nuc in seq):
            raise ValueError("Invalid nucleotides found in {} sequence for {}".format(name, mirna_id))
    
    return mature_seq, pre_seq, start, end, available_pre_ids

def prompt_pre_id_selection(available_pre_ids):
    """Prompt user to select a precursor ID from available options."""
    print("Multiple precursor sequences available for this miRNA:")
    for i, pre_id in enumerate(available_pre_ids, 1):
        print("{}. {}".format(i, pre_id))
    
    while True:
        try:
            selection = input("Enter the number of your choice: ")
            idx = int(selection) - 1
            if 0 <= idx < len(available_pre_ids):
                return available_pre_ids[idx]
            print("Please enter a number between 1 and {}".format(len(available_pre_ids)))
        except ValueError:
            print("Please enter a valid number")

def validate_modification(mod, max_pos):
    """Validate and parse a modification string."""
    try:
        position_str, change = mod.split(':')
        original, new = change.split('|')
        position = int(position_str)
        
        if position < 1 or position > max_pos:
            raise ValueError("Position {} is out of range (1-{})".format(position, max_pos))
        if original.upper() not in VALID_NUCLEOTIDES or new.upper() not in VALID_NUCLEOTIDES:
            raise ValueError("Nucleotides must be A, T, C, G, or U")
            
        return position - 1, original.upper(), new.upper()  # Convert to 0-based index
    except ValueError as e:
        raise ValueError("Invalid modification '{}': {}".format(mod, str(e)))

def apply_modifications(sequence, modifications):
    """Apply multiple nucleotide modifications to a sequence."""
    modified_seq = sequence
    successful_mods = []
    
    for mod in modifications:
        try:
            pos, original, new = validate_modification(mod, len(sequence))
            if modified_seq[pos] != original:
                raise ValueError("Expected '{}' at position {}, found '{}'".format(original, pos + 1, modified_seq[pos]))

            modified_seq = modified_seq[:pos] + new + modified_seq[pos + 1:]
            successful_mods.append(mod)
        except ValueError as e:
            print("Warning: {}. Skipping modification.".format(str(e)))
    
    return modified_seq, successful_mods

def apply_shift(pre_seq, mature_seq, start, end, shift):
    """Apply sequence shift based on precursor coordinates."""
    try:
        left_shift, right_shift = map(int, shift.split('|'))
        new_start = start + left_shift
        new_end = end + right_shift
        
        if new_start < 0 or new_end > len(pre_seq) or new_start >= new_end:
            raise ValueError("Shift would result in invalid sequence coordinates")
            
        shifted_seq = pre_seq[new_start:new_end]
        return shifted_seq, "{}|{}".format(left_shift, right_shift)
    except ValueError as e:
        raise ValueError("Invalid shift '{}': {}".format(shift, str(e)))

def handle_combined_operation(pre_seq, mature_seq, start, end, modifications, shift):
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
                raise ValueError("Adjusted position {} out of bounds".format(adjusted_pos + 1))
            # Handle T/U conversion case
            current_base = modified_seq[adjusted_pos]
            expected_base = original.upper()
            if current_base == "U" and expected_base == "T":
                pass
            elif current_base != expected_base:          
                raise ValueError("Expected '{}' at position {}, found '{}'".format(expected_base, original_pos, current_base))
            # Apply modification    
            modified_seq = modified_seq[:adjusted_pos] + new.upper() + modified_seq[adjusted_pos + 1:]
            successful_mods.append(mod)
        except ValueError as e:
            print("Warning: {}. Skipping modification: {}".format(str(e), mod))
    
    if successful_mods:
        return modified_seq, "{},{}".format('&'.join(successful_mods), shift_info)
    return shifted_seq, shift_info

def validate_sequence_length(sequence, name):
    """Check if sequence length is within recommended range."""
    length = len(sequence)
    if length < MIN_LENGTH or length > MAX_LENGTH:
        print("Warning: {} sequence length {} is outside recommended range ({}-{})".format(name, length, MIN_LENGTH, MAX_LENGTH))
        return False
    return True

def write_fasta(output_path, sequences):
    """Write sequences to FASTA file with validation."""
    try:
        with open(output_path, 'w') as f:
            for header, seq in sequences:
                if not header.startswith('>'):
                    header = ">{}".format(header)
                f.write("{}\n{}\n".format(header, seq))
    except IOError as e:
        raise IOError("Failed to write output file: {}".format(str(e)))

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
                raise ValueError("Specified precursor ID '{}' not found. Available options: {}".format(args.pre_id, ', '.join(available_pre_ids)))
            
            if not sys.stdin.isatty():
                raise ValueError("Multiple precursor sequences available for {}. Please specify one with --pre-id: {}".format(args.mirna_id, ', '.join(available_pre_ids)))
            
            selected_pre_id = prompt_pre_id_selection(available_pre_ids)
            mature_seq, pre_seq, start, end, _ = find_mirna_sequence(data, args.mirna_id, selected_pre_id)
        
        if not mature_seq or pre_seq is None or start is None or end is None:
            raise ValueError("No sequence found for {}".format(args.mirna_id))
        
        validate_sequence_length(mature_seq, "Mature")
        
        # Prepare sequences
        sequences = [("{},WT".format(args.mirna_id), mature_seq)]
        
        # Handle operations based on arguments
        if args.both:
            mod_shift_seq, info = handle_combined_operation(
                pre_seq, mature_seq, start, end, args.modification, args.shift
            )
            if validate_sequence_length(mod_shift_seq, "Modified-shifted"):
                sequences.append(("{},{},modified_shifted".format(args.mirna_id, info), mod_shift_seq))
        else:
            if args.modification:
                mod_seq, mods = apply_modifications(mature_seq, args.modification)
                if mods and validate_sequence_length(mod_seq, "Modified"):
                    sequences.append(("{},{},modified".format(args.mirna_id, '&'.join(mods)), mod_seq))
            
            if args.shift:
                shifted_seq, shift_info = apply_shift(pre_seq, mature_seq, start, end, args.shift)
                if validate_sequence_length(shifted_seq, "Shifted"):
                    sequences.append(("{},{},shifted".format(args.mirna_id, shift_info), shifted_seq))
        
        # Write output
        write_fasta(args.output, sequences)
        print("Successfully wrote {} sequences to {}".format(len(sequences), args.output))
        
    except Exception as e:
        print("Error: {}".format(str(e)), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
