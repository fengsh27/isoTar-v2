#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import json
import csv
import sys

# Matches Ensembl (ENST00000284637) or RefSeq mRNA (NM_001164664), stripping version suffix
_TRANSCRIPT_RE = re.compile(r'(ENST\d+|NM_\d+)(?:\.\d+)?')

def _extract_transcript_id(text):
    """Return ENST or NM_ transcript ID from text, without version suffix. None if not found."""
    m = _TRANSCRIPT_RE.search(text)
    return m.group(1) if m else None


def read_sequences_from_json(json_file):
    # Read sequences from a JSON file.

    if not os.path.exists(json_file):
        raise FileNotFoundError("JSON file not found: {}".format(json_file))
        
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            
        # Handle different possible JSON structures
        if isinstance(data, dict):
            if 'sequences' in data:
                sequences = data['sequences']
            else:
                sequences = [data]
        elif isinstance(data, list):
            sequences = data
        else:
            raise ValueError("Invalid JSON structure. Expected a dictionary with 'sequences' key or a list of sequences")
            
        if not isinstance(sequences, list):
            raise ValueError("Sequences must be a list")
            
        # Validate each sequence
        valid_sequences = []
        for seq in sequences:
            valid_sequences.append(seq)
            
        if not valid_sequences:
            raise ValueError("No valid sequences found in JSON file")
            
        return valid_sequences
        
    except json.JSONDecodeError as e:
        raise ValueError("Invalid JSON file: {}".format(str(e)))
    except Exception as e:
        raise Exception("Error reading JSON file: {}".format(str(e)))

def parseTargetScanResults(output_f_path, result_dict):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:        
            handler = csv.reader(f, delimiter='\t')
            next(handler) # Skip header row
            for line in handler:  
                if len(line) == 14:
                    if line[2] =="9606" and line[8] != '6mer':
                        # Target - remove version number if present
                        tar = re.sub(r'\.[0-9]+$', '', line[0])
                        if tar and tar not in results:
                            results.append(tar)      
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['TargetScan'] = results

    return result_dict

def parsePITAResults(output_f_path, result_dict):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:        
            handler = csv.reader(f, delimiter='\t')
            for line in handler:
                if len(line) == 13:
                    tar = _extract_transcript_id(line[0])
                    if tar:
                        ddG = float(line[12])
                        if ddG <= -10.0:
                            if tar not in results:
                                results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['PITA'] = results

    return result_dict

def parseRnahybridResults(output_f_path, result_dict):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:        
            handler = csv.reader(f, delimiter=':')
            for line in handler:
                if len(line) == 11:
                    tar = _extract_transcript_id(line[0])
                    if tar:
                        # Get the seed region 2-7
                        target_seq = line[8][-8:-1]
                        mirna_seq = line[9][-8:-1]
                        if not re.search(r'\s', mirna_seq):
                            seq_check = False
                            for i in range(0, len(mirna_seq)):  # Changed xrange to range
                                # miRNA nucleotide
                                mir_nt = mirna_seq[i]
                                # Target nucleotide
                                tar_nt = target_seq[i]
                                if (mir_nt.upper() == 'U' and tar_nt.upper() == 'G') or \
                                    (mir_nt.upper() == 'G' and tar_nt.upper() == 'U'):
                                    seq_check = True
                                    break
                                elif (mir_nt.upper() == 'T' and tar_nt.upper() == 'G') or \
                                    (mir_nt.upper() == 'G' and tar_nt.upper() == 'T'):
                                    seq_check = True
                                    break
                            if not seq_check and tar not in results:
                                results.append(tar)          
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['RNAhybrid'] = results

    return result_dict

def parseMirmapResults(output_f_path, result_dict):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                # miRNA - Target
                matchObj = re.match(r'^>[^,]+,.*?\s+(\S+)\s*$', lines[i], re.M|re.I)
                if matchObj:
                    tar = _extract_transcript_id(matchObj.group(1))
                    i += 2
                    if i < len(lines):
                        matchObj = re.match(r'.*[0-9]+.*', lines[i], re.M|re.I) # check this line contain any number
                        if matchObj:
                            i += 6
                            if i < len(lines):
                                matchObj = re.match(r'^\s*ΔG binding \(kcal/mol\)\s+([-+]?\d+\.?\d*)\s*$', lines[i])
                                if matchObj:
                                    dg_binding = float(matchObj.group(1))
                                    # Only add if ΔG binding < -20
                                    if dg_binding is not None and dg_binding <= -20:
                                        if tar not in results:
                                            results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['miRmap'] = results

    return result_dict

def parseMirandaResults(output_f_path, result_dict):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                # Query
                matchObj = re.match(r'^\s+Query:\s+3\'\s+([^\s]+)\s+5\'$', lines[i], re.M|re.I)
                if matchObj:
                    i += 1
                    # Pairing
                    if i < len(lines):
                        # Get the seed region 2-7
                        seed_region = lines[i][-9:-2]
                        if len(seed_region.replace('|', '')) == 0:
                            i += 1
                            if i < len(lines):
                                matchObj1 = re.match(r'^\s+Ref:\s+5\'\s+([^\s]+)\s+3\'$', lines[i], re.M|re.I)
                                i += 5
                                # miRNA - Target
                                if i < len(lines):
                                    matchObj2 = re.match(r'^>([^\s]+)\s+(\S+)', lines[i], re.M|re.I)
                                    if matchObj2:
                                        tar = _extract_transcript_id(matchObj2.group(2))
                                        if tar and tar not in results:
                                            results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['miRanda'] = results

    return result_dict

def parseDMISOResults(output_f_path, result_dict):
    return result_dict
def process_sequence(sequence, result_dir):
    # Process a single sequence and generate prediction results. 
    try:
        output_f_path_miRanda = os.path.join(result_dir, "miRanda", "{}_miRanda_results.txt".format(sequence['header']))
        output_f_path_miRmap = os.path.join(result_dir, "miRmap", "{}_miRmap_results.txt".format(sequence['header']))
        output_f_path_RNAhybrid = os.path.join(result_dir, "RNAhybrid", "{}_RNAhybrid_results.txt".format(sequence['header']))
        output_f_path_PITA = os.path.join(result_dir, "PITA", "{}_PITA_results.tab".format(sequence['header']))
        output_f_path_TargetScan = os.path.join(result_dir, "TargetScan", "{}_Targetscan_output_sort.txt".format(sequence['header']))

        prediction_results = sequence.copy()
        if os.path.exists(output_f_path_miRanda):
            prediction_results = parseMirandaResults(output_f_path_miRanda, prediction_results)
        if os.path.exists(output_f_path_miRmap):
            prediction_results = parseMirmapResults(output_f_path_miRmap, prediction_results)
        if os.path.exists(output_f_path_RNAhybrid):
            prediction_results = parseRnahybridResults(output_f_path_RNAhybrid, prediction_results)
        if os.path.exists(output_f_path_PITA):
            prediction_results = parsePITAResults(output_f_path_PITA, prediction_results)
        if os.path.exists(output_f_path_TargetScan):
            prediction_results = parseTargetScanResults(output_f_path_TargetScan, prediction_results)
        
        return prediction_results
        
    except Exception as e:
        # Add more context to the error
        error_context = {
            'sequence_header': sequence.get('header', 'UNKNOWN'),
            'sequence_type': sequence.get('type', 'UNKNOWN'),
            'error_type': type(e).__name__,
            'error_message': str(e)
        }
        raise Exception("Error processing sequence: {}".format(json.dumps(error_context, indent=2)))

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Process miRNA sequences and generate prediction results')
    parser.add_argument('result_dir', help='Path to miRNA prediction results')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    args = parser.parse_args()
    result_dir = args.result_dir
    output_dir = os.path.join(result_dir, "miRNA_prediction_results")
    try:
        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Read sequences from JSON file
        print("Reading sequences from {}...".format(result_dir))
        json_file = os.path.join(result_dir, "mirna_prediction_parameters.json")
        sequences = read_sequences_from_json(json_file)
        print("Found {} valid sequences".format(len(sequences)))
        
        # Process each sequence
        successful = 0
        failed = 0
        
        for sequence in sequences:
            try:
                prediction_results = process_sequence(sequence, result_dir)
                
                # Generate output filename
                output_filename = "{}_results.json".format(sequence['header'])
                output_path = os.path.join(output_dir, output_filename)
                
                # Write results to file
                with open(output_path, 'w') as file:
                    json.dump(prediction_results, file, indent=4)

                print("Processed sequence {} - Results saved to {}".format(sequence['header'], output_path))
                successful += 1
                
            except Exception as e:
                failed += 1
                print("Error processing sequence {}:".format(sequence['header']))
                if args.verbose:
                    print("Error details:\n{}".format(str(e)))
                    print("Sequence data: {}".format(json.dumps(sequence, indent=2)))
                else:
                    print("Error: {}".format(str(e)))
                print()
        
        # Print summary
        print("\nProcessing Summary:")
        print("Total sequences: {}".format(len(sequences)))
        print("Successfully processed: {}".format(successful))
        print("Failed to process: {}".format(failed))
        
        if failed > 0:
            sys.exit(1)
            
    except Exception as e:
        print("Fatal error: {}".format(str(e)))
        sys.exit(1)

if __name__ == "__main__":
    main()
