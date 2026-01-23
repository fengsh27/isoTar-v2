import os, re
import sys
import argparse
import multiprocessing
import subprocess
import json
import shutil
from typing import Dict, List

# Predefined list of tools
ALLOWED_TOOLS = ["miRanda", "miRmap", "Targetscan", "RNAhybrid", "PITA", "DMISO"]
# miRanda path
MIRANDA = '/opt/miRanda/bin/miranda'
# RNAhybrid path
RNAHYBRID = "/opt/RNAhybrid/src/RNAhybrid"
# PITA path
PITA = "/opt/PITA64bit/pita_prediction.pl"
# TargetScan path
TARGETSCAN = "/opt/TargetScan/"
# DMISO path
DMISO = "/usr/local/bin/dmiso"

# 3 UTR PATH
HUMAN_HG19_3UTR = '/opt/human/hg19/3utr.fa'
HUMAN_HG38_3UTR = '/opt/human/hg38/3utr.fasta'

# Ensure miRmap modules can be imported
MIRMAP_SRC = "/opt/miRmap/src"
if MIRMAP_SRC not in sys.path:
    sys.path.append(MIRMAP_SRC)

def parse_fasta(fasta_file: str) -> List[Dict[str, str]]:
    """
    Parse a FASTA file and extract headers, sequences, their lengths, and types.
    Remove sequences with lengths smaller than 17 or greater than 30.
    """
    sequences = []
    with open(fasta_file, 'r') as file:
        header = ""
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if header:  # Save the previous sequence
                    # Determine the type based on the header
                    if "WT" in header:
                        seq_type = "WT"
                    elif "modified" in header:
                        seq_type = "modified"
                    elif "shifted" in header:
                        seq_type = "shifted"
                    elif "modified_shifted" in header:
                        seq_type = "modified_shifted"
                    else:
                        seq_type = "unknown"  # Default type if no match
                    
                    # Check sequence length
                    if 17 <= len(sequence) <= 30:
                        sequences.append({
                            "header": header,
                            "sequence": sequence,
                            "length": len(sequence),
                            "type": seq_type
                        })
                    else:
                        print(f"Removed sequence '{header}' (type: {seq_type}) with length {len(sequence)}.")
                
                header = line[1:]  # Remove ">" from the header
                sequence = ""  # Reset sequence
            else:  # Sequence line
                sequence += line

        # Save the last sequence
        if header:
            # Determine the type based on the header
            if "WT" in header:
                seq_type = "WT"
            elif "modified" in header:
                seq_type = "modified"
            elif "shifted" in header:
                seq_type = "shifted"
            elif "modified_shifted" in header:
                seq_type = "modified_shifted"
            else:
                seq_type = "unknown"  # Default type if no match
            
            # Check sequence length
            if 17 <= len(sequence) <= 30:
                sequences.append({
                    "header": header,
                    "sequence": sequence,
                    "length": len(sequence),
                    "type": seq_type
                })
            else:
                print(f"Removed sequence '{header}' (type: {seq_type}) with length {len(sequence)}.")
    return sequences

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

def save_to_json(data: List[Dict[str, str]], tools: List[str], num_cores: int, output_path: str):
    """
    Save the parsed FASTA data and tools list to a JSON file.
    """
    output_data = {
        "sequences": data,
        "tools": tools,
        "num_cores": num_cores
    }
    with open(output_path, 'w') as file:
        json.dump(output_data, file, indent=4)
    print(f"Input FASTA data and tools saved to: {output_path}")

def process_3utr_fasta(utr_file: str, num_cores: int, temp_folder: str):
    """
    Process a 3' UTR FASTA file by splitting it into multiple subfiles based on the number of cores.
    Count sequences by the number of lines starting with ">".
    """
    print(f"Processing 3' UTR FASTA file: {utr_file}")
    
    # Create temp directory
    os.makedirs(temp_folder+"/utr", exist_ok=True)
    
    # Count the number of sequences by counting lines starting with ">"
    total_sequences = 0
    with open(utr_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                total_sequences += 1
    print(f"Total sequences in 3' UTR file: {total_sequences}")
    
    # Calculate the number of sequences per subfile
    sequences_per_file = total_sequences // num_cores
    remainder = total_sequences % num_cores
    print(f"Splitting into {num_cores} subfiles with {sequences_per_file} sequences each (remainder: {remainder})")
    
    # Split sequences into subfiles
    current_sequence_count = 0
    current_file_index = 1
    subfile_name = f"{temp_folder}/utr/temp_3utr_part{current_file_index}.fasta"
    subfile = open(subfile_name, 'w')
    print(f"Subfile {current_file_index} created: {subfile_name}")
    
    with open(utr_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                current_sequence_count += 1
                # Check if we need to switch to a new subfile
                if current_sequence_count > sequences_per_file + (1 if current_file_index <= remainder else 0):
                    subfile.close()
                    current_file_index += 1
                    subfile_name = f"{temp_folder}/utr/temp_3utr_part{current_file_index}.fasta"
                    subfile = open(subfile_name, 'w')
                    print(f"Subfile {current_file_index} created: {subfile_name}")
                    current_sequence_count = 1
            subfile.write(line)
    
    subfile.close()
    print("3' UTR FASTA file splitting complete.")

def sanitize_output_path(output_file: str) -> str:
    """Sanitize output filename to avoid illegal characters."""
    output_dir, filename = os.path.split(output_file)
    safe_filename = re.sub(r"[^A-Za-z0-9._-]+", "_", filename)
    if safe_filename in {"", ".", ".."}:
        safe_filename = "miranda_output.txt"
    return os.path.join(output_dir, safe_filename)


def run_miranda(mirna_file: str, utr_file: str, output_file: str):
    """Run miRanda on a given miRNA and UTR file."""
    safe_output_file = sanitize_output_path(output_file)
    cmd = [
        MIRANDA,
        mirna_file,
        utr_file,
        "-en", "-20",
        "-out", safe_output_file,
        "-quiet"
    ]
    subprocess.run(cmd, check=True)

def run_rnahybrid(mirna_file: str, utr_file: str, output_file: str, mirna_length: int, utr_length: int):
    """Run RNAhybrid on a given miRNA and UTR file."""
    cmd = [
        RNAHYBRID,
        "-c",
        "-e", "-20",
        "-s", "3utr_human",
        "-q", mirna_file,
        "-t", utr_file,
        "-m", str(utr_length),
        "-n", str(mirna_length)
    ]
    with open(output_file, 'w') as outfile:
        subprocess.run(cmd, check=True, stdout=outfile)

def run_pita(mirna_file: str, utr_file: str, output_prefix: str):
    """Run PITA on a given miRNA and UTR file."""
    cmd = [
        "perl",
        PITA,
        "-mir", mirna_file,
        "-utr", utr_file,
        "-prefix", output_prefix,
        "-l", "7-8",
        "-gxp",
        "-gu", "7;0;8;0",
        "-m", "7;0;8;0"
    ]
    # Run PITA and redirect output to /dev/null
    with open('/dev/null', 'w') as devnull:
        subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)

def run_targetscan(targetscan_input: str, utr_input: str, output_file_1: str, bln_bins_file: str, output_file_2: str):
    """Run TargetScan Script 1"""
    cmd1 = ["perl", TARGETSCAN+"TargetScan_70/targetscan_70.pl", targetscan_input, utr_input, output_file_1]
    subprocess.run(cmd1, check=True)
    
    """Run TargetScan Script 2"""
    cmd2 = ['perl', TARGETSCAN+"TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl", targetscan_input, output_file_1, bln_bins_file]
    with open(output_file_2, 'w') as outfile, open(os.devnull, 'w') as devnull:
        subprocess.run(cmd2, stdout=outfile, stderr=devnull, check=True)
        
def run_dmiso(mirna_file: str, utr_file: str, output_file: str):
    """Run DMISO on a given miRNA and UTR file."""
    cmd = [
        "python3",
        DMISO,
        "-m",mirna_file,
        "-t",utr_file,
        "-o", output_file
    ]
    # Run DMISO and redirect output to /dev/null
    with open('/dev/null', 'w') as devnull:
        subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)

def targetscan_prep(sequence: str, header: str, out_dir: str):
    """TargetScan_prep"""
    # load mirR_Family_Info
    mirna_family_info_path = '/opt/TargetScan/Datasets/miR_Family_Info.json'
    mirna_family_info = load_json(mirna_family_info_path)
    # Prepare TargetScan miRNA file
    mirna_fasta_path = f"{out_dir}/{header}_targetscan.txt"
    with open(mirna_fasta_path, 'w') as f:
        # seed region
        seed = sequence[1:8]        
        # Check U-T nucleotides
        mirna_u = seed.find('U')        
        # Replace miRNA's T with U
        if mirna_u == -1:
            seed = seed.replace('T', 'U')
        seed = seed.upper()        
        # Default
        species_id = '9606'
        # If the seed exists into miR Family Info
        if seed in mirna_family_info:
            species_id = ';'.join(mirna_family_info[seed])
        # Format the output line
        identifier_clean = header.split(",")[0].replace('hsa-', '')
        line = f"{identifier_clean}\t{seed}\t{species_id}\n"  # Added newline character
        f.write(line)

def run_mirmap(mirna_seq: str, mirna_header: str, utr_file: str, output_file: str):
    """Run miRmap on a given miRNA sequence and UTR file."""
    import mirmap.target
    import mirmap.if_lib_spatt
    import mirmap.scores
    
    # Read UTR sequences from file
    utr_sequences = []
    utr_headers = []
    with open(utr_file, 'r') as f:
        current_seq = ""
        for line in f:
            if line.startswith(">"):
                utr_header = line.split(" ")[0]
                utr_header = utr_header.split("_")[2]
                utr_headers.append(utr_header)
                if current_seq:
                    utr_sequences.append(current_seq)
                    current_seq = ""
            else:
                current_seq += line.strip()
        if current_seq:
            utr_sequences.append(current_seq)
    # Check length of utr file
    if len(utr_headers) != len(utr_sequences):
        print(utr_file+": The number of utr header and utr sequence not match.")
    # Prepare output file
    with open(output_file, 'w') as out_f:
        # Process each UTR sequence
        for i in range(len(utr_headers)):
            out_f.write(">"+mirna_header+" "+utr_headers[i]+"\n\n")
            # Convert miRNA sequence (U to T)
            mirna_seq_t = mirna_seq.replace("U", "T")
                
            # Find targets
            targets = mirmap.target.find_targets_with_seed(utr_sequences[i], mirna_seq_t)
            if not targets:
                out_f.write("\n")
                continue           
            # Initialize Spatt and Calculate scores
            spatt_path = "/opt/miRmap/libs/default/libspatt2.so"  # Adjust path as needed
            if_spatt = mirmap.if_lib_spatt.Spatt(spatt_path)
            scores = mirmap.scores.calc_scores(targets[0], if_spatt=if_spatt)              
            
            # Write scores
            out_f.write(targets[0].report() + "\n")
            out_f.write(mirmap.scores.report_scores(scores) + "\n\n")

def parse_dmiso_results(dmiso_file: str, output_file: str):
    """Parse DMISO results."""
    filtered_results = []   
    with open(dmiso_file, 'r') as f:
        # Read header line to get column indices
        header = f.readline().strip().split('\t')
        
        # Find column indices
        target_id_idx = header.index('Target ID')
        target_seq_idx = header.index('Target Sequence')
        pred_score_idx = header.index('Prediction Score')
        
        # Process each data line
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            columns = line.split('\t')
            
            # Check if we have enough columns
            if len(columns) <= max(target_id_idx, target_seq_idx, pred_score_idx):
                continue
            
            try:
                # Extract prediction score and filter
                pred_score = float(columns[pred_score_idx])
                if pred_score <= 0.99:
                    continue
                
                # Extract ENST ID from Target ID
                target_id_full = columns[target_id_idx]
                # Pattern to match ENST ID: ENST followed by digits and optional .version
                enst_match = re.search(r'(ENST\d+\.?\d*)', target_id_full)
                if enst_match:
                    enst_id = enst_match.group(1)
                else:
                    # If no ENST pattern found, use the original target ID
                    enst_id = target_id_full
                
                target_sequence = columns[target_seq_idx]
                
                filtered_results.append([enst_id, target_sequence, str(pred_score)])
                
            except (ValueError, IndexError) as e:
                print(f"Error processing line: {line}. {e}")
                continue
    
    # Write filtered results to TSV file
    with open(output_file, 'w') as out_f:
        # Write header
        out_f.write("Target ID\tTarget Sequence\tPrediction Score\n")
        
        # Write filtered data
        for result in filtered_results:
            out_f.write('\t'.join(result) + '\n')
    
    # Remove the original DMISO file
    os.remove(dmiso_file)
                
def get_longest_utr_length(utr_file: str) -> int:
    """Get the length of the longest sequence in a UTR file."""
    max_length = 0
    current_length = 0
    with open(utr_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_length > max_length:
                    max_length = current_length
                current_length = 0
            else:
                current_length += len(line.strip())
    return max(max_length, current_length)
    
def process_tools(sequences: List[Dict[str, str]], tools: List[str], utr_file: str, output_folder: str, temp_folder: str):
    seq_num = 0
    for seq in sequences:
        # Create a temporary FASTA file for the single sequence
        temp_fasta = f"{temp_folder}/seq_{seq_num}.fasta"
        name_fasta = f"{temp_folder}/{seq['sequence']}.fasta"
        with open(temp_fasta, 'w') as file:
            file.write(f">{seq['header']}\n{seq['sequence']}\n")
       
        for tool in tools:
            if tool == "miRanda":
                # Output directory for miRanda
                miranda_out_dir = os.path.join(output_folder, "miRanda")
                # Define the output file path
                output_file = f"{miranda_out_dir}/{seq['header']}_miRanda_results.txt"
                # Run the tool
                print(f"miRanda is processing {name_fasta}")
                run_miranda(temp_fasta, utr_file, output_file)
            elif tool == "RNAhybrid":
                # Output directory for miRanda
                rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
                # Define the output file path
                output_file = f"{rnahybrid_out_dir}/{seq['header']}_RNAhybrid_results.txt"
                # Prepare parameters
                max_utr_length = get_longest_utr_length(utr_file)
                # Run the tool
                print(f"RNAhybrid is processing {name_fasta}")
                run_rnahybrid(temp_fasta, utr_file, output_file, seq['length'], max_utr_length)
            elif tool == "miRmap":
                # Create output directory for miRmap
                mirmap_out_dir = os.path.join(output_folder, "miRmap")
                # Define the output file path
                output_file = f"{mirmap_out_dir}/{seq['header']}_miRmap_results.txt"
                # Run the tool
                print(f"miRmap is processing {name_fasta}")
                run_mirmap(seq['sequence'], seq['header'], utr_file, output_file) 
            elif tool == "DMISO":         
                # Output directory for DMISO
                dmiso_out_dir = os.path.join(output_folder, "DMISO")
                # Define the output file path
                temp_output_file = f"{dmiso_out_dir}/{seq['header']}_DMISO.txt"
                output_file = f"{dmiso_out_dir}/{seq['header']}_DMISO_results.txt"
                # Run the tool
                print(f"DMISO is processing {name_fasta}")
                run_dmiso(temp_fasta, utr_file, temp_output_file)
                # Parse DMISO results
                parse_dmiso_results(temp_output_file, output_file)
            elif tool == "PITA":
                # Output directory for PITA
                pita_out_dir = os.path.join(output_folder, "PITA")
                # Define the output file path
                output_file_prefix = f"{pita_out_dir}/{seq['header']}"       
                # Run the tool
                print(f"PITA is processing {name_fasta}")
                run_pita(temp_fasta, utr_file, output_file_prefix)
                # remove temp file
                if os.path.exists("tmp_seqfile1"):
                    os.remove("tmp_seqfile1")
                if os.path.exists("tmp_seqfile2"):
                    os.remove("tmp_seqfile2")
                if os.path.exists(output_file_prefix+"_pita_results.gxp"):
                    os.remove(output_file_prefix+"_pita_results.gxp")                
            elif tool == "Targetscan":
                                # Output directory for PITA
                targetscan_out_dir = os.path.join(output_folder, "Targetscan")
                # Define the output file path
                output_file1 = f"{targetscan_out_dir}/{seq['header']}_Targetscan_results1.txt"
                output_file2 = f"{targetscan_out_dir}/{seq['header']}_Targetscan_results2.txt"       
                # Run the tool
                print(f"Targetscan is processing {name_fasta}")
                targetscan_prep(seq['sequence'], seq['header'], targetscan_out_dir)
                # TargetScan Input File path
                targetscan_input = f"{targetscan_out_dir}/{seq['header']}_targetscan.txt"
                # utr path
                utr_path = "/opt/TargetScan/Datasets/3utr"
                bln_bins_path = "/opt/TargetScan/Datasets/bln_bins"
                # Process Targetscan
                for i in range(64):
                    utr_file = os.path.join(utr_path, f'targetscan_utr_part_{i}.txt')
                    output_file_1 = f"{targetscan_out_dir}/{seq['header']}_part_{i}_out1.txt"
                    bln_bins_file = os.path.join(bln_bins_path, f'targetscan_median_bls_bins_part_{i}.txt')
                    output_file_2 = f"{targetscan_out_dir}/{seq['header']}_part_{i}_out2.txt"
                    # Run targetscan
                    run_targetscan(targetscan_input, utr_file, output_file_1, bln_bins_file, output_file_2)
                # Merge results for first output file
                with open(output_file1, 'w') as merged:
                    # Write header from first file (assuming all have same header)
                    first_file = f"{targetscan_out_dir}/{seq['header']}_part_0_out1.txt"
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged.write(header)
                    
                    # Apeend content from all files
                    for i in range(64):
                        part_file = f"{targetscan_out_dir}/{seq['header']}_part_{i}_out1.txt"
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)
                
                # Merge results for second output file
                with open(output_file2, 'w') as merged:
                    # Write header from first file (assuming all have same header)
                    first_file = f"{targetscan_out_dir}/{seq['header']}_part_0_out2.txt"
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged.write(header)
                    
                    # Apeend content from all files
                    for i in range(64):
                        part_file = f"{targetscan_out_dir}/{seq['header']}_part_{i}_out2.txt"
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)               
            else:
                # Handle other tools
                print(f"Tool {tool} is processing {name_fasta}")                              
        seq_num += 1
        
def process_tools_in_parallel(sequences: List[Dict[str, str]], tools: List[str], num_cores: int, output_folder: str, temp_folder: str):
    # Get all UTR subfiles
    utr_subfiles = [os.path.join(temp_folder+"/utr", f) for f in os.listdir(temp_folder+"/utr") if f.startswith("temp_3utr_part")]
    
    # Create a pool of workers
    pool = multiprocessing.Pool(processes=num_cores)
    seq_num = 0
    for seq in sequences:
        # Create a temporary FASTA file for the single sequence
        temp_fasta = f"{temp_folder}/seq_{seq_num}.fasta"
        name_fasta = f"{temp_folder}/{seq['sequence']}.fasta"
        with open(temp_fasta, 'w') as file:
            file.write(f">{seq['header']}\n{seq['sequence']}\n")

        for tool in tools:
            if tool == "miRanda":
                # Create output directory for miRanda
                miranda_out_dir = os.path.join(output_folder, "miRanda")
                # Define the output file path
                output_file = f"{miranda_out_dir}/{seq['header']}_miRanda_results.txt"
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(miranda_out_dir, f"Seq_{seq_num}_miRanda_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print(f"miRanda is processing {name_fasta}")
                pool.starmap(run_miranda, args)
                
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(miranda_out_dir, f"Seq_{seq_num}_miRanda_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
            
            elif tool == "RNAhybrid":
                # Create output directory for miRanda
                rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
                # Define the output file path
                output_file = f"{rnahybrid_out_dir}/{seq['header']}_RNAhybrid_results.txt"
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(rnahybrid_out_dir, f"Seq_{seq_num}_RNAhybrid_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                    max_utr_length = get_longest_utr_length(utr_file)
                    args.append((temp_fasta, utr_file, temp_output_file, seq['length'], max_utr_length))
                
                # Run in parallel
                print(f"RNAhybrid is processing {name_fasta}")
                pool.starmap(run_rnahybrid, args)
                
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(rnahybrid_out_dir, f"Seq_{seq_num}_RNAhybrid_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                            
            elif tool == "miRmap":
                # Create output directory for miRmap
                mirmap_out_dir = os.path.join(output_folder, "miRmap")
                # Define the output file path
                output_file = f"{mirmap_out_dir}/{seq['header']}_miRmap_results.txt"
    
                # Prepare arguments for parallel processing
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(mirmap_out_dir, f"Seq_{seq_num}_miRmap_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                    args.append((seq['sequence'], seq['header'], utr_file, temp_output_file))
    
                # Run in parallel
                print(f"miRmap is processing {name_fasta}")
                pool.starmap(run_mirmap, args)
    
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(mirmap_out_dir, f"Seq_{seq_num}_miRmap_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                
            elif tool == "DMISO":
                # Create output directory for DMISO
                dmiso_out_dir = os.path.join(output_folder, "DMISO")
                # Define the output file path
                output_file_before = f"{dmiso_out_dir}/{seq['header']}_DMISO.txt"
                output_file = f"{dmiso_out_dir}/{seq['header']}_DMISO_results.txt"
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(dmiso_out_dir, f"Seq_{seq_num}_DMISO_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print(f"DMISO is processing {name_fasta}")
                pool.starmap(run_dmiso, args)
                
                # Merge results
                with open(output_file_before, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(dmiso_out_dir, f"Seq_{seq_num}_DMISO_results_{os.path.basename(utr_file).replace('.fasta', '')}.out")
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                # Parse DMISO results
                parse_dmiso_results(output_file_before, output_file)
                                                     
            elif tool == "PITA":
                # Create output directory for PITA
                pita_out_dir = os.path.join(output_folder, "PITA")
                # Define the output file path
                output_file_prefix = f"{pita_out_dir}/{seq['header']}_PITA_results"
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(pita_out_dir, f"Seq_{seq_num}_{os.path.basename(utr_file).replace('.fasta', '')}")
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print(f"PITA is processing {name_fasta}")
                pool.starmap(run_pita, args)
                
                # Merge results
                with open(output_file_prefix+".tab", 'w') as merged_file:
                    first_file = f"{pita_out_dir}/Seq_{seq_num}_temp_3utr_part1_pita_results.tab"
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged_file.write(header)
                    # Apeend content from all files
                    for i in range(num_cores):
                        part_file = f"{pita_out_dir}/Seq_{seq_num}_temp_3utr_part{i+1}_pita_results.tab"
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged_file.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)
                        part_file = f"{pita_out_dir}/Seq_{seq_num}_temp_3utr_part{i+1}_pita_results.gxp"
                        if os.path.exists(part_file):
                            os.remove(part_file)
       
                # Merge results
                with open(output_file_prefix+"_targets.tab", 'w') as merged_file:
                    first_file = f"{pita_out_dir}/Seq_{seq_num}_temp_3utr_part1_pita_results_targets.tab"
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged_file.write(header)
                    # Apeend content from all files
                    for i in range(num_cores):
                        part_file = f"{pita_out_dir}/Seq_{seq_num}_temp_3utr_part{i+1}_pita_results_targets.tab"
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged_file.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)
                
                # remove temp file
                if os.path.exists("tmp_seqfile1"):
                    os.remove("tmp_seqfile1")
                if os.path.exists("tmp_seqfile2"):
                    os.remove("tmp_seqfile2")           
            elif tool == "Targetscan":         
                # Run the tool
                print(f"Targetscan is processing {name_fasta}")                
            else:
                # Handle other tools
                print(f"Tool {tool} is processing {name_fasta}")                
    seq_num += 1
    pool.close()
    pool.join()

def cleanup_temp_folder(temp_folder):
    """
    Delete the temp folder and its contents after processing is complete.
    """
    if os.path.exists(temp_folder):
        shutil.rmtree(temp_folder)
        print("Temp folder and its contents deleted.")
    else:
        print("Temp folder does not exist.")

def main():
    parser = argparse.ArgumentParser(
        description="MicroRNA Sequence Prediction Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-c", "--cores", type=int, required=True, help="number of CPU cores")
    parser.add_argument("-i", "--input", type=str, required=True, help="path to mirna sequence input file")
    parser.add_argument("-t", "--tools", type=str, nargs="+", required=True, 
                        help=f"List of tools to run. Choose from: {', '.join(ALLOWED_TOOLS)}")
    parser.add_argument("-g", "--genome", type=str, choices=["hg19", "hg38"], help="Reference genome (hg19 or hg38)", required=True)
    parser.add_argument("-o", "--output", type=str, required=True, help="output folder name")

    args = parser.parse_args()
    # Parse arguments
    num_cores = args.cores
    mirna = args.input
    tools = args.tools
    output_folder = args.output
    genome = args.genome

    # Validate number of cores
    available_cores = multiprocessing.cpu_count()
    if num_cores > available_cores:
        print(f"Warning: Requested {num_cores} cores, but only {available_cores} are available.")
        num_cores = available_cores

    # Validate tools
    invalid_tools = [tool for tool in tools if tool not in ALLOWED_TOOLS]
    if invalid_tools:
        raise ValueError(f"Invalid tools selected: {', '.join(invalid_tools)}. Allowed tools are: {', '.join(ALLOWED_TOOLS)}")

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    temp_folder = output_folder+"/temp"
    os.makedirs(temp_folder, exist_ok=True)

    # Parse FASTA files
    sequences = parse_fasta(mirna)
    save_to_json(sequences, tools, num_cores, f"{output_folder}/mirna_prediction_parameters.json")  # Save JSON to output folder
    
    # Determine the 3' UTR file based on the genome
    utr_file = HUMAN_HG19_3UTR if genome == "hg19" else HUMAN_HG38_3UTR

    # Create output directory
    for tool in tools:
        if tool == "miRanda":
            miranda_out_dir = os.path.join(output_folder, "miRanda")
            os.makedirs(miranda_out_dir, exist_ok=True)
        elif tool == "RNAhybrid":
            rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
            os.makedirs(rnahybrid_out_dir, exist_ok=True)
        elif tool == "miRmap":
            mirmap_out_dir = os.path.join(output_folder, "miRmap")
            os.makedirs(mirmap_out_dir, exist_ok=True)
        elif tool == "DMISO":
            dmiso_out_dir = os.path.join(output_folder, "DMISO")
            os.makedirs(dmiso_out_dir, exist_ok=True)
        elif tool == "PITA":
            pita_out_dir = os.path.join(output_folder, "PITA")
            os.makedirs(pita_out_dir, exist_ok=True)
        elif tool == "Targetscan":
            targetscan_out_dir = os.path.join(output_folder, "Targetscan")
            os.makedirs(targetscan_out_dir, exist_ok=True)
        else:
            other_out_dir = os.path.join(output_folder, "Other")
            os.makedirs(other_out_dir, exist_ok=True)
              
    # Run Prediction for single or mutiple cores
    if num_cores == 1:
        process_tools(sequences, tools, utr_file, output_folder, temp_folder)
    else:
        process_3utr_fasta(utr_file, num_cores, temp_folder)
        process_tools_in_parallel(sequences, tools, num_cores, output_folder ,temp_folder)

    # Clean up the temp folder
    cleanup_temp_folder(temp_folder)

    print("Processing complete.")

if __name__ == "__main__":
    main()
