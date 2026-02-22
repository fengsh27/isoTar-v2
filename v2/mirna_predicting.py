import os, re
import sys
import argparse
import multiprocessing
import subprocess

if not hasattr(subprocess, "run"):
    def _run(cmd, check=False, stdout=None, stderr=None):
        result = subprocess.call(cmd, stdout=stdout, stderr=stderr)
        if check and result != 0:
            raise RuntimeError("Command failed with exit status {}: {}".format(result, cmd))
        return result

    subprocess.run = _run
import json
import shutil


# Predefined list of tools
ALLOWED_TOOLS = ["miRanda", "miRmap", "Targetscan", "RNAhybrid", "PITA", "DMISO"]
# miRanda path
MIRANDA = '/usr/local/bin/miranda'
# RNAhybrid path
RNAHYBRID = "/usr/local/bin/RNAhybrid"
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

def parse_fasta(fasta_file):
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
                        print("Removed sequence '{}' (type: {}) with length {}.".format(header, seq_type, len(sequence)))
                
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
                print("Removed sequence '{}' (type: {}) with length {}.".format(header, seq_type, len(sequence)))
    return sequences

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

def save_to_json(data, tools, num_cores, output_path):
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
    print("Input FASTA data and tools saved to: {}".format(output_path))

def process_3utr_fasta(utr_file, num_cores, temp_folder):
    """
    Process a 3' UTR FASTA file by splitting it into multiple subfiles based on the number of cores.
    Count sequences by the number of lines starting with ">".
    """
    print("Processing 3' UTR FASTA file: {}".format(utr_file))
    
    # Create temp directory
    if not os.path.exists(temp_folder + "/utr"):
        os.makedirs(temp_folder + "/utr")
    
    # Count the number of sequences by counting lines starting with ">"
    total_sequences = 0
    with open(utr_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                total_sequences += 1
    print("Total sequences in 3' UTR file: {}".format(total_sequences))
    
    # Calculate the number of sequences per subfile
    sequences_per_file = total_sequences // num_cores
    remainder = total_sequences % num_cores
    print("Splitting into {} subfiles with {} sequences each (remainder: {})".format(num_cores, sequences_per_file, remainder))
    
    # Split sequences into subfiles
    current_sequence_count = 0
    current_file_index = 1
    subfile_name = "{}/utr/temp_3utr_part{}.fasta".format(temp_folder, current_file_index)
    subfile = open(subfile_name, 'w')
    print("Subfile {} created: {}".format(current_file_index, subfile_name))
    
    with open(utr_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                current_sequence_count += 1
                # Check if we need to switch to a new subfile
                if current_sequence_count > sequences_per_file + (1 if current_file_index <= remainder else 0):
                    subfile.close()
                    current_file_index += 1
                    subfile_name = "{}/utr/temp_3utr_part{}.fasta".format(temp_folder, current_file_index)
                    subfile = open(subfile_name, 'w')
                    print("Subfile {} created: {}".format(current_file_index, subfile_name))
                    current_sequence_count = 1
            subfile.write(line)
    
    subfile.close()
    print("3' UTR FASTA file splitting complete.")

def sanitize_output_path(output_file):
    """Sanitize output filename to avoid illegal characters."""
    output_dir, filename = os.path.split(output_file)
    safe_filename = re.sub(r"[^A-Za-z0-9._-]+", "_", filename)
    if safe_filename in {"", ".", ".."}:
        safe_filename = "miranda_output.txt"
    return os.path.join(output_dir, safe_filename)


def run_miranda(mirna_file, utr_file, output_file):
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

def run_rnahybrid(mirna_file, utr_file, output_file, mirna_length, utr_length):
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

def run_pita(mirna_file, utr_file, output_prefix):
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

def run_targetscan(targetscan_input, utr_input, output_file_1, bln_bins_file, output_file_2):
    """Run TargetScan Script 1"""
    cmd1 = ["perl", TARGETSCAN+"TargetScan_70/targetscan_70.pl", targetscan_input, utr_input, output_file_1]
    subprocess.run(cmd1, check=True)
    
    """Run TargetScan Script 2"""
    cmd2 = ['perl', TARGETSCAN+"TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl", targetscan_input, output_file_1, bln_bins_file]
    with open(output_file_2, 'w') as outfile, open(os.devnull, 'w') as devnull:
        subprocess.run(cmd2, stdout=outfile, stderr=devnull, check=True)
        
def run_dmiso(mirna_file, utr_file, output_file):
    """Run DMISO on a given miRNA and UTR file."""
    cmd = [
        "python3.6",
        "/opt/DMISO/DMISO-main/dmiso.py",
        "-m",mirna_file,
        "-t",utr_file,
        "-o", output_file
    ]
    # Run DMISO and redirect output to /dev/null
    with open('/dev/null', 'w') as devnull:
        subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)


def run_miranda_with_params(params):
    return run_miranda(*params)


def run_rnahybrid_with_params(params):
    return run_rnahybrid(*params)


def run_mirmap_with_params(params):
    return run_mirmap(*params)


def run_dmiso_with_params(params):
    return run_dmiso(*params)


def run_pita_with_params(params):
    return run_pita(*params)

def targetscan_prep(sequence, header, out_dir):
    """TargetScan_prep"""
    # load mirR_Family_Info
    mirna_family_info_path = '/opt/TargetScan/Datasets/miR_Family_Info.json'
    mirna_family_info = load_json(mirna_family_info_path)
    # Prepare TargetScan miRNA file
    mirna_fasta_path = "{}/{}_targetscan.txt".format(out_dir, header)
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
        line = "{}\t{}\t{}\n".format(identifier_clean, seed, species_id)
        f.write(line)

def run_mirmap(mirna_seq, mirna_header, utr_file, output_file):
    """Run miRmap on a given miRNA sequence and UTR file."""
    import mirmap
    import mirmap.library_link
    
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
                
            # Initialize miRmap object
            mm_obj = mirmap.mm(utr_sequences[i], mirna_seq_t)
            spatt_path = os.environ.get("MIRMAP_SPATT_LIB", "/opt/miRmap/libs/default/libspatt2.so")
            if os.path.exists(spatt_path):
                mm_obj.libs = mirmap.library_link.LibraryLink(os.path.dirname(spatt_path))
            else:
                print("Warning: Spatt library not found at {}. Continuing without Spatt.".format(spatt_path))

            # Find targets with seed
            mm_obj.find_potential_targets_with_seed()
            if len(mm_obj.end_sites) == 0:
                out_f.write("\n")
                continue

            # Evaluate scores (best effort)
            mm_obj.eval_tgs_au()
            mm_obj.eval_tgs_position()
            mm_obj.eval_tgs_pairing3p()
            mm_obj.eval_tgs_score()
            mm_obj.eval_dg_duplex()
            mm_obj.eval_dg_open()
            mm_obj.eval_dg_total()
            mm_obj.eval_prob_exact()
            mm_obj.eval_prob_binomial()
            mm_obj.cons_blss = [0.0] * len(mm_obj.end_sites)
            mm_obj.selec_phylops = [1.0] * len(mm_obj.end_sites)
            mm_obj.eval_score()

            # Write report
            out_f.write(mm_obj.report() + "\n\n")

def parse_dmiso_results(dmiso_file, output_file):
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
                print("Error processing line: {}. {}".format(line, e))
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
                
def get_longest_utr_length(utr_file):
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
    
def process_tools(sequences, tools, utr_file, output_folder, temp_folder):
    seq_num = 0
    for seq in sequences:
        # Create a temporary FASTA file for the single sequence
        temp_fasta = "{}/seq_{}.fasta".format(temp_folder, seq_num)
        name_fasta = "{}/{}.fasta".format(temp_folder, seq['sequence'])
        with open(temp_fasta, 'w') as file:
            file.write(">{}\n{}\n".format(seq['header'], seq['sequence']))
       
        for tool in tools:
            if tool == "miRanda":
                # Output directory for miRanda
                miranda_out_dir = os.path.join(output_folder, "miRanda")
                # Define the output file path
                output_file = "{}/{}_miRanda_results.txt".format(miranda_out_dir, seq['header'])
                # Run the tool
                print("miRanda is processing {}".format(name_fasta))
                run_miranda(temp_fasta, utr_file, output_file)
            elif tool == "RNAhybrid":
                # Output directory for miRanda
                rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
                # Define the output file path
                output_file = "{}/{}_RNAhybrid_results.txt".format(rnahybrid_out_dir, seq['header'])
                # Prepare parameters
                max_utr_length = get_longest_utr_length(utr_file)
                # Run the tool
                print("RNAhybrid is processing {}".format(name_fasta))
                run_rnahybrid(temp_fasta, utr_file, output_file, int(seq['length']), max_utr_length)
            elif tool == "miRmap":
                # Create output directory for miRmap
                mirmap_out_dir = os.path.join(output_folder, "miRmap")
                # Define the output file path
                output_file = "{}/{}_miRmap_results.txt".format(mirmap_out_dir, seq['header'])
                # Run the tool
                print("miRmap is processing {}".format(name_fasta))
                run_mirmap(seq['sequence'], seq['header'], utr_file, output_file) 
            elif tool == "DMISO":         
                # Output directory for DMISO
                dmiso_out_dir = os.path.join(output_folder, "DMISO")
                # Define the output file path
                temp_output_file = "{}/{}_DMISO.txt".format(dmiso_out_dir, seq['header'])
                output_file = "{}/{}_DMISO_results.txt".format(dmiso_out_dir, seq['header'])
                # Run the tool
                print("DMISO is processing {}".format(name_fasta))
                run_dmiso(temp_fasta, utr_file, temp_output_file)
                # Parse DMISO results
                parse_dmiso_results(temp_output_file, output_file)
            elif tool == "PITA":
                # Output directory for PITA
                pita_out_dir = os.path.join(output_folder, "PITA")
                # Define the output file path
                output_file_prefix = "{}/{}".format(pita_out_dir, seq['header'])
                # Run the tool
                print("PITA is processing {}".format(name_fasta))
                run_pita(temp_fasta, utr_file, output_file_prefix)
                # remove temp file
                if os.path.exists("tmp_seqfile1"):
                    os.remove("tmp_seqfile1")
                if os.path.exists("tmp_seqfile2"):
                    os.remove("tmp_seqfile2")
                if os.path.exists(output_file_prefix + "_pita_results.gxp"):
                    os.remove(output_file_prefix + "_pita_results.gxp")
            elif tool == "Targetscan":
                                # Output directory for PITA
                targetscan_out_dir = os.path.join(output_folder, "Targetscan")
                # Define the output file path
                output_file1 = "{}/{}_Targetscan_results1.txt".format(targetscan_out_dir, seq['header'])
                output_file2 = "{}/{}_Targetscan_results2.txt".format(targetscan_out_dir, seq['header'])
                # Run the tool
                print("Targetscan is processing {}".format(name_fasta))
                targetscan_prep(seq['sequence'], seq['header'], targetscan_out_dir)
                # TargetScan Input File path
                targetscan_input = "{}/{}_targetscan.txt".format(targetscan_out_dir, seq['header'])
                # utr path
                utr_path = "/opt/TargetScan/Datasets/3utr"
                bln_bins_path = "/opt/TargetScan/Datasets/bln_bins"
                # Process Targetscan
                for i in range(64):
                    utr_file = os.path.join(utr_path, 'targetscan_utr_part_{}.txt'.format(i))
                    output_file_1 = "{}/{}_part_{}_out1.txt".format(targetscan_out_dir, seq['header'], i)
                    bln_bins_file = os.path.join(bln_bins_path, 'targetscan_median_bls_bins_part_{}.txt'.format(i))
                    output_file_2 = "{}/{}_part_{}_out2.txt".format(targetscan_out_dir, seq['header'], i)
                    # Run targetscan
                    run_targetscan(targetscan_input, utr_file, output_file_1, bln_bins_file, output_file_2)
                # Merge results for first output file
                with open(output_file1, 'w') as merged:
                    # Write header from first file (assuming all have same header)
                    first_file = "{}/{}_part_0_out1.txt".format(targetscan_out_dir, seq['header'])
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged.write(header)
                    
                    # Apeend content from all files
                    for i in range(64):
                        part_file = "{}/{}_part_{}_out1.txt".format(targetscan_out_dir, seq['header'], i)
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
                    first_file = "{}/{}_part_0_out2.txt".format(targetscan_out_dir, seq['header'])
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged.write(header)
                    
                    # Apeend content from all files
                    for i in range(64):
                        part_file = "{}/{}_part_{}_out2.txt".format(targetscan_out_dir, seq['header'], i)
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)               
            else:
                # Handle other tools
                print("Tool {} is processing {}".format(tool, name_fasta))
        seq_num += 1
        
def process_tools_in_parallel(sequences, tools, num_cores, output_folder, temp_folder):
    # Get all UTR subfiles
    utr_subfiles = [os.path.join(temp_folder+"/utr", f) for f in os.listdir(temp_folder+"/utr") if f.startswith("temp_3utr_part")]
    
    # Create a pool of workers
    pool = multiprocessing.Pool(processes=num_cores)
    seq_num = 0
    for seq in sequences:
        # Create a temporary FASTA file for the single sequence
        temp_fasta = "{}/seq_{}.fasta".format(temp_folder, seq_num)
        name_fasta = "{}/{}.fasta".format(temp_folder, seq['sequence'])
        with open(temp_fasta, 'w') as file:
            file.write(">{}\n{}\n".format(seq['header'], seq['sequence']))

        for tool in tools:
            if tool == "miRanda":
                # Create output directory for miRanda
                miranda_out_dir = os.path.join(output_folder, "miRanda")
                # Define the output file path
                output_file = "{}/{}_miRanda_results.txt".format(miranda_out_dir, seq['header'])
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(miranda_out_dir, "Seq_{}_miRanda_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print("miRanda is processing {}".format(name_fasta))
                pool.map(run_miranda_with_params, args)
                
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(miranda_out_dir, "Seq_{}_miRanda_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
            
            elif tool == "RNAhybrid":
                # Create output directory for miRanda
                rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
                # Define the output file path
                output_file = "{}/{}_RNAhybrid_results.txt".format(rnahybrid_out_dir, seq['header'])
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(rnahybrid_out_dir, "Seq_{}_RNAhybrid_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                    max_utr_length = get_longest_utr_length(utr_file)
                    args.append((temp_fasta, utr_file, temp_output_file, seq['length'], max_utr_length))
                
                # Run in parallel
                print("RNAhybrid is processing {}".format(name_fasta))
                pool.map(run_rnahybrid_with_params, args)
                
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(rnahybrid_out_dir, "Seq_{}_RNAhybrid_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                            
            elif tool == "miRmap":
                # Create output directory for miRmap
                mirmap_out_dir = os.path.join(output_folder, "miRmap")
                # Define the output file path
                output_file = "{}/{}_miRmap_results.txt".format(mirmap_out_dir, seq['header'])
    
                # Prepare arguments for parallel processing
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(mirmap_out_dir, "Seq_{}_miRmap_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                    args.append((seq['sequence'], seq['header'], utr_file, temp_output_file))
    
                # Run in parallel
                print("miRmap is processing {}".format(name_fasta))
                pool.map(run_mirmap_with_params, args)
    
                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(mirmap_out_dir, "Seq_{}_miRmap_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                
            elif tool == "DMISO":
                # Create output directory for DMISO
                dmiso_out_dir = os.path.join(output_folder, "DMISO")
                # Define the output file path
                output_file_before = "{}/{}_DMISO.txt".format(dmiso_out_dir, seq['header'])
                output_file = "{}/{}_DMISO_results.txt".format(dmiso_out_dir, seq['header'])
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(dmiso_out_dir, "Seq_{}_DMISO_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print("DMISO is processing {}".format(name_fasta))
                pool.map(run_dmiso_with_params, args)
                
                # Merge results
                with open(output_file_before, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(dmiso_out_dir, "Seq_{}_DMISO_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
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
                output_file_prefix = "{}/{}_PITA_results".format(pita_out_dir, seq['header'])
                # Prepare arguments for each parallel run
                args = []
                for utr_file in utr_subfiles:
                    temp_output_file = os.path.join(pita_out_dir, "Seq_{}_{}".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                    args.append((temp_fasta, utr_file, temp_output_file))
                
                # Run in parallel
                print("PITA is processing {}".format(name_fasta))
                pool.map(run_pita_with_params, args)
                
                # Merge results
                with open(output_file_prefix+".tab", 'w') as merged_file:
                    first_file = "{}/Seq_{}_temp_3utr_part1_pita_results.tab".format(pita_out_dir, seq_num)
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged_file.write(header)
                    # Apeend content from all files
                    for i in range(num_cores):
                        part_file = "{}/Seq_{}_temp_3utr_part{}_pita_results.tab".format(pita_out_dir, seq_num, i + 1)
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                # Skip header for all files
                                next(pf)
                                merged_file.write(pf.read())
                            # Remove the part file after merging
                            os.remove(part_file)
                        part_file = "{}/Seq_{}_temp_3utr_part{}_pita_results.gxp".format(pita_out_dir, seq_num, i + 1)
                        if os.path.exists(part_file):
                            os.remove(part_file)
       
                # Merge results
                with open(output_file_prefix+"_targets.tab", 'w') as merged_file:
                    first_file = "{}/Seq_{}_temp_3utr_part1_pita_results_targets.tab".format(pita_out_dir, seq_num)
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            header = first.readline()
                            merged_file.write(header)
                    # Apeend content from all files
                    for i in range(num_cores):
                        part_file = "{}/Seq_{}_temp_3utr_part{}_pita_results_targets.tab".format(pita_out_dir, seq_num, i + 1)
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
                print("Targetscan is processing {}".format(name_fasta))
            else:
                # Handle other tools
                print("Tool {} is processing {}".format(tool, name_fasta))
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
                        help="List of tools to run. Choose from: {}".format(', '.join(ALLOWED_TOOLS)))
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
        print("Warning: Requested {} cores, but only {} are available.".format(num_cores, available_cores))
        num_cores = available_cores

    # Validate tools
    invalid_tools = [tool for tool in tools if tool not in ALLOWED_TOOLS]
    if invalid_tools:
        raise ValueError("Invalid tools selected: {}. Allowed tools are: {}".format(', '.join(invalid_tools), ', '.join(ALLOWED_TOOLS)))

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    temp_folder = output_folder + "/temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # Parse FASTA files
    sequences = parse_fasta(mirna)
    save_to_json(sequences, tools, num_cores, "{}/mirna_prediction_parameters.json".format(output_folder))
    
    # Determine the 3' UTR file based on the genome
    utr_file = HUMAN_HG19_3UTR if genome == "hg19" else HUMAN_HG38_3UTR

    # Create output directory
    for tool in tools:
        if tool == "miRanda":
            miranda_out_dir = os.path.join(output_folder, "miRanda")
            if not os.path.exists(miranda_out_dir):
                os.makedirs(miranda_out_dir)
        elif tool == "RNAhybrid":
            rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
            if not os.path.exists(rnahybrid_out_dir):
                os.makedirs(rnahybrid_out_dir)
        elif tool == "miRmap":
            mirmap_out_dir = os.path.join(output_folder, "miRmap")
            if not os.path.exists(mirmap_out_dir):
                os.makedirs(mirmap_out_dir)
        elif tool == "DMISO":
            dmiso_out_dir = os.path.join(output_folder, "DMISO")
            if not os.path.exists(dmiso_out_dir):
                os.makedirs(dmiso_out_dir)
        elif tool == "PITA":
            pita_out_dir = os.path.join(output_folder, "PITA")
            if not os.path.exists(pita_out_dir):
                os.makedirs(pita_out_dir)
        elif tool == "Targetscan":
            targetscan_out_dir = os.path.join(output_folder, "Targetscan")
            if not os.path.exists(targetscan_out_dir):
                os.makedirs(targetscan_out_dir)
        else:
            other_out_dir = os.path.join(output_folder, "Other")
            if not os.path.exists(other_out_dir):
                os.makedirs(other_out_dir)
              
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
