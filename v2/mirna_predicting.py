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
import time


# Predefined list of tools
ALLOWED_TOOLS = ["miRanda", "miRmap", "Targetscan", "RNAhybrid", "PITA", "DMISO"]
# Tools that cannot run against an arbitrary lncRNA target pool: TargetScan
# ignores the target FASTA and reads its own precomputed 3' UTR + conservation
# datasets (so it would silently return gene results). PITA is compatible -- it
# scores whatever FASTA is passed to its -utr flag, so the lncRNA reference is
# fed there directly. Restrict TargetScan to gene target-type.
LNCRNA_INCOMPATIBLE_TOOLS = ["Targetscan"]
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
HUMAN_HG19_3UTR = '/opt/reference_files/hsa_HG19_only_3UTRs.fasta'
HUMAN_HG38_3UTR = '/opt/reference_files/hsa_HG38_only_3UTRs.fasta'

ROUNDWORM_3UTR = '/opt/reference_files/cel_WBcel235_3UTRs.fasta'
DOG_3UTR = '/opt/reference_files/cfa_CanFam3.1_3UTRs.fasta'
FRUITFLY_3UTR = '/opt/reference_files/dme_Release6_3UTRs.fasta'
ZEBRAFISH_3UTR = '/opt/reference_files/dre_GRCz11_3UTRs.fasta'
GRAYSHORTTAILEDOPOSSUM_3UTR = '/opt/reference_files/mdo_MonDom5_3UTRs.fasta'
RHESUSMACAQUE_3UTR = '/opt/reference_files/mml_Mmul_8.0.1_3UTRs.fasta'
HOUSEMOUSE_3UTR = '/opt/reference_files/mmu_GRCm38_3UTRs.fasta'
CHIMPANZEE_3UTR = '/opt/reference_files/ptr_Pan_tro3.0_3UTRs.fasta'
NORWAYRAT_3UTR = '/opt/reference_files/rno_RGSC6_rn6_3UTRs.fasta'

# lncRNA PATH (Ensembl ncrna dumps filtered to lncRNA biotypes;
# built by scripts/build_lncrna_references.sh)
HUMAN_HG19_LNCRNA = '/opt/reference_files/hsa_HG19_lncRNAs.fasta'
HUMAN_HG38_LNCRNA = '/opt/reference_files/hsa_HG38_lncRNAs.fasta'

ROUNDWORM_LNCRNA = '/opt/reference_files/cel_WBcel235_lncRNAs.fasta'
DOG_LNCRNA = '/opt/reference_files/cfa_CanFam3.1_lncRNAs.fasta'
FRUITFLY_LNCRNA = '/opt/reference_files/dme_Release6_lncRNAs.fasta'
ZEBRAFISH_LNCRNA = '/opt/reference_files/dre_GRCz11_lncRNAs.fasta'
GRAYSHORTTAILEDOPOSSUM_LNCRNA = '/opt/reference_files/mdo_MonDom5_lncRNAs.fasta'
RHESUSMACAQUE_LNCRNA = '/opt/reference_files/mml_Mmul_8.0.1_lncRNAs.fasta'
HOUSEMOUSE_LNCRNA = '/opt/reference_files/mmu_GRCm38_lncRNAs.fasta'
CHIMPANZEE_LNCRNA = '/opt/reference_files/ptr_Pan_tro3.0_lncRNAs.fasta'
NORWAYRAT_LNCRNA = '/opt/reference_files/rno_RGSC6_rn6_lncRNAs.fasta'

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

_HDR_ACCESSION_RE = re.compile(r'([A-Z]{2,3}_\d+)')


def filter_utr_fasta(utr_file, targets, output_path):
    """Write a filtered copy of utr_file keeping only records whose RefSeq accession
    (version suffix stripped) or 3rd '_'-delimited header token is in targets.

    UTR FASTA headers look like '>ce11_ncbiRefSeq_NM_059873.7 range=...'."""
    target_set = set()
    for t in targets:
        t = t.strip()
        if not t:
            continue
        target_set.add(t)
        m = re.match(r'^([A-Z]{2,3}_\d+)\.\d+$', t)   # strip trailing .version
        if m:
            target_set.add(m.group(1))
    kept = 0
    write = False
    with open(utr_file, 'r') as fin, open(output_path, 'w') as fout:
        for line in fin:
            if line.startswith(">"):
                first = line.split(" ")[0]            # ">ce11_ncbiRefSeq_NM_059873.7"
                ids = set()
                m = _HDR_ACCESSION_RE.search(first)
                if m:
                    ids.add(m.group(1))               # "NM_059873" (no version)
                parts = first.lstrip(">").split("_")
                if len(parts) > 2:
                    ids.add(parts[2])                 # legacy gene-symbol convention
                write = bool(ids & target_set)
                if write:
                    kept += 1
            if write:
                fout.write(line)
    print("Target filter: kept {} UTR records (target set size: {})".format(kept, len(target_set)))
    return output_path


def sanitize_output_path(output_file):
    """Sanitize output filename to avoid illegal characters."""
    # output_dir, filename = os.path.split(output_file)
    # safe_filename = re.sub(r"[^A-Za-z0-9._-]+", "_", filename)
    # if safe_filename in {"", ".", ".."}:
    #     safe_filename = "miranda_output.txt"
    return output_file


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

# RNAhybrid ships only three empirically-calibrated length/energy distributions
# for its -s p-value model: 3utr_human, 3utr_worm and 3utr_fly. There is no
# per-species set beyond these, so map each supported genome to the closest one.
_RNAHYBRID_DATASET_MAP = {
    "dme": "3utr_fly",    # Drosophila melanogaster
    "cel": "3utr_worm",   # Caenorhabditis elegans
}


def rnahybrid_dataset_for_genome(genome):
    """Return the RNAhybrid -s distribution set for a genome/species code.

    Mammals, fish, etc. have no dedicated calibration, so they fall back to
    3utr_human (the closest available proxy)."""
    return _RNAHYBRID_DATASET_MAP.get(genome, "3utr_human")


def run_rnahybrid(mirna_file, utr_file, output_file, mirna_length, utr_length, species_set="3utr_human"):
    """Run RNAhybrid on a given miRNA and UTR file."""
    cmd = [
        RNAHYBRID,
        "-c",
        "-e", "-20",
        "-s", species_set,
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


def run_targetscan_with_params(params):
    return run_targetscan(*params)

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
                # UTR FASTA headers look like ">hg38_ncbiRefSeqCurated_NM_000051.4 range=..."
                # Splitting on "_" and taking [2] only yields "NM" because the
                # RefSeq accession itself contains an underscore. Extract the
                # full RefSeq/Ensembl accession instead so parse_result.py can
                # match it via _TRANSCRIPT_RE.
                first_tok = line.split(" ")[0]
                m = re.search(r'(NM_\d+(?:\.\d+)?|ENST\d+(?:\.\d+)?)', first_tok)
                utr_header = m.group(1) if m else first_tok.lstrip(">")
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
    
def _write_progress(output_folder, tools, tool_statuses):
    progress_path = os.path.join(output_folder, "progress.json")
    # Prediction runs in two phases that share this file: python3.6 for most
    # tools, then python2.7 for miRmap. Merge with any existing on-disk status
    # so the second phase does not clobber the first -- otherwise a finished
    # job's progress.json would list only miRmap and the UI would show the
    # other tools as "pending" after a refresh.
    merged = {}
    try:
        with open(progress_path, 'r') as f:
            prev = json.load(f)
        prev_status = prev.get("tools_status")
        if isinstance(prev_status, dict):
            merged.update(prev_status)
    except (IOError, OSError, ValueError):
        pass
    merged.update(tool_statuses)

    completed = sum(1 for t in merged if merged[t].get("status") == "done")
    current = None
    for t in tools:
        if tool_statuses[t]["status"] == "running":
            current = t
            break
    data = {
        "total_tools": len(merged),
        "completed_tools": completed,
        "current_tool": current,
        "tools_status": merged,
        "updated_at": int(time.time())
    }
    with open(progress_path, 'w') as f:
        json.dump(data, f)


def process_tools(sequences, tools, utr_file, output_folder, temp_folder, rnahybrid_set="3utr_human"):
    tool_statuses = {}
    for t in tools:
        tool_statuses[t] = {"status": "pending", "started_at": None, "finished_at": None}
    _write_progress(output_folder, tools, tool_statuses)
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
                tool_statuses["miRanda"]["status"] = "running"
                tool_statuses["miRanda"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                run_miranda(temp_fasta, utr_file, output_file)
                tool_statuses["miRanda"]["status"] = "done"
                tool_statuses["miRanda"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "RNAhybrid":
                # Output directory for miRanda
                rnahybrid_out_dir = os.path.join(output_folder, "RNAhybrid")
                # Define the output file path
                output_file = "{}/{}_RNAhybrid_results.txt".format(rnahybrid_out_dir, seq['header'])
                # Prepare parameters
                max_utr_length = get_longest_utr_length(utr_file)
                # Run the tool
                print("RNAhybrid is processing {}".format(name_fasta))
                tool_statuses["RNAhybrid"]["status"] = "running"
                tool_statuses["RNAhybrid"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                run_rnahybrid(temp_fasta, utr_file, output_file, int(seq['length']), max_utr_length, rnahybrid_set)
                tool_statuses["RNAhybrid"]["status"] = "done"
                tool_statuses["RNAhybrid"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "miRmap":
                # Create output directory for miRmap
                mirmap_out_dir = os.path.join(output_folder, "miRmap")
                # Define the output file path
                output_file = "{}/{}_miRmap_results.txt".format(mirmap_out_dir, seq['header'])
                # Run the tool
                print("miRmap is processing {}".format(name_fasta))
                tool_statuses["miRmap"]["status"] = "running"
                tool_statuses["miRmap"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                run_mirmap(seq['sequence'], seq['header'], utr_file, output_file)
                tool_statuses["miRmap"]["status"] = "done"
                tool_statuses["miRmap"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "DMISO":         
                # Output directory for DMISO
                dmiso_out_dir = os.path.join(output_folder, "DMISO")
                # Define the output file path
                temp_output_file = "{}/{}_DMISO.txt".format(dmiso_out_dir, seq['header'])
                output_file = "{}/{}_DMISO_results.txt".format(dmiso_out_dir, seq['header'])
                # Run the tool
                print("DMISO is processing {}".format(name_fasta))
                tool_statuses["DMISO"]["status"] = "running"
                tool_statuses["DMISO"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                run_dmiso(temp_fasta, utr_file, temp_output_file)
                # Parse DMISO results
                parse_dmiso_results(temp_output_file, output_file)
                tool_statuses["DMISO"]["status"] = "done"
                tool_statuses["DMISO"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "PITA":
                # Output directory for PITA
                pita_out_dir = os.path.join(output_folder, "PITA")
                # Define the output file path
                output_file_prefix = "{}/{}".format(pita_out_dir, seq['header'])
                # Run the tool
                print("PITA is processing {}".format(name_fasta))
                tool_statuses["PITA"]["status"] = "running"
                tool_statuses["PITA"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                run_pita(temp_fasta, utr_file, output_file_prefix)
                # remove temp file
                if os.path.exists("tmp_seqfile1"):
                    os.remove("tmp_seqfile1")
                if os.path.exists("tmp_seqfile2"):
                    os.remove("tmp_seqfile2")
                if os.path.exists(output_file_prefix + "_pita_results.gxp"):
                    os.remove(output_file_prefix + "_pita_results.gxp")
                tool_statuses["PITA"]["status"] = "done"
                tool_statuses["PITA"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "Targetscan":
                                # Output directory for PITA
                targetscan_out_dir = os.path.join(output_folder, "Targetscan")
                # Define the output file path
                output_file1 = "{}/{}_Targetscan_results1.txt".format(targetscan_out_dir, seq['header'])
                output_file2 = "{}/{}_Targetscan_results2.txt".format(targetscan_out_dir, seq['header'])
                # Run the tool
                print("Targetscan is processing {}".format(name_fasta))
                tool_statuses["Targetscan"]["status"] = "running"
                tool_statuses["Targetscan"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
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
                tool_statuses["Targetscan"]["status"] = "done"
                tool_statuses["Targetscan"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            else:
                # Handle other tools
                print("Tool {} is processing {}".format(tool, name_fasta))
        seq_num += 1

def process_tools_in_parallel(sequences, tools, num_cores, output_folder, temp_folder, rnahybrid_set="3utr_human"):
    # Get all UTR subfiles
    utr_subfiles = [os.path.join(temp_folder+"/utr", f) for f in os.listdir(temp_folder+"/utr") if f.startswith("temp_3utr_part")]

    # Initialize progress tracking
    tool_statuses = {}
    for t in tools:
        tool_statuses[t] = {"status": "pending", "started_at": None, "finished_at": None}
    _write_progress(output_folder, tools, tool_statuses)

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
                tool_statuses["miRanda"]["status"] = "running"
                tool_statuses["miRanda"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                pool.map(run_miranda_with_params, args)

                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(miranda_out_dir, "Seq_{}_miRanda_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                tool_statuses["miRanda"]["status"] = "done"
                tool_statuses["miRanda"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)

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
                    args.append((temp_fasta, utr_file, temp_output_file, seq['length'], max_utr_length, rnahybrid_set))
                
                # Run in parallel
                print("RNAhybrid is processing {}".format(name_fasta))
                tool_statuses["RNAhybrid"]["status"] = "running"
                tool_statuses["RNAhybrid"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                pool.map(run_rnahybrid_with_params, args)

                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(rnahybrid_out_dir, "Seq_{}_RNAhybrid_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                tool_statuses["RNAhybrid"]["status"] = "done"
                tool_statuses["RNAhybrid"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)

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
                tool_statuses["miRmap"]["status"] = "running"
                tool_statuses["miRmap"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
                pool.map(run_mirmap_with_params, args)

                # Merge results
                with open(output_file, 'w') as merged_file:
                    for utr_file in utr_subfiles:
                        part_file = os.path.join(mirmap_out_dir, "Seq_{}_miRmap_results_{}.out".format(seq_num, os.path.basename(utr_file).replace('.fasta', '')))
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                merged_file.write(pf.read())
                            os.remove(part_file)  # Remove temporary part file
                tool_statuses["miRmap"]["status"] = "done"
                tool_statuses["miRmap"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)

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
                tool_statuses["DMISO"]["status"] = "running"
                tool_statuses["DMISO"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
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
                tool_statuses["DMISO"]["status"] = "done"
                tool_statuses["DMISO"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)

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
                tool_statuses["PITA"]["status"] = "running"
                tool_statuses["PITA"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
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
                tool_statuses["PITA"]["status"] = "done"
                tool_statuses["PITA"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
            elif tool == "Targetscan":
                # TargetScan uses its own pre-split 64-part reference at
                # /opt/TargetScan/Datasets/3utr/ -- we parallelize across parts,
                # not across the user UTR file (which TargetScan ignores).
                targetscan_out_dir = os.path.join(output_folder, "Targetscan")
                output_file1 = "{}/{}_Targetscan_results1.txt".format(targetscan_out_dir, seq['header'])
                output_file2 = "{}/{}_Targetscan_results2.txt".format(targetscan_out_dir, seq['header'])

                print("Targetscan is processing {}".format(name_fasta))
                tool_statuses["Targetscan"]["status"] = "running"
                tool_statuses["Targetscan"]["started_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)

                targetscan_prep(seq['sequence'], seq['header'], targetscan_out_dir)
                targetscan_input = "{}/{}_targetscan.txt".format(targetscan_out_dir, seq['header'])
                ts_utr_dir = "/opt/TargetScan/Datasets/3utr"
                ts_bln_dir = "/opt/TargetScan/Datasets/bln_bins"

                ts_args = []
                for i in range(64):
                    utr_part = os.path.join(ts_utr_dir, 'targetscan_utr_part_{}.txt'.format(i))
                    bln_part = os.path.join(ts_bln_dir, 'targetscan_median_bls_bins_part_{}.txt'.format(i))
                    out1 = "{}/{}_part_{}_out1.txt".format(targetscan_out_dir, seq['header'], i)
                    out2 = "{}/{}_part_{}_out2.txt".format(targetscan_out_dir, seq['header'], i)
                    ts_args.append((targetscan_input, utr_part, out1, bln_part, out2))

                pool.map(run_targetscan_with_params, ts_args)

                # Merge per-part outputs (header from part 0, body from all parts).
                with open(output_file1, 'w') as merged:
                    first_file = "{}/{}_part_0_out1.txt".format(targetscan_out_dir, seq['header'])
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            merged.write(first.readline())
                    for i in range(64):
                        part_file = "{}/{}_part_{}_out1.txt".format(targetscan_out_dir, seq['header'], i)
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                next(pf)
                                merged.write(pf.read())
                            os.remove(part_file)

                with open(output_file2, 'w') as merged:
                    first_file = "{}/{}_part_0_out2.txt".format(targetscan_out_dir, seq['header'])
                    if os.path.exists(first_file):
                        with open(first_file, 'r') as first:
                            merged.write(first.readline())
                    for i in range(64):
                        part_file = "{}/{}_part_{}_out2.txt".format(targetscan_out_dir, seq['header'], i)
                        if os.path.exists(part_file):
                            with open(part_file, 'r') as pf:
                                next(pf)
                                merged.write(pf.read())
                            os.remove(part_file)

                tool_statuses["Targetscan"]["status"] = "done"
                tool_statuses["Targetscan"]["finished_at"] = int(time.time())
                _write_progress(output_folder, tools, tool_statuses)
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
    parser.add_argument("-g", "--genome", type=str,
                        choices=["hg19", "hg38", "cel", "cfa", "dme", "dre", "mdo", "mml", "mmu", "ptr", "rno"],
                        help="Reference genome/species code", required=True)
    parser.add_argument("-o", "--output", type=str, required=True, help="output folder name")
    parser.add_argument("-tf", "--target_file", type=str, default=None,
                        help="optional path to a text file with one target gene symbol per line; only 3' UTRs matching these genes will be scanned (gene target-type only)")
    parser.add_argument("-tt", "--target-type", dest="target_type",
                        choices=["gene", "lncrna"], default="gene",
                        help="target pool to scan: 'gene' = 3' UTRs (default, "
                             "legacy behavior), 'lncrna' = lncRNA transcripts")

    args = parser.parse_args()
    # Parse arguments
    num_cores = args.cores
    mirna = args.input
    tools = args.tools
    output_folder = args.output
    genome = args.genome
    target_type = args.target_type

    # Validate number of cores
    available_cores = multiprocessing.cpu_count()
    if num_cores > available_cores:
        print("Warning: Requested {} cores, but only {} are available.".format(num_cores, available_cores))
        num_cores = available_cores

    # Validate tools
    invalid_tools = [tool for tool in tools if tool not in ALLOWED_TOOLS]
    if invalid_tools:
        raise ValueError("Invalid tools selected: {}. Allowed tools are: {}".format(', '.join(invalid_tools), ', '.join(ALLOWED_TOOLS)))

    # Reject tools that can't honor a lncRNA target pool (see LNCRNA_INCOMPATIBLE_TOOLS).
    if target_type == "lncrna":
        bad = [tool for tool in tools if tool in LNCRNA_INCOMPATIBLE_TOOLS]
        if bad:
            raise ValueError(
                "Tools {} are not supported for --target-type lncrna; "
                "they only work on 3' UTR (gene) targets. Drop them or use "
                "--target-type gene.".format(', '.join(bad))
            )

    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    temp_folder = output_folder + "/temp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    # Parse FASTA files
    sequences = parse_fasta(mirna)
    save_to_json(sequences, tools, num_cores, "{}/mirna_prediction_parameters.json".format(output_folder))
    
    # Reference sequence pool, keyed by (target_type, genome/species code).
    # "gene"   -> 3' UTRs of protein-coding genes (legacy default)
    # "lncrna" -> lncRNA transcripts (Ensembl, biotype-filtered)
    _TARGET_FILE_MAP = {
        "gene": {
            "hg19": HUMAN_HG19_3UTR,
            "hg38": HUMAN_HG38_3UTR,
            "cel":  ROUNDWORM_3UTR,
            "cfa":  DOG_3UTR,
            "dme":  FRUITFLY_3UTR,
            "dre":  ZEBRAFISH_3UTR,
            "mdo":  GRAYSHORTTAILEDOPOSSUM_3UTR,
            "mml":  RHESUSMACAQUE_3UTR,
            "mmu":  HOUSEMOUSE_3UTR,
            "ptr":  CHIMPANZEE_3UTR,
            "rno":  NORWAYRAT_3UTR,
        },
        "lncrna": {
            "hg19": HUMAN_HG19_LNCRNA,
            "hg38": HUMAN_HG38_LNCRNA,
            "cel":  ROUNDWORM_LNCRNA,
            "cfa":  DOG_LNCRNA,
            "dme":  FRUITFLY_LNCRNA,
            "dre":  ZEBRAFISH_LNCRNA,
            "mdo":  GRAYSHORTTAILEDOPOSSUM_LNCRNA,
            "mml":  RHESUSMACAQUE_LNCRNA,
            "mmu":  HOUSEMOUSE_LNCRNA,
            "ptr":  CHIMPANZEE_LNCRNA,
            "rno":  NORWAYRAT_LNCRNA,
        },
    }
    utr_file = _TARGET_FILE_MAP[target_type][genome]

    if args.target_file:
        # The target-gene filter keys on RefSeq accessions in 3' UTR headers;
        # lncRNA targets use Ensembl headers and are not gene-symbol filterable.
        if target_type != "gene":
            raise ValueError(
                "--target_file is only supported for --target-type gene, "
                "not '{}'.".format(target_type)
            )
        with open(args.target_file, 'r') as f:
            target_genes = [line.strip() for line in f if line.strip()]
        if target_genes:
            filtered_utr_path = os.path.join(temp_folder, "filtered_3utr.fasta")
            utr_file = filter_utr_fasta(utr_file, target_genes, filtered_utr_path)
            print("Scanning {} target genes from: {}".format(len(target_genes), args.target_file))

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
              
    # Select the RNAhybrid distribution set for this species (see rnahybrid_dataset_for_genome)
    rnahybrid_set = rnahybrid_dataset_for_genome(genome)

    # Run Prediction for single or mutiple cores
    if num_cores == 1:
        process_tools(sequences, tools, utr_file, output_folder, temp_folder, rnahybrid_set)
    else:
        process_3utr_fasta(utr_file, num_cores, temp_folder)
        process_tools_in_parallel(sequences, tools, num_cores, output_folder ,temp_folder, rnahybrid_set)

    # Clean up the temp folder
    cleanup_temp_folder(temp_folder)

    print("Processing complete.")

if __name__ == "__main__":
    main()
