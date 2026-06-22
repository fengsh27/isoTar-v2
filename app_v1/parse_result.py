#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import re
import json
import csv
import sys

# Matches Ensembl (ENST00000284637) or RefSeq mRNA (NM_001164664), stripping version suffix
_TRANSCRIPT_RE = re.compile(r'(ENST\d+|NM_\d+)(?:\.\d+)?')
_REFSEQ_RE = re.compile(r'^[A-Z]{2,3}_\d+(\.\d+)?$')

def _extract_transcript_id(text):
    """Return ENST or NM_ transcript ID from text, without version suffix. None if not found."""
    m = _TRANSCRIPT_RE.search(text)
    return m.group(1) if m else None


# lncRNA reference FASTAs come from Ensembl /ncrna/ dumps, so their target IDs
# are species-specific Ensembl transcripts (human ENST, mouse ENSMUST, dog
# ENSCAFT, rat ENSRNOT, ...), FlyBase transcripts (FBtr0347114) or WormBase
# sequence names (Y51H4A.27) -- none of which the human-centric ENST/NM_ regex
# above matches. This anchored pattern strips a trailing Ensembl/FlyBase
# '.<version>' suffix; WormBase isoform suffixes are left intact (the '.27' is
# part of the name, not a version).
_LNCRNA_TRANSCRIPT_RE = re.compile(r'^(ENS[A-Z]{0,6}T\d+|FBtr\d+)(?:\.\d+)?')


def _extract_lncrna_transcript_id(text):
    """Transcript ID for a lncRNA target token. Returns the version-stripped
    Ensembl/FlyBase accession, or the first whitespace-delimited token as-is
    (e.g. WormBase Y51H4A.27). None for empty input.

    Some tools (DMISO) keep the full FASTA description in the ID field, so the
    pattern is anchored at the start of the token rather than searched."""
    if not text:
        return None
    t = text.strip()
    if not t:
        return None
    m = _LNCRNA_TRANSCRIPT_RE.match(t)
    if m:
        return m.group(1)
    return t.split()[0]


def _default_reference_db():
    """Resolve reference_mapping.db path. Tries env var, then app_v1/, then Docker path."""
    env = os.environ.get("ISOTAR_REFERENCE_MAPPING_DB")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    # app_v1/parse_result.py -> app_v1/reference_mapping.db
    local = os.path.join(here, "reference_mapping.db")
    if os.path.exists(local):
        return local
    # v2/parse_result.py -> ../app_v1/reference_mapping.db
    sibling = os.path.join(here, "..", "app_v1", "reference_mapping.db")
    if os.path.exists(sibling):
        return sibling
    return "/app_v1/reference_mapping.db"


def build_enst_to_refseq_map(ref_db_path=None):
    """Return dict mapping unversioned ENST -> set of unversioned RefSeq IDs.

    Two indexed scans (ensembl_mapping, gene_mapping) joined on UPPER(symbol)
    in Python. Avoids `COLLATE NOCASE` on the JOIN which defeats the symbol
    index and turns this into a multi-billion-comparison nested loop."""
    import sqlite3
    if ref_db_path is None:
        ref_db_path = _default_reference_db()
    mp = {}
    if not os.path.exists(ref_db_path):
        return mp
    conn = sqlite3.connect(ref_db_path)
    try:
        c = conn.cursor()
        # symbol(uppercased) -> set of unversioned RefSeq IDs.
        # Restrict to human species — TargetScan output is filtered to 9606
        # in parseTargetScanResults, so cross-species symbol collisions
        # (e.g. mouse "Trp53" sharing a HGNC bucket) would only add noise.
        symbol_to_refseqs = {}
        c.execute(
            "SELECT symbol, raw_id FROM gene_mapping "
            "WHERE species LIKE 'hsa_%' AND symbol IS NOT NULL AND raw_id IS NOT NULL"
        )
        for sym, raw_id in c.fetchall():
            symbol_to_refseqs.setdefault(sym.upper(), set()).add(raw_id.split(".")[0])
        # walk ensembl_mapping, project each ENST onto its symbol's RefSeq set
        c.execute("SELECT ensembl_id, symbol FROM ensembl_mapping WHERE ensembl_id IS NOT NULL AND symbol IS NOT NULL")
        for ensembl_id, sym in c.fetchall():
            refseqs = symbol_to_refseqs.get(sym.upper())
            if not refseqs:
                continue
            enst = ensembl_id.split(".")[0]
            bucket = mp.setdefault(enst, set())
            bucket.update(refseqs)
    finally:
        conn.close()
    return mp


def load_targets_file(targets_file):
    """Read targets.txt (one ID per line). Return set of unversioned RefSeq IDs,
    or None if the file is absent (meaning 'no target filter')."""
    if not targets_file or not os.path.exists(targets_file):
        return None
    targets = set()
    with open(targets_file, "r") as f:
        for line in f:
            t = line.strip()
            if not t:
                continue
            if _REFSEQ_RE.match(t):
                targets.add(t.split(".")[0])
    return targets if targets else None


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

def parseTargetScanResults(output_f_path, result_dict, enst_to_refseq=None, targets=None):
    """Parse TargetScan output, optionally converting ENST IDs to RefSeq.

    enst_to_refseq: dict[ENST -> set[RefSeq]]. If provided, each ENST hit is
        expanded to its RefSeq IDs; ENSTs not in the map are dropped.
    targets: set[RefSeq]. If provided, only RefSeq IDs in this set are kept.
        Has no effect when enst_to_refseq is None (we'd be filtering ENSTs
        against a RefSeq set — would drop everything)."""
    raw_hits = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            handler = csv.reader(f, delimiter='\t')
            next(handler) # Skip header row
            for line in handler:
                if len(line) == 14:
                    if line[2] =="9606" and line[8] != '6mer':
                        # Target - remove version number if present
                        tar = re.sub(r'\.[0-9]+$', '', line[0])
                        if tar and tar not in raw_hits:
                            raw_hits.append(tar)

    if enst_to_refseq is not None:
        results = []
        seen = set()
        for enst in raw_hits:
            for refseq in enst_to_refseq.get(enst, ()):
                if targets is not None and refseq not in targets:
                    continue
                if refseq not in seen:
                    seen.add(refseq)
                    results.append(refseq)
    else:
        results = raw_hits

    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['TargetScan'] = results

    return result_dict

def parsePITAResults(output_f_path, result_dict, id_extractor=_extract_transcript_id):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            handler = csv.reader(f, delimiter='\t')
            for line in handler:
                if len(line) == 13:
                    tar = id_extractor(line[0])
                    if tar:
                        ddG = float(line[12])
                        if ddG <= -10.0:
                            if tar not in results:
                                results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['PITA'] = results

    return result_dict

def parseRnahybridResults(output_f_path, result_dict, id_extractor=_extract_transcript_id):
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            handler = csv.reader(f, delimiter=':')
            for line in handler:
                if len(line) == 11:
                    tar = id_extractor(line[0])
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

def parseMirmapResults(output_f_path, result_dict, id_extractor=_extract_transcript_id):
    results = []
    if os.path.exists(output_f_path):
        # miRmap writes "DeltaG binding (kcal/mol)" with a literal Greek Delta (UTF-8 0xCE 0x94),
        # so reading without an explicit encoding fails on containers without a UTF-8 locale.
        with open(output_f_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                # miRNA - Target
                matchObj = re.match(r'^>[^,]+,.*?\s+(\S+)\s*$', lines[i], re.M|re.I)
                if matchObj:
                    tar = id_extractor(matchObj.group(1))
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

def parseMirandaResults(output_f_path, result_dict, id_extractor=_extract_transcript_id):
    """Extract target accessions from a miRanda -out file.

    miRanda -quiet writes only the summary lines: '>' per individual hit and
    '>>' per (miRNA, target) total, both tab-separated with the UTR accession
    (e.g. 'hg38_ncbiRefSeqCurated_NM_000051.4') in column 1. The verbose
    (non-quiet) format also emits these same summary lines after each
    alignment block, so this parser handles both invocation styles.
    """
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    continue
                parts = line.rstrip('\r\n').split('\t')
                if len(parts) < 2:
                    continue
                tar = id_extractor(parts[1])
                if tar and tar not in results:
                    results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['miRanda'] = results

    return result_dict

# DNA/IUPAC complement table for the DMISO seed-match filter. We complement a
# short (7 nt) string, so a dict lookup is fine and keeps us off str.maketrans
# (absent in Python 2) for 2.7/3.5 compatibility.
_DMISO_COMPLEMENT = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
    'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S', 'W': 'W',
    'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
}


def _mirna_seed_match_pattern(mirna_seq):
    """Reverse-complement of the miRNA seed (positions 2-8) as a DNA 7-mer.

    This is the motif a target site must contain for a canonical seed match.
    Returns '' when the sequence is missing or shorter than 8 nt, which the
    caller treats as 'apply no seed filter'."""
    if not mirna_seq:
        return ''
    s = mirna_seq.upper().replace('U', 'T')
    if len(s) < 8:
        return ''
    seed = s[::-1][-8:-1]  # 7 nt corresponding to miRNA positions 2-8
    return ''.join(_DMISO_COMPLEMENT.get(b, 'N') for b in seed)


def parseDMISOResults(output_f_path, result_dict, mirna_sequence=None, id_extractor=_extract_transcript_id):
    """Parse DMISO output, keeping only targets with a perfect seed match.

    DMISO output is tab-separated: Target ID / Target Sequence / Prediction
    Score (score already filtered > 0.99 upstream by mirna_predicting.py). A
    target is kept only when its sequence contains the reverse-complement of
    the miRNA seed (positions 2-8). When mirna_sequence is missing or too
    short, no seed filter is applied and all listed targets are kept."""
    results = []
    seen = set()  # O(1) dedup: DMISO files can hold ~400k rows
    pattern = _mirna_seed_match_pattern(mirna_sequence)
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            handler = csv.reader(f, delimiter='\t')
            next(handler)  # Skip header row
            for line in handler:
                if len(line) < 2:
                    continue
                if pattern and pattern not in line[1].upper():
                    continue
                tar = id_extractor(line[0])
                if tar and tar not in seen:
                    seen.add(tar)
                    results.append(tar)
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['DMISO'] = results

    return result_dict


def process_sequence(sequence, result_dir, enst_to_refseq=None, targets=None):
    # Process a single sequence and generate prediction results.
    try:
        output_f_path_miRanda = os.path.join(result_dir, "miRanda", "{}_miRanda_results.txt".format(sequence['header']))
        output_f_path_miRmap = os.path.join(result_dir, "miRmap", "{}_miRmap_results.txt".format(sequence['header']))
        output_f_path_RNAhybrid = os.path.join(result_dir, "RNAhybrid", "{}_RNAhybrid_results.txt".format(sequence['header']))
        output_f_path_PITA = os.path.join(result_dir, "PITA", "{}_PITA_results.tab".format(sequence['header']))
        output_f_path_TargetScan = os.path.join(result_dir, "Targetscan", "{}_Targetscan_results1.txt".format(sequence['header']))
        output_f_path_DMISO = os.path.join(result_dir, "DMISO", "{}_DMISO_results.txt".format(sequence['header']))

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
            prediction_results = parseTargetScanResults(
                output_f_path_TargetScan, prediction_results,
                enst_to_refseq=enst_to_refseq, targets=targets,
            )
        if os.path.exists(output_f_path_DMISO):
            prediction_results = parseDMISOResults(
                output_f_path_DMISO, prediction_results,
                mirna_sequence=sequence.get('sequence'),
            )

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
    parser.add_argument('--targets-file', default=None,
                        help='Optional targets.txt (one RefSeq per line). '
                             'Defaults to <parent of result_dir>/targets.txt if present.')
    args = parser.parse_args()
    result_dir = args.result_dir
    output_dir = os.path.join(result_dir, "miRNA_prediction_results")

    # Resolve targets file: explicit arg > sibling of result_dir
    targets_file = args.targets_file
    if targets_file is None:
        candidate = os.path.join(os.path.dirname(os.path.abspath(result_dir)), "targets.txt")
        if os.path.exists(candidate):
            targets_file = candidate
    targets = load_targets_file(targets_file)
    enst_to_refseq = build_enst_to_refseq_map()
    if targets:
        print("Loaded {} target RefSeq IDs from {}".format(len(targets), targets_file))
    print("ENST->RefSeq map: {} ENST entries".format(len(enst_to_refseq)))

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
                prediction_results = process_sequence(
                    sequence, result_dir,
                    enst_to_refseq=enst_to_refseq, targets=targets,
                )
                
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
