#!/usr/bin/env python
# -*- coding: utf-8 -*-
import io
import os
import re
import json
import csv
import sys

# Matches Ensembl (ENST00000284637) or RefSeq mRNA (NM_001164664), stripping version suffix
_TRANSCRIPT_RE = re.compile(r'(ENST\d+|NM_\d+)(?:\.\d+)?')
_REFSEQ_RE = re.compile(r'^[A-Z]{2,3}_\d+(\.\d+)?$')
# miRmap report: one header per (miRNA, transcript) pair, then one block per
# predicted site. The transcript accession is the header's last field. The
# header historically looked like ">mir,WT NM_000051", but the ",WT" suffix is
# absent for some run types, so the comma is not required here.
_MIRMAP_HEADER_RE = re.compile(r'^>.*\s+(\S+)\s*$')
# Emitted exactly once per site by both the miRmap 1.x report layout and v2's
# _build_mirmap2_block, so it is a reliable site delimiter. Group 1 is the site's
# ΔG binding, which _MIRMAP_DG_BINDING_MAX gates on.
_MIRMAP_SITE_RE = re.compile(r'^\s*ΔG binding \(kcal/mol\)\s+([-+]?\d+\.?\d*)\s*$')
# Keep a transcript only if some site binds at least this strongly (kcal/mol).
# Set to None to accept every seed match (isoTar v1's behaviour).
_MIRMAP_DG_BINDING_MAX = -20.0

def _extract_transcript_id(text):
    """Return ENST or NM_ transcript ID from text, without version suffix. None if not found."""
    m = _TRANSCRIPT_RE.search(text)
    return m.group(1) if m else None


def _default_reference_db():
    """Resolve reference_mapping.db path. Tries env var, then app_v1/, then Docker path."""
    env = os.environ.get("ISOTAR_REFERENCE_MAPPING_DB")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    # v2/parse_result.py -> ../app_v1/reference_mapping.db
    sibling = os.path.join(here, "..", "app_v1", "reference_mapping.db")
    if os.path.exists(sibling):
        return sibling
    # Fallback: co-located
    local = os.path.join(here, "reference_mapping.db")
    if os.path.exists(local):
        return local
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
        Has no effect when enst_to_refseq is None."""
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

def parsePITAResults(output_f_path, result_dict):
    results = []
    rows_seen = 0
    dgopen_seen = 0
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            handler = csv.reader(f, delimiter='\t')
            for line in handler:
                if len(line) == 13:
                    tar = _extract_transcript_id(line[0])
                    if tar:
                        rows_seen += 1
                        # col 11 = dGopen (site accessibility); empty when the
                        # RNAddG4 binary fails to emit output (see guard below).
                        if line[11].strip() != '':
                            dgopen_seen += 1
                        ddG = float(line[12])
                        if ddG <= -10.0:
                            if tar not in results:
                                results.append(tar)
    # Accessibility guard: PITA's ddG = dGduplex - dGopen. If the RNAddG4
    # binary fails (it segfaults in some images), dGopen comes back empty and
    # ddG silently collapses to dGduplex, so PITA over-predicts ~20-40x. Fail
    # loudly here rather than emit invalid targets. run_pita discards PITA's
    # stderr, so a die() inside RNAddG_compute.pl would be swallowed -- this
    # wrapper-side check is the reliable place to catch it.
    if rows_seen > 0 and dgopen_seen == 0:
        raise RuntimeError(
            "PITA: dGopen is empty for all {} site rows in {} -- the RNAddG4 "
            "accessibility binary produced no output, so ddG collapsed to "
            "dGduplex and PITA predictions are invalid (massive over-prediction). "
            "Replace the broken RNAddG4 binary at /opt/PITA64bit/Bin/ViennaRNA/"
            "ViennaRNA-1.6/Progs/RNAddG4 and re-run.".format(rows_seen, output_f_path)
        )
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
                # A compact-format record is 11 colon-separated fields, but
                # field 2 is the miRNA name and isomiR names carry colons
                # (e.g. "hsa-let-7a-5p,10:A|C,modified"), so the count grows
                # with them. Both ends are stable -- the target is first and
                # the four alignment rows are last -- so anchor from the right.
                # On an 11-field line line[-3]/line[-2] are line[8]/line[9].
                if len(line) >= 11:
                    tar = _extract_transcript_id(line[0])
                    if tar:
                        # Get the seed region 2-7
                        target_seq = line[-3][-8:-1]
                        mirna_seq = line[-2][-8:-1]
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
    """Collect transcripts with at least one site binding at or below
    _MIRMAP_DG_BINDING_MAX (kcal/mol).

    isoTar v1 applied NO energy threshold -- it kept every transcript with a
    seed match -- so this gate is a deliberate divergence from v1, kept because
    the weak tail is not worth reporting. It is only meaningful now that the
    prediction side asks for 7/8mer seeds (mirna_predicting.run_mirmap). Under
    miRmap 2's [6,7] default the seed pool was dominated by 6mers, which bind
    far more weakly: ΔG binding ran min/median/max -24.94/-7.96/-3.62 and -20
    sat past the end of that distribution, discarding ~99.9% of the output.
    With [7,8] the same measurement gives -24.02/-15.81/-11.48, where -20 is a
    stringent but real tail cutoff (~10% of hit transcripts on a 2428-UTR
    sample; the exact share is miRNA-dependent, since a GC-rich miRNA clears
    -20 more often). Do NOT restore this gate without the seed-length fix.

    A transcript qualifies on its BEST site, not its first. The pre-2026-08
    implementation read a fixed offset (+8) from the header, so it only ever
    tested a transcript's first site -- and miRmap 2 returns sites sorted by
    seed end descending (nearest the UTR 3' end first), which has nothing to do
    with binding energy. Scanning line by line also drops the assumption that
    every block has a fixed number of feature lines: the 1.x and 2.x layouts
    differ, so those offsets were only accidentally correct.
    """
    results = []
    seen = set()
    if os.path.exists(output_f_path):
        # miRmap writes "ΔG binding (kcal/mol)" with a literal Greek Delta (UTF-8 0xCE 0x94),
        # so reading without an explicit encoding fails on containers without a UTF-8 locale.
        # io.open is used (not builtin open) so the encoding= kwarg works under Py2.7 too.
        with io.open(output_f_path, 'r', encoding='utf-8') as f:
            tar = None
            for line in f:
                header = _MIRMAP_HEADER_RE.match(line)
                if header:
                    tar = _extract_transcript_id(header.group(1))
                    continue
                if tar is None:
                    continue
                site = _MIRMAP_SITE_RE.match(line)
                if not site:
                    continue
                if (_MIRMAP_DG_BINDING_MAX is not None
                        and float(site.group(1)) > _MIRMAP_DG_BINDING_MAX):
                    # Too weak -- keep scanning this transcript's other sites.
                    continue
                if tar not in seen:
                    seen.add(tar)
                    results.append(tar)
                # Recorded -- skip this transcript's remaining sites.
                tar = None
    if 'prediction' not in result_dict:
        result_dict["prediction"] = {}
    result_dict["prediction"]['miRmap'] = results

    return result_dict

def parseMirandaResults(output_f_path, result_dict):
    """Extract target accessions from a miRanda -out file, keeping a target
    only when the miRNA seed is fully Watson-Crick paired -- the canonical
    seed-match filter isoTar v1 applied. (v1's own comment says "seed region
    2-7", but its slice actually spans miRNA positions 2-8, i.e. 7mer-m8. The
    behaviour here is kept identical to v1; only the label is corrected.)

    Reading every '>' summary line (all hits above the -sc/-en thresholds)
    over-predicts ~7x relative to v1, because miRanda reports many alignments
    with an imperfect seed. Each hit emits an alignment block: a 'Query:' line,
    the pairing string on the next line, a 'Ref:' line, then (5 lines further)
    the '>' summary line whose 2nd tab field is the target id. miranda v3.3a
    emits these blocks even under -quiet. A truly summary-only output has no
    alignment block and yields nothing here -- acceptable, since the seed
    cannot be verified without it.
    """
    results = []
    if os.path.exists(output_f_path):
        with open(output_f_path, 'r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if not re.match(r"^\s+Query:\s+3'\s+\S+\s+5'$", lines[i], re.I):
                continue
            # The pairing string is the next line. miRanda prints the miRNA
            # 3'->5', so position 1 is rightmost; the line ends with a trailing
            # space (pos 1) then '\n'. Slice [-9:-2] therefore spans miRNA
            # positions 2-8. It must be all '|' (Watson-Crick): a ' ' (mismatch)
            # or ':' (G:U wobble) anywhere in that window rejects the hit.
            if i + 1 >= len(lines) or lines[i + 1][-9:-2].replace('|', '') != '':
                continue
            # Query(+0) -> pairing(+1) -> Ref(+2) -> +5 lands on the '>' line.
            j = i + 7
            if j >= len(lines) or not lines[j].startswith('>'):
                continue
            parts = lines[j].rstrip('\r\n').split('\t')
            if len(parts) < 2:
                continue
            tar = _extract_transcript_id(parts[1])
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


def parseDMISOResults(output_f_path, result_dict, mirna_sequence=None):
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
                tar = _extract_transcript_id(line[0])
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

def _build_label_map(ref_db_path):
    """Return dict mapping raw_id -> symbol from reference_mapping.db, or {} if unavailable."""
    import sqlite3
    label_map = {}
    candidates = [ref_db_path]
    # Fallback: Docker path
    if not os.path.exists(ref_db_path):
        candidates.append("/app_v1/reference_mapping.db")
    for path in candidates:
        if os.path.exists(path):
            conn = sqlite3.connect(path)
            try:
                c = conn.cursor()
                c.execute("SELECT raw_id, symbol FROM gene_mapping WHERE symbol IS NOT NULL")
                for raw_id, symbol in c.fetchall():
                    label_map[raw_id] = symbol
            finally:
                conn.close()
            break
    return label_map


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Process miRNA sequences and generate prediction results')
    parser.add_argument('result_dir', help='Path to miRNA prediction results')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    parser.add_argument('--gene-label', action='store_true',
                        help='Replace gene IDs with gene labels (symbols) from reference_mapping.db')
    parser.add_argument('--targets-file', default=None,
                        help='Optional targets.txt (one RefSeq per line). '
                             'Defaults to <parent of result_dir>/targets.txt if present.')
    args = parser.parse_args()
    result_dir = args.result_dir
    output_dir = os.path.join(result_dir, "miRNA_prediction_results")

    # Build label map if --gene-label requested
    label_map = {}
    if args.gene_label:
        ref_db = _default_reference_db()
        label_map = _build_label_map(ref_db)
        if not label_map:
            print("Warning: --gene-label requested but reference_mapping.db not found or empty")

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

                # Replace gene IDs with gene labels if requested
                if label_map and "prediction" in prediction_results:
                    for tool in prediction_results["prediction"]:
                        prediction_results["prediction"][tool] = [
                            label_map.get(gid, gid)
                            for gid in prediction_results["prediction"][tool]
                        ]

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
