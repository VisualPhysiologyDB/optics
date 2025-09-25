import subprocess
import os
import tempfile
import pandas as pd
import json
import argparse
import io
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
from joblib import Parallel, delayed

# --- Configuration & Constants ---
# Predefined reference sequences for common use cases.
REFERENCE_SEQUENCES = {
    "bovine": '>Bovine_Rhodopsin\nMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA',
    "squid": '>Squid_Rhodopsin\nMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA',
    "microbe": '>Acetabularia_Rhodopsin_1\nMSNPNPFQTTLGTDAQWVVFAVMALAAIVFSIAVQFRPLPLRLTYYVNIAICTIAATAYYAMAVNGGDNKPTAGTGADERQVIYARYIDWVFTTPLLLLDLVLLTNMPATMIAWIMGADIAMIAFGIIGAFTVGSYKWFYFVVGCIMLAVLAWGMINPIFKEELQKHKEYTGAYTTLLIYLIVLWVIYPIVWGLGAGGHIIGVDVEIIAMGILDLLAKPLYAIGVLITVEVVYGKLGKEEAQPLTA'
}

# --- BLAST Caching Functions ---

def get_cache_path(database_path, wrk_dir):
    """Generates a cache filename based on the database name."""
    db_name = Path(database_path).stem        
    return Path(f"{wrk_dir}/data/cached_blastp_analysis/blast_cache_{db_name}.json")

def load_blast_cache(cache_path):
    """Loads the BLAST results cache from a JSON file."""
    if cache_path.exists():
        with open(cache_path, 'r') as f:
            return json.load(f)
    return {}

def save_blast_cache(cache, cache_path):
    """Saves the BLAST results cache to a JSON file."""
    with open(cache_path, 'w') as f:
        json.dump(cache, f, indent=4)

# --- Core Bioinformatic Functions ---

def run_blastp_bulk(query_file, database):
    """
    Performs a single BLASTp search for all sequences in a query file.
    """
    print("INFO: Running blastp for all query sequences...")
    try:
        # We request specific columns for easier parsing later.
        outfmt_str = "6 qseqid sseqid pident length mismatch gapopen evalue bitscore"
        blast_cmd = ["blastp", "-query", query_file, "-db", database, "-outfmt", outfmt_str]
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        return result.stdout
    except FileNotFoundError:
        print("ERROR: 'blastp' command not found. Make sure BLAST+ is installed and in your system's PATH.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BLASTp failed with the following error:\n{e.stderr}")
        return None

def filter_best_blast_hits(blast_output_string):
    """
    Parses BLAST output and filters for the single best hit per query sequence.
    The best hit is determined by the lowest e-value, with bitscore as a tie-breaker.
    """
    if not blast_output_string.strip():
        print("WARNING: BLASTp returned no hits.")
        return pd.DataFrame()

    col_names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore"]
    try:
        # Use an in-memory text buffer to read the string output into pandas
        blast_results_df = pd.read_csv(io.StringIO(blast_output_string), sep="\t", header=None, names=col_names)
    except Exception as e:
        print(f"ERROR: Failed to parse BLAST output. {e}")
        return pd.DataFrame()

    # Sort by e-value (ascending) and bitscore (descending) to rank hits
    blast_results_df.sort_values(by=['evalue', 'bitscore'], ascending=[True, False], inplace=True)
    
    # Keep only the first occurrence of each query ID, which is now the best hit
    best_hits_df = blast_results_df.drop_duplicates(subset='qseqid', keep='first').reset_index(drop=True)
    
    #print(f"INFO: Filtered BLAST results down to {len(best_hits_df)} best hits.")
    return best_hits_df


def align_sequences(query_record, target_record, reference_record, temp_dir, wrk_dir):
    """Aligns query, target, and reference sequences using MAFFT, with fallbacks."""
    with tempfile.NamedTemporaryFile(mode="w", dir=temp_dir, suffix=".fasta", delete=False) as temp_input_file:
        temp_input_path = temp_input_file.name
        temp_input_file.write(f"{reference_record}\n")
        temp_input_file.write(f"{query_record}\n")
        temp_input_file.write(f"{target_record}\n")

    with tempfile.NamedTemporaryFile(mode="w", dir=temp_dir, suffix=".fasta", delete=False) as temp_output_file:
        temp_output_path = temp_output_file.name
    
    mafft_executables = [
        'mafft',
        str(Path(wrk_dir) / 'optics_scripts/mafft/mafft-win/mafft.bat'),
        str(Path(wrk_dir) / 'optics_scripts/mafft/mafft-mac/mafft.bat')
    ]

    alignment_successful = False
    last_error = ""
    for exe in mafft_executables:
        try:
            mafft_cmd = [exe, '--auto', temp_input_path]
            with open(temp_output_path, 'w') as f_out:
                subprocess.run(mafft_cmd, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)
            alignment_successful = True
            break
        except FileNotFoundError:
            last_error = f"Executable not found at '{exe}'."
        except subprocess.CalledProcessError as e:
            last_error = e.stderr
            continue
    
    os.remove(temp_input_path)

    if not alignment_successful:
        print(f"ERROR: MAFFT alignment failed for all options. Last error:\n{last_error}")
        os.remove(temp_output_path)
        return None

    return temp_output_path

def analyze_differences(aligned_query_seq, aligned_target_seq, aligned_ref_seq):
    """
    Analyzes differences between two aligned sequences and maps them to a
    reference sequence's non-gap positions.
    """
    differences = []
    ref_pos_counter = 0
    for i, (ref_char, query_char, target_char) in enumerate(zip(aligned_ref_seq, aligned_query_seq, aligned_target_seq)):
        is_ref_gap = (ref_char == '-')
        
        if not is_ref_gap:
            ref_pos_counter += 1

        if query_char != target_char:
            position_label = ref_pos_counter if not is_ref_gap else 'N/A'
            if position_label != 'N/A':
                diff_string = f"{query_char}{position_label}{target_char}"
                differences.append(diff_string)
            
    return "; ".join(differences) if differences else "No Differences"

# --- Data Extraction and Processing ---

def get_reference_record(ref_id, ref_file_path):
    """Retrieves the appropriate reference sequence record."""
    if ref_id in REFERENCE_SEQUENCES:
        return REFERENCE_SEQUENCES[ref_id]
    elif ref_file_path and Path(ref_file_path).exists():
        try:
            record = next(SeqIO.parse(ref_file_path, "fasta"))
            return f">{record.description}\n{str(record.seq)}"
        except StopIteration:
            print(f"ERROR: Custom reference file '{ref_file_path}' is empty.")
            return None
    else:
        print(f"ERROR: Could not find reference '{ref_id}' or file '{ref_file_path}'.")
        return None

def extract_record_from_fasta(fasta_file, seq_id):
    """Efficiently retrieves a single sequence record from a FASTA file."""
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Match the ID up to the first space to handle diverse FASTA headers
            if record.id.split()[0] == seq_id.split()[0]:
                return f">{record.description}\n{str(record.seq)}"
    except FileNotFoundError:
        print(f"ERROR: FASTA database file not found at {fasta_file}")
    return None

def align_and_analyze_one(query_record, closest_match_id, blast_metrics, ref_seq_id, ref_file, opsin_fasta, metadata_df, temp_dir, wrk_dir):
    """
    Worker function to process a single query.
    This function performs alignment, analysis, and metadata retrieval.
    It's designed to be called in parallel.
    """
    query_id = query_record.id
    
    # 1. Extract sequences for alignment
    query_fasta_record = f">{query_record.description}\n{str(query_record.seq)}"
    closest_match_record = extract_record_from_fasta(opsin_fasta, closest_match_id)
    reference_record = get_reference_record(ref_seq_id, ref_file)

    if not all([closest_match_record, reference_record]):
        return { "query_id": query_id, "status": "Failed to extract sequences for alignment" }

    # 2. Alignment
    alignment_file = align_sequences(query_fasta_record, closest_match_record, reference_record, temp_dir, wrk_dir)
    if not alignment_file:
        return { "query_id": query_id, "status": "MAFFT Alignment failed" }

    # 3. Parse alignment and analyze differences
    try:
        aligned_records = list(SeqIO.parse(alignment_file, "fasta"))
        aligned_ref = aligned_records[0].seq
        aligned_query = aligned_records[1].seq
        aligned_target = aligned_records[2].seq
        differences = analyze_differences(aligned_query, aligned_target, aligned_ref)
    finally:
        os.remove(alignment_file) # Clean up alignment file
    
    # 4. Retrieve metadata
    try:
        # Match the ID up to the first space to handle diverse metadata files
        match_meta = metadata_df.loc[closest_match_id.split()[0]]
        species = match_meta.get("Species", "-")
        opsin_family = match_meta.get("Opsin_Family", "-")
        lmax = match_meta.get("Lambda_Max", "-")
        accession = match_meta.get("Accession", closest_match_id)
    except KeyError:
        species, opsin_family, lmax, accession = "-", "-", "-", closest_match_id

    # 5. Assemble results
    result_data = {
        "query_id": query_id,
        "status": "Success",
        "closest_match_id": closest_match_id,
        "match_species": species,
        "match_opsin_family": opsin_family,
        "match_accession": accession,
        "match_lambda_max": lmax,
        **blast_metrics,
        "differences_from_match": differences
    }
    return result_data


# --- Analysis Module Function ---
def run_blastp_analysis(query_sequences, query_ids, opsin_database, opsin_db_fasta, opsin_db_meta, 
                 output_file="sequence_analysis_report.csv", 
                 ref_seq_id="bovine", reffile=None, wrk_dir=".", n_jobs=1):
    """
    Main function to run the complete BLAST, alignment, and analysis pipeline.
    This function is designed to be called from other Python scripts.
    """
    # --- Setup ---
    temp_dir = Path(wrk_dir) / "temp_align_files"
    temp_dir.mkdir(exist_ok=True)
    
    try:
        delimiter = "\t" if '.tsv' in opsin_db_meta else ","
        metadata_df = pd.read_csv(opsin_db_meta, delimiter=delimiter, index_col=0)
    except FileNotFoundError:
        print(f"ERROR: Metadata file not found at {opsin_db_meta}")
        return

    # --- Step 1: Check cache for existing results ---
    cache_path = get_cache_path(opsin_database, wrk_dir)
    blast_cache = load_blast_cache(cache_path)
    
    cached_results = {}
    uncached_records = []
    query_id_to_seq = {qid: seq for qid, seq in zip(query_ids, query_sequences)}

    for qid, seq in zip(query_ids, query_sequences):
        
        qid = qid.replace('\n', ' ').replace('\r', '').strip()

        if seq in blast_cache:
            result = blast_cache[seq].copy()
            result['query_id'] = qid # Ensure query_id matches current run
            cached_results[qid] = result
        else:
            record = SeqRecord(Seq(seq), id=qid, description="")
            uncached_records.append(record)
    
    print(f"INFO: Found {len(cached_results)} sequences in cache. Processing {len(uncached_records)} new sequences.")

    # --- Step 2: Process uncached sequences ---
    newly_processed_results = []
    if uncached_records:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, suffix=".fasta") as tfile:
            temp_query_path = tfile.name
            SeqIO.write(uncached_records, tfile, 'fasta')

        blast_output = run_blastp_bulk(temp_query_path, opsin_database)
        os.remove(temp_query_path)

        if blast_output:
            best_hits_df = filter_best_blast_hits(blast_output)
            if not best_hits_df.empty:
                uncached_query_records = {rec.id: rec for rec in uncached_records}
                tasks_args = []
                for _, row in best_hits_df.iterrows():
                    query_id = row['qseqid']
                    if query_id in uncached_query_records:
                        blast_metrics = {
                            "percent_identity": row['pident'],
                            "alignment_length": row['length'],
                            "mismatches": row['mismatch'],
                            "gap_opens": row['gapopen'],
                            "e_value": row['evalue'],
                        }
                        args = (
                            uncached_query_records[query_id], row['sseqid'], blast_metrics,
                            ref_seq_id, reffile, opsin_db_fasta, metadata_df, temp_dir, wrk_dir
                        )
                        tasks_args.append(args)

                print(f"INFO: Starting parwise analysis to closest VPOD match for {len(tasks_args)} sequences...")
                newly_processed_results = Parallel(n_jobs=n_jobs, verbose=5)(
                    delayed(align_and_analyze_one)(*args) for args in tasks_args
                )

                # Update cache with new successful results
                for result in newly_processed_results:
                    if result.get('status') == 'Success':
                        original_sequence = query_id_to_seq[result['query_id']]
                        blast_cache[original_sequence] = result
            
            # Identify and log uncached sequences that had no BLAST hits
            queries_with_hits = set(best_hits_df['qseqid'])
            all_uncached_ids = {rec.id for rec in uncached_records}
            queries_no_hits = all_uncached_ids - queries_with_hits
            for q_id in queries_no_hits:
                newly_processed_results.append({"query_id": q_id, "status": "No BLAST hits found"})

    # --- Step 3: Finalization ---
    save_blast_cache(blast_cache, cache_path)
    
    all_results_dict = cached_results.copy()
    for result in newly_processed_results:
        all_results_dict[result['query_id']] = result
        
    final_ordered_list = [all_results_dict.get(qid, {'query_id': qid, 'status': 'Processing Error'}) for qid in query_ids]
    
    if final_ordered_list:
        results_df = pd.DataFrame(final_ordered_list)
        results_df.to_csv(output_file, index=False)
        print(f"\n✅ Analysis complete. BLASTp analysis results saved to '{output_file}'\n")
        return results_df
    else:
        print("\n❌ Analysis finished, but no results were generated.")
        
    try:
        if not any(temp_dir.iterdir()):
             os.rmdir(temp_dir)
    except OSError:
        print(f"INFO: Temporary files are in '{temp_dir}'. You may delete this folder.")


# --- Main Execution (for command-line use) ---
def main():
    parser = argparse.ArgumentParser(
        description="""
        BLASTp, align, and analyze protein sequences against a reference database.
        Caches results to avoid re-computing for previously seen sequences.
        NOTE: This script requires 'joblib' and 'biopython'. 
        Install them using: pip install joblib biopython
        """,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("query_file", help="Path to the input FASTA file with one or more query sequences.")
    parser.add_argument("opsin_database", help="Path to the BLAST database (e.g., created with makeblastdb).")
    parser.add_argument("opsin_db_fasta", help="Path to the FASTA file corresponding to the BLAST database.")
    parser.add_argument("opsin_db_meta", help="Path to the metadata CSV/TSV file for the database.")
    parser.add_argument("-o", "--output_file", default="sequence_analysis_report.csv", help="Name of the output CSV file.")
    parser.add_argument("-r", "--ref_seq_id", default="bovine", choices=["bovine", "squid", "microbe", "custom"], help="Reference sequence ID to use for position mapping.")
    parser.add_argument("-rf", "--reffile", help="Path to a custom reference FASTA file (required if --ref_seq_id is 'custom').")
    parser.add_argument("-wd", "--wrk_dir", default=".", help="Working directory for temporary files.")
    parser.add_argument("-j", "--n_jobs", type=int, default=1, help="Number of parallel jobs for the alignment step. Use -1 to use all available CPU cores.")

    args = parser.parse_args()

    # Parse query file into lists for the main function
    try:
        query_records = list(SeqIO.parse(args.query_file, "fasta"))
        query_ids = [rec.id for rec in query_records]
        query_sequences = [str(rec.seq) for rec in query_records]
    except FileNotFoundError:
        print(f"ERROR: Query file not found at {args.query_file}")
        return
    
    if not query_records:
        print(f"ERROR: No sequences found in query file: {args.query_file}")
        return

    run_blastp_analysis(
        query_sequences=query_sequences,
        query_ids=query_ids,
        opsin_database=args.opsin_database,
        opsin_db_fasta=args.opsin_db_fasta,
        opsin_db_meta=args.opsin_db_meta,
        output_file=args.output_file,
        ref_seq_id=args.ref_seq_id,
        reffile=args.reffile,
        wrk_dir=args.wrk_dir,
        n_jobs=args.n_jobs
    )

if __name__ == "__main__":
    main()