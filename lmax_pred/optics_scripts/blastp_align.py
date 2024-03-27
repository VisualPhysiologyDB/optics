#blastp_align.py 

import subprocess
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# --- BLAST search ---
def run_blastp(query_file, database):
    """Performs BLASTp search and returns the top hit."""
    blastp_path = "blastp"  # Replace with the actual path to blastp
    blast_cmd = [blastp_path, "-query", query_file, "-db", database, "-outfmt", "6"]
    blast_result = subprocess.run(blast_cmd, capture_output=True, text=True).stdout
    blast_metrics = parse_blast_result(blast_result)
    closest_match_id = blast_result.splitlines()[1].split("\t")[1]

    return closest_match_id, blast_metrics

# --- Alignment ---
def align_sequences(query_seq, target_seq, reference_seq):
    """Aligns query, target, and reference sequences using Clustal Omega."""
    sequences = [query_seq, target_seq, reference_seq]
    aligner = MultipleSeqAlignment(sequences)
    alignment = aligner.align()
    return alignment

# --- Extract sequences ---
def extract_from_fasta(fasta_file, seq_id):
    """Retrieves a sequence from a FASTA file by its ID."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == seq_id:
            return record

# --- BLAST result parsing ---
def parse_blast_result(blast_output):
    """Extracts key metrics from BLAST tabular output (outfmt 6)."""
    lines = blast_output.splitlines()
    top_hit_data = lines[1].split("\t")  # Data from the second line
    return {
        "percent_identity": float(top_hit_data[2]),
        "alignment_length": int(top_hit_data[3]),
        "mismatches": int(top_hit_data[4]),
        "gap_opens": int(top_hit_data[5]),
        "e_value": float(top_hit_data[10]),
    }

# --- Difference analysis with bovine position mapping ---
def analyze_differences(query_seq, closest_match_seq, reference_seq):
    differences = []
    bovine_pos = 1  # Initialize the bovine position counter
    for i, (q, c, r) in enumerate(zip(query_seq, closest_match_seq, reference_seq)):
        if q != c:
            differences.append((i, q, c, bovine_pos))
        if r != '-':  # Increment bovine position counter only if not a gap 
            bovine_pos += 1
    return differences

# --- Main script ---
def seq_sim_report(query_file, opsin_database, opsin_db_fasta, opsin_db_meta, ouput_file):
    # BLAST search
    closest_match_id, blast_metrics = parse_blast_result(query_file, opsin_database)
    
    # Extract sequences
    query_record = extract_from_fasta(query_file, seq_id="your_query_id")  # Replace with query ID
    closest_match_record = extract_from_fasta(opsin_database, closest_match_id)
    reference_record = extract_from_fasta(opsin_db_fasta, "Bovine")

    # Alignment
    alignment = align_sequences(query_record, closest_match_record, reference_record)

    # Difference analysis with bovine mapping
    query_seq = str(alignment[0].seq)
    closest_match_seq = str(alignment[1].seq)
    reference_seq = str(alignment[2].seq)
    differences = analyze_differences(query_seq, closest_match_seq, reference_seq)

    # Retrieving Phenotype Meta-data for the closest match sequence
    # --- Load metadata ---
    metadata_df = pd.read_csv(opsin_db_meta, delimiter="\t", index_col="ID")
    closest_match_data = metadata_df.loc[closest_match_id]
    species = closest_match_data["Species"]
    opsin_family = closest_match_data["Opsin_Family"]
    lmax = closest_match_data["Lambda_Max"]
    accession = closest_match_data["Accession"]

    # Report generation
    with open(ouput_file, 'a+') as f:
        f.write(f"Closest Match to Query Sequence: {accession}\n")
        f.write(f"Species\tOpsin Family\tλmax")
        f.write(f"{species}\t{opsin_family}\t{lmax}")
        f.write("*** BLAST Metrics ***")
        for metric, value in blast_metrics.items():
            f.write(f"{metric.capitalize()}: {value}\n")
        f.write("\nAlignment:")
        f.write(alignment)
        f.write("\nDifferences between query and closest match (with bovine position):")
        for pos, query_aa, closest_aa, bovine_pos in differences:
         f.write(f"Position (Bovine): {bovine_pos}, Query: {query_aa}, Closest Match: {closest_aa}\n")

    print(f"Closest Match to Query Sequence: {accession}\n")
    print(f"Species\tOpsin Family\tλmax")
    print(f"{species}\t{opsin_family}\t{lmax}")
    print("*** BLAST Metrics ***")
    for metric, value in blast_metrics.items():
        print(f"{metric.capitalize()}: {value}\n")

    print("\nAlignment:")
    print(alignment)

    print("\nDifferences between query and closest match (with bovine position):")
    for pos, query_aa, closest_aa, bovine_pos in differences:
        print(f"Position (Bovine): {bovine_pos}, Query: {query_aa}, Closest Match: {closest_aa}\n")

    