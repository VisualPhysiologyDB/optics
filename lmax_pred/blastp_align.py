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
    print(closest_match_id)

    return closest_match_id, blast_metrics

# --- Alignment ---
def align_sequences(query_seq, target_seq, reference_seq):
    """Aligns query, target, and reference sequences using MAFFT"""
    #sequences = [query_seq, target_seq, reference_seq]
    #aligner = MultipleSeqAlignment(sequences)
    #alignment = aligner.align()
    temp_seqs = '/home/PIA/galaxy/tools/optics/tmp/blastp_temp_seqs.fasta'  
    with open(temp_seqs, "w") as temp_file:  # Key change
        temp_file.write(f'{reference_seq}\n')
        temp_file.write(f'{query_seq}\n')
        temp_file.write(f'{target_seq}\n')
    #with open(temp_seqs, "r") as temp_file:
        #for lines in temp_file:
            #print(lines)
    new_ali = '/home/PIA/galaxy/tools/optics/tmp/blastp_temp_ali.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    mafft_exe ='mafft' #change to your own directory for mafft.bat or mafft execution file
    cmd = [mafft_exe,'--auto',temp_seqs]
    with open(new_ali, 'w') as f:
        try:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
            print('Alignment Successful!')
        except subprocess.CalledProcessError as e:
            raise Exception(f'MAFFT alignment failed.\n{e.stderr.decode()}')
    return new_ali

# --- Extract sequences ---
def extract_from_fasta(fasta_file, seq_id):
    """Retrieves a sequence from a FASTA file by its ID."""
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == seq_id:
            seq = ">" + str(record.id).replace(' ','_') + "\n" + str(record.seq)
            #print(f"The Sequence is: {seq}") 
            return (seq)  # Return the sequence as a string


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
    try:
        query_seq_list = query_seq.split('\n')[1:]
        query_seq = ''
        closest_match_seq_list = closest_match_seq.split('\n')[1:]
        closest_match_seq = ''
        reference_seq_list = reference_seq.split('\n')[1:]
        reference_seq = ''

        for x in range(len(query_seq_list)):
            query_seq+=query_seq_list[x]
            closest_match_seq+=closest_match_seq_list[x]
            reference_seq+=reference_seq_list[x]

        print(f'{query_seq}\n{closest_match_seq}\n{reference_seq}')
        differences = []
        bovine_pos = 1  # Initialize the bovine position counter
        for i, (q, c, r) in enumerate(zip(query_seq, closest_match_seq, reference_seq)):
            if q != c and r != '-':
                differences.append((i, q, c, bovine_pos))
                bovine_pos += 1 # Increment bovine position counter only if not a gap 
            elif q != c and r == '-':
                differences.append((i, q, c, 'N/A'))
            elif q == c and r == '-':
                pass
            else:
                bovine_pos+=1 
        return differences
    except:
        raise Exception

# --- Main script ---
def seq_sim_report(query_file, name, ref_seq_id, opsin_database, opsin_db_fasta, opsin_db_meta, ouput_file):
    # BLAST search
    closest_match_id, blast_metrics = run_blastp(query_file, opsin_database)
    
    # Extract sequences
    with open(query_file, 'r') as q:
        for lines in q:
            if '>' in lines:
                name = lines.replace('>','').split(' ')[0]
                print(name)
                break
            else:
                pass
    query_record = extract_from_fasta(query_file, seq_id=name) # Replace with query ID
    print(f"Query record is still: {query_record}")
    closest_match_record = extract_from_fasta(opsin_db_fasta, seq_id=closest_match_id)
    if ref_seq_id == "Bovine":
        reference_record = '>Bovine\nMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA'
    else:
        reference_record = '>Squid\nMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA'
    #print(reference_record)

    # Alignment
    alignment = align_sequences(query_record, closest_match_record, reference_record)

    with open(alignment, 'r') as f:
        aligned_sequences = []
        i = 0
        line_count = 0
        entry = ""
        lines = f.readlines()
        #print(lines) #to check that the alignment actually worked
        num_lines = len(lines)
        #print(num_lines)

        for line in lines:
            if '>' in line:
                if i == 1:
                    aligned_sequences.append(entry)
                    entry = ""
                    entry += f'{line.strip()}\n'
                    #print(sequences)
                    line_count+=1
                else:
                    entry += f'{line.strip()}\n'
                    i+=1
                    line_count+=1
            else:
                entry += line.strip()
                line_count+=1
                if line_count >= num_lines:
                        aligned_sequences.append(entry)
                        print(aligned_sequences)
#print(sequences)
    # Difference analysis with bovine mapping
    reference_seq = str(aligned_sequences[0])
    query_seq = str(aligned_sequences[1])
    closest_match_seq = str(aligned_sequences[2])
    differences = analyze_differences(query_seq, closest_match_seq, reference_seq)

    # Retrieving Phenotype Meta-data for the closest match sequence
    # --- Load metadata ---
    metadata_df = pd.read_csv(opsin_db_meta, delimiter="\t", index_col=0)
    closest_match_data = metadata_df.loc[closest_match_id]
    species = closest_match_data["Species"]
    opsin_family = closest_match_data["Opsin_Family"]
    lmax = closest_match_data["Lambda_Max"]
    accession = closest_match_data["Accession"]

    # Report generation
    with open(ouput_file, 'a+') as f:
        f.write(f"Closest Match to Query Sequence - {name}:\n")
        f.write(f"Species\tOpsin Family\tAccession\tλmax\n")
        f.write(f"{species}\t{opsin_family}\t{accession}\t{lmax}\n")
        f.write("*** BLAST Metrics ***\n")
        for metric, value in blast_metrics.items():
            f.write(f"{metric.capitalize()}: {value}\n")
        f.write("Alignment:\n")
        with open(alignment, 'r') as infile:
            for line in infile:
                f.write(line)
        f.write("\nDifferences between query and closest match (with Bovine position):\n")

        for pos, query_aa, closest_aa, bovine_pos in differences:
         f.write(f"Position: {bovine_pos}\tQuery: {query_aa}\tClosest Match: {closest_aa}\n")
        f.write("\n")

    #print(f"Closest Match to Query Sequence: {accession}\n")
    #print(f"Species\tOpsin Family\tλmax")
    #print(f"{species}\t{opsin_family}\t{lmax}")
    #print("*** BLAST Metrics ***")
    #for metric, value in blast_metrics.items():
        #print(f"{metric.capitalize()}: {value}\n")

    #print("\nAlignment:\n")
    #print(alignment)

    #print("\nDifferences between query and closest match (with bovine position):")
    #for pos, query_aa, closest_aa, bovine_pos in differences:
        #print(f"Position (Bovine): {bovine_pos}, Query: {query_aa}, Closest Match: {closest_aa}\n")

    