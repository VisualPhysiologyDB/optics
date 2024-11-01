#blastp_align.py 

import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
import random
import pandas as pd
# --- BLAST search ---
def run_blastp(query_file, database):
    """Performs BLASTp search and returns the top hit."""
    blastp_path = "blastp"  # Replace with the actual path to blastp
    blast_cmd = [blastp_path, "-query", query_file, "-db", database, "-outfmt", "6"]
    run_blastp = subprocess.run(blast_cmd, capture_output=True, text=True)

    # Access standard output (stdout)
    blast_result = run_blastp.stdout

    # Access standard error (stderr)
    stderr_data = run_blastp.stderr

    blast_metrics = parse_blast_result(blast_result)
    try:    
        closest_match_id = blast_result.splitlines()[1].split("\t")[1]
    except:
        closest_match_id = '-'
        print(f'This is the blast result which raised the following error:\n{blast_result}\n')
        print(f'This is the error message from blastp:\n{stderr_data}\n')

    return closest_match_id, blast_metrics

# --- Alignment ---
def align_sequences(query_seq, target_seq, reference_seq):
    """Aligns query, target, and reference sequences using MAFFT"""
    #get working directory to make temp files
    wrk_dir = os.getcwd().replace('\\','/')
    #print(wrk_dir)
    random_number = str(random.randint(1, 10000))
    temp_seqs = f'{wrk_dir}/tmp/blastp_temp_seqs_{random_number}.fasta'  
    with open(temp_seqs, "w") as temp_file:  # Key change
        temp_file.write(f'{reference_seq}\n')
        temp_file.write(f'{query_seq}\n')
        temp_file.write(f'{target_seq}\n')
    #with open(temp_seqs, "r") as temp_file:
        #for lines in temp_file:
            #print(lines)
    new_ali = f'{wrk_dir}/tmp/blastp_temp_ali_{random_number}.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    try:
        mafft_exe ='mafft' #change to your own directory for mafft.bat or mafft execution file
        cmd = [mafft_exe,'--auto',temp_seqs]
        with open(new_ali, 'w') as f:
            #try:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
                #print('Alignment Successful!')
            #except subprocess.CalledProcessError as e:
            #    raise Exception(f'MAFFT alignment failed.\n{e.stderr.decode()}')
    except Exception as e:
        #print(f'Exception {e} has occured. Tring Windows MAFFT')
        try:
            mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-win/mafft.bat'
            cmd = [mafft_exe,'--auto',temp_seqs]
            with open(new_ali, 'w') as f:
                #try:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
                    #print('Alignment Successful!')
                #except subprocess.CalledProcessError as e:
                #    raise Exception(f'MAFFT alignment failed.\n{e.stderr.decode()}')`
        except Exception as e: 
            #print(f'Exception {e} has occured. Tring Max MAFFT')
            try:
                mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-mac/mafft.bat'
                cmd = [mafft_exe,'--auto',temp_seqs]
                with open(new_ali, 'w') as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
                    #print('Alignment Successful!')
            except subprocess.CalledProcessError as e:
                raise Exception(f'MAFFT alignment failed for all three options. Check your FASTA file to make sure it is formatted correctly\n{e.stderr.decode()}')
    try:
        os.remove(temp_seqs)
    except FileNotFoundError:
        raise Exception("File does not exist")
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
    try:
        top_hit_data = lines[1].split("\t")  # Data from the second line
        return {
        "percent_identity": float(top_hit_data[2]),
        "alignment_length": int(top_hit_data[3]),
        "mismatches": int(top_hit_data[4]),
        "gap_opens": int(top_hit_data[5]),
        "e_value": float(top_hit_data[10]),
        }
    except:
        try:
            top_hit_data = lines[0].split("\t")  # Data from the second line
            return {
            "percent_identity": float(top_hit_data[2]),
            "alignment_length": int(top_hit_data[3]),
            "mismatches": int(top_hit_data[4]),
            "gap_opens": int(top_hit_data[5]),
            "e_value": float(top_hit_data[10]),
            }
        except:
            print("\nAn error occured during blastp but moving on with blank '-' entries for this sequence\n")
            #print(lines)
            return None


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

        #print(f'{query_seq}\n{closest_match_seq}\n{reference_seq}')
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
def seq_sim_report(query_file, name, ref_seq_id, opsin_database, opsin_db_fasta, opsin_db_meta, ouput_file, reffile):
    # BLAST search
    closest_match_id, blast_metrics = run_blastp(query_file, opsin_database)
    if closest_match_id == '-':
        return('blastp unsuccessful')
    # Extract sequences
    query_record = ''
    with open(query_file, 'r') as q:
        for lines in q:
            entry = lines.strip().replace('\n','') + '\n'
            query_record+=entry
    #print(f'This is the query record:\n{query_record}\n')
    
    closest_match_record = extract_from_fasta(opsin_db_fasta, seq_id=closest_match_id)
    
    if ref_seq_id == "bovine":
        reference_record = '>Bovine\nMNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA'
    elif ref_seq_id == "squid":
        reference_record = '>Squid\nMGRDLRDNETWWYNPSIVVHPHWREFDQVPDAVYYSLGIFIGICGIIGCGGNGIVIYLFTKTKSLQTPANMFIINLAFSDFTFSLVNGFPLMTISCFLKKWIFGFAACKVYGFIGGIFGFMSIMTMAMISIDRYNVIGRPMAASKKMSHRRAFIMIIFVWLWSVLWAIGPIFGWGAYTLEGVLCNCSFDYISRDSTTRSNILCMFILGFFGPILIIFFCYFNIVMSVSNHEKEMAAMAKRLNAKELRKAQAGANAEMRLAKISIVIVSQFLLSWSPYAVVALLAQFGPLEWVTPYAAQLPVMFAKASAIHNPMIYSVSHPKFREAISQTFPWVLTCCQFDDKETEDDKDAETEIPAGESSDAAPSADAAQMKEMMAMMQKMQQQQAAYPPQGYAPPPQGYPPQGYPPQGYPPQGYPPQGYPPPPQGAPPQGAPPAAPPQGVDNQAYQA'
    elif ref_seq_id == "microbe":
        reference_record = '>AcetR1\nMSNPNPFQTTLGTDAQWVVFAVMALAAIVFSIAVQFRPLPLRLTYYVNIAICTIAATAYYAMAVNGGDNKPTAGTGADERQVIYARYIDWVFTTPLLLLDLVLLTNMPATMIAWIMGADIAMIAFGIIGAFTVGSYKWFYFVVGCIMLAVLAWGMINPIFKEELQKHKEYTGAYTTLLIYLIVLWVIYPIVWGLGAGGHIIGVDVEIIAMGILDLLAKPLYAIGVLITVEVVYGKLGKEEAQPLTA'
    else:
        #this is for extracting the custom reference sequence from a fasta file of presumable 1 sequence.
        reference_record = ''
        with open(reffile, 'r') as f:
            lines = f.readlines()
            i = 0
            for line in lines:
                if '>' in line:
                    if i == 0:
                        entry = line.strip().replace('\n','') + '\n'
                        reference_record+=entry
                        ref_seq_id = line.strip().replace('\n','').replace('>','').split(' ')[0]
                        i+=1
                    else:
                        break
                else:
                    reference_record+=line
    #print(f'This is the reference record:\n{reference_record}')

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
                        #print(aligned_sequences)
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
    with open(ouput_file, 'a+', encoding="utf-8") as f:
        f.write(f"Closest Match to Query Sequence - {name}:\n")
        f.write(f"VPOD_ID\tSpecies\tOpsin Family\tAccession\tλmax\n")
        f.write(f"{closest_match_id}\t{species}\t{opsin_family}\t{accession}\t{lmax}\n")
        f.write("*** BLAST Metrics ***\n")
        for metric, value in blast_metrics.items():
            if 'percent' in metric:
                percent_iden = value
            f.write(f"{metric.capitalize()}: {value}\n")
        f.write("Alignment:\n")
        with open(alignment, 'r') as infile:
            for line in infile:
                f.write(line)
        f.write(f"\nDifferences between query and closest match (with {ref_seq_id} position):\n")

        for pos, query_aa, closest_aa, bovine_pos in differences:
         f.write(f"Position: {bovine_pos}\tQuery: {query_aa}\tClosest Match: {closest_aa}\n")
        f.write("\n")
    try:
        os.remove(alignment)
    except FileNotFoundError:
        raise Exception("File does not exist")
    return(percent_iden)
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

    