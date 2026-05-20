"""
A script for performing in-silico mutagenesis on protein sequences.

This script can generate mutated protein sequences based on a set of specified
mutations against a wild-type sequence. It can fetch sequences from NCBI or use
a predefined set of ancestral sequences.

The main functionalities are:
1.  Generate all combinatorial mutants from a list of mutations against a
    wild-type accession.
2.  Take a predefined list of mutant accessions (in the format
    'ACCESSION_M1,M2,...') and generate the corresponding mutated protein
    sequences.

The script aligns the target sequence to a reference sequence (e.g., Bos taurus)
to ensure mutation positions are correct, even with insertions or deletions.

Command-Line Usage Examples:

1.  Generate all combinations from a list of mutations:
    python mutagenesis.py \
        --wt_accession AncBovine \
        --mutations "A116S,S119A,G121A" \
        --output_file combined_mutants.fasta \
        --reference_accession NM_001014890.2

2.  Generate sequences for specific mutants from a file:
    python mutagenesis.py \
        --mutant_file my_mutants.txt \
        --output_file specific_mutants.fasta

3.  Generate a sequence for a single mutant:
    python mutagenesis.py \
        --mutant_accession "AncBovine_A116S,G121A" \
        --output_file single_mutant.fasta
"""

import re
import argparse
import itertools
import pandas as pd
import os
import json
from Bio import Entrez, SeqIO
from skbio import Protein
from skbio.alignment import global_pairwise_align_protein
from Bio.Align import substitution_matrices

# --- Global Configuration ---

# It's good practice for the user of the script to provide their email
# for NCBI Entrez.
ENTREZ_EMAIL = 'your.email@example.com'
Entrez.email = ENTREZ_EMAIL

cache_dir = f"./data/cached_seqs"
os.makedirs(cache_dir, exist_ok=True)
cache_file = f"{cache_dir}/wt_opsin_seq_dict.json"
try:
    with open(cache_file, 'r') as f:
        CACHED_SEQUENCES = json.load(f)
    original_key_count = len(CACHED_SEQUENCES.keys())
    print('\nWT sequence cache file successfully loaded.\n')
except (json.JSONDecodeError, FileNotFoundError):
    CACHED_SEQUENCES = {}
    original_key_count = 0
    print('\nCache file not found or invalid. A new cache will be created.\n')
    

def extract_fasta_entries(file):
    with open(file, 'r') as f:
        sequences = []
        names = []
        first_entry = True
        line_count = 0
        entry = ""
        lines = f.readlines()
        num_lines = len(lines)

        for line in lines:
            if '>' in line:
                if first_entry == False:
                    # Append completed entry
                    sequences.append(entry)
                    # Append name of new entry
                    names.append(line.replace('>','').strip().replace(' ','_').replace('\n',''))
                    # Restart with new entry
                    entry = ""
                    line_count+=1
                else:
                    # First entry - must declare entry outside the loop, so this is just a neccessary artifact
                    names.append(line.replace('>','').strip().replace(' ','_').replace('\n',''))
                    first_entry = False
                    line_count+=1
            else:
                # This should be adding the new lines of aa data 
                entry += line.strip().replace('\n','')
                line_count+=1
                if line_count >= num_lines:
                     sequences.append(entry)
    return names,sequences

def fetch_protein_sequence(accession, sequence_type='target', db_order=['nucleotide','protein']):
    """
    Fetches a protein sequence from various sources.

    Order of operations:
    1. Checks for manually provided sequence string.
    2. Checks the hardcoded ancestral sequence dictionary (cache).
    3. Queries NCBI databases (iterating through db_order).
    4. If all else fails, prompts the user for a manual sequence string.

    Args:
        accession (str): The accession name, number, or 'manual'.
        sequence_type (str): Type of sequence ('reference' or 'target').
        db_order (list): List of databases to query.

    Returns:
        str: The fetched amino acid sequence.
    """
    if db_order is None:
        db_order = ['nucleotide', 'protein']

    if accession.lower() == "manual":
        return input(f"Enter {sequence_type.capitalize()} Sequence: ")

    for db in db_order:
        try:
            print(f"Fetching '{accession}' from NCBI {db.capitalize()}...")
            handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "gb")
            handle.close()
            
            if db == "nucleotide":
                for feature in record.features:
                    if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                        return feature.qualifiers['translation'][0]
                raise ValueError(f"No CDS with translation found in {db} record")
            
            elif db == "protein":
                return str(record.seq)
                
        except Exception as e:
            print(f"  - Warning: Could not fetch '{accession}' from {db} database. Error: {e}")
            continue # Attempt next database

    # If all searches fail
    manual_seq = input(f"Please enter the {sequence_type} sequence for '{accession}' manually: ")
    return manual_seq


def get_mutant_combinations(wt_accession, mutations_list):
    """
    Generates a list of mutant accession strings from all combinations of mutations.

    For a wt_accession 'WT' and mutations ['A1B', 'C2D'], it will generate:
    ['WT_A1B', 'WT_C2D', 'WT_A1B,C2D']

    Args:
        wt_accession (str): The accession name of the wild-type sequence.
        mutations_list (list): A list of mutations as strings (e.g., ['A123G', 'F45S']).

    Returns:
        list: A list of formatted mutant accession strings.
    """
    mutant_accessions = []
    for i in range(1, len(mutations_list) + 1):
        # Get all combinations of mutations of length i
        combinations = itertools.combinations(mutations_list, i)
        for combo in combinations:
            # Format the accession string: e.g., "MyProtein_A1B,C2D"
            mut_str = ",".join(combo)
            mutant_accessions.append(f"{wt_accession}_{mut_str}")
    return mutant_accessions


def parse_mutant_accession(mutant_accession):
    """
    Parses a mutant accession string into its components.

    Args:
        mutant_accession (str): The string to parse (e.g., "MyProtein_A1B,C2D").

    Returns:
        tuple: A tuple containing (wild_type_accession, list_of_mutations).
               Returns (mutant_accession, []) if no mutations are found.
    """
    mutant_accession = mutant_accession.strip()
    count = mutant_accession.count('_')
    
    if count == 1 and len(mutant_accession.split('_')[0]) > 3:
        parts = mutant_accession.split('_')
        wt_accession = parts[0]
        mutations_str = parts[1]
        # Handle single or multiple mutations
        mutations = mutations_str.split(',')
        return wt_accession, mutations
    elif count >= 2:
        parts = mutant_accession.split('_')
        wt_accession = ''
        for x in range(len(parts)-1):
            wt_accession+=parts[x]+'_'
        wt_accession = wt_accession[:-1]
        mutations_str = parts[-1]
        # Handle single or multiple mutations
        mutations = mutations_str.split(',')
        return wt_accession, mutations
    else:
        # It's a wild-type sequence with no mutations
        return mutant_accession, []


def get_mutant_seqs(
    mutant_accessions,
    db_order=['nucleotide','protein'],
    wt_seq_file=None,
    output_file=None,
    reference_accession=None,
    output_format='fasta',
    allow_wt=True,
    email='your.email@example.com'
):
    """
    Generates mutated sequences and writes them to an output file.

    For each mutant accession string, this function will:
    1. Parse the wild-type accession and the required mutations.
    2. Fetch the wild-type and reference protein sequences.
    3. Align the WT to the reference to correctly map mutation sites.
    4. Apply the mutations.
    5. Write the final sequence to the output file.

    Args:
        mutant_accessions (list): List of mutant accession strings to process.
        output_file (str): Path to the output file.
        reference_accession (str): Accession of the sequence for numbering.
        output_format (str): 'fasta' or 'tsv'.
        allow_wt (bool): If True, allows processing of sequences with no mutations.
    """
    # First we need to check if the user provided a FASTA file of WT sequences and add those to the cached sequences dictionary
    if wt_seq_file != None:
        wt_accs, sequences = extract_fasta_entries(wt_seq_file)
        for wt_acc, seq in zip(wt_accs,sequences):
          if wt_acc not in CACHED_SEQUENCES.keys():
            CACHED_SEQUENCES[wt_acc] = seq

    if email != 'your.email@example.com':
        Entrez.email = email
    CURRENT_WT = ''
    
    if reference_accession not in CACHED_SEQUENCES.keys():
        print(f"\nFetching reference sequence '{reference_accession}'...")
        CACHED_SEQUENCES[reference_accession] = (fetch_protein_sequence(reference_accession, 'reference', db_order))

    reference_protein = Protein(CACHED_SEQUENCES[reference_accession])

    # Prepare output file
    if output_format == 'tsv':
        with open(output_file, 'w') as f:
            f.write('Accession\tSequence\n')
    else:
        with open(output_file, 'w') as f:
            f.write('') # Create an empty file to overwrite previous runs

    print("\nProcessing mutations...")
    for mutant_acc in mutant_accessions:
        wt_acc, mutations = parse_mutant_accession(mutant_acc)

        if not mutations:
            if allow_wt:
                print(f"\nProcessing Wild-Type: {wt_acc}")
                if wt_acc not in CACHED_SEQUENCES.keys():
                    CACHED_SEQUENCES[wt_acc] = fetch_protein_sequence(wt_acc, 'target', db_order)
                wt_seq = Protein(CACHED_SEQUENCES[wt_acc])
                final_seq = wt_seq
                final_acc = wt_acc
            else:
                print(f"\nSkipping Wild-Type sequence '{wt_acc}' as per settings.")
                continue
        else:
            print(f"\nProcessing: {mutant_acc}")
            if wt_acc not in CACHED_SEQUENCES.keys():
                CACHED_SEQUENCES[wt_acc] = fetch_protein_sequence(wt_acc, 'target', db_order)
            
            wt_protein = Protein(CACHED_SEQUENCES[wt_acc])

            # Align the WT to the reference to get the correct numbering
            if wt_protein != CURRENT_WT:
                substitution_matrix = substitution_matrices.load("BLOSUM62")
                alignment, _, _ = global_pairwise_align_protein(
                    reference_protein, wt_protein,
                    gap_open_penalty=11,
                    gap_extend_penalty=1,
                    substitution_matrix=substitution_matrix
                )
                
                aligned_ref = str(alignment[0])
                aligned_wt = str(alignment[1])
                CURRENT_WT = wt_protein
            
            mutated_seq_list = list(aligned_wt)

            for mutation in mutations:
                # Regex to parse mutations like 'A123G'
                match = re.match(r'([A-Z])(\d+)([A-Z])', mutation, re.IGNORECASE)
                if not match:
                    print(f"  - WARNING: Skipping invalid mutation format '{mutation}'")
                    continue
                
                original_aa, pos, new_aa = match.groups()
                pos = int(pos) - 1  # Convert to 0-based index

                # Find the equivalent position in the aligned sequence
                ref_pos_count = 0
                aligned_idx = -1
                for i, char in enumerate(aligned_ref):
                    if char != '-':
                        if ref_pos_count == pos:
                            aligned_idx = i
                            break
                        ref_pos_count += 1

                if aligned_idx == -1:
                    print(f"  - WARNING: Position {pos+1} not found in reference alignment. Skipping.")
                    continue
                
                # Check if the original AA matches the sequence at that position
                if aligned_wt[aligned_idx].upper() != original_aa.upper():
                    print(f"  - WARNING: Mismatch at position {pos+1}. "
                          f"Expected '{original_aa}', found '{aligned_wt[aligned_idx]}'. "
                          "This can happen if the reference and target are very divergent. Skipping mutation.")
                    continue
                
                # Apply mutation
                mutated_seq_list[aligned_idx] = new_aa

            final_seq = "".join(mutated_seq_list).replace('-', '')
            final_acc = mutant_acc

        # Write to output file
        with open(output_file, 'a+') as f:
            if output_format == 'fasta':
                f.write(f'>{final_acc}\n{final_seq}\n')
            else:
                f.write(f'{final_acc}\t{final_seq}\n')
        print(f"  -> Saved as {final_acc}")

    print(f"\nProcessing complete. Output saved to '{output_file}'.")
    
    if len(CACHED_SEQUENCES.keys()) > original_key_count:
        try:
            with open(cache_file, 'w') as f:
                json.dump(CACHED_SEQUENCES, f, indent=4)
        except Exception as e:
            print(f"Error: Could not save cache file: {e}")


def main():
    """Main function to parse arguments and run the mutagenesis workflow."""
    parser = argparse.ArgumentParser(
        description='A script for in-silico site-directed mutagenesis.',
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
            This script allows two main modes of operation:
            1. Generate all combinations of mutants from a base wild-type sequence.
            (Use --wt_accession and --mutations)
            2. Generate sequences for a pre-defined list of mutants.
            (Use --mutant_file or --mutant_accession)
        """
    )

    # --- Input Modes (Mutually Exclusive) ---
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--mutant_file',
        help='Path to a text file containing mutant accessions, one per line.'
    )
    input_group.add_argument(
        '--mutant_accession',
        help="A single mutant accession string (e.g., 'MyProtein_A123G,F45S')."
    )
    input_group.add_argument(
        '--wt_accession',
        help='Wild-type accession to generate mutants from (used with --mutations).'
    )
    input_group.add_argument(
        '--wt_file',
        help='Path to a FASTA file contaning the wildtype sequence(s)...'
    )
    parser.add_argument(
        '--mutations',
        help='Comma-separated list of mutations for combinatorial generation...'
    )
    parser.add_argument(
        '-o', '--output_file',
        required=True,
        help='Path to the output file for the generated sequences.'
    )
    parser.add_argument(
        '-ra', '--reference_accession',
        default='NM_001014890',
        help='Reference accession for sequence numbering.'
    )
    parser.add_argument(
        '--db_preference',
        choices=['nucleotide', 'protein'],
        default='nucleotide',
        help='Which NCBI database to search first. Defaults to nucleotide, falling back to protein.'
    )
    parser.add_argument(
        '--output_format',
        choices=['fasta', 'tsv'],
        default='fasta',
        help="Format for the output file. Default: 'fasta'."
    )
    parser.add_argument(
        '--no_wt',
        action='store_true',
        help='Flag to prevent processing of wild-type sequences (those with no mutations).'
    )
    parser.add_argument(
        '--email',
        default=ENTREZ_EMAIL,
        help='Your email address for NCBI Entrez queries.'
    )

    args = parser.parse_args()

    # --- Argument Validation ---
    if args.email != 'your.email@example.com':
        Entrez.email = args.email
    else:
        print("Warning: Using default Entrez email. Please provide your own with the --email flag.")

    if args.wt_accession and not args.mutations:
        parser.error('--mutations is required when using --wt_accession.')

    mutant_list = []
    db_order = ['nucleotide', 'protein'] if args.db_preference == 'nucleotide' else ['protein', 'nucleotide']

    # --- Workflow Selection ---
    if args.wt_accession:
        # Mode 1: Generate combinations
        print("Mode: Generating mutant combinations.")
        mutations_list = [m.strip() for m in args.mutations.split(',')]
        mutant_list = get_mutant_combinations(args.wt_accession, mutations_list)
        # Also add the WT to the list if desired
        if not args.no_wt:
            mutant_list.insert(0, args.wt_accession)
        
    elif args.mutant_file:
        # Mode 2a: Process mutants from a file
        print(f"Mode: Processing mutants from file '{args.mutant_file}'.")
        with open(args.mutant_file, 'r') as f:
            mutant_list = [line.strip() for line in f if line.strip()]

    elif args.mutant_accession:
        # Mode 2b: Process a single mutant
        print(f"Mode: Processing single mutant '{args.mutant_accession}'.")
        mutant_list = [args.mutant_accession]
            
    wt_seq_file = args.wt_file if args.wt_file else None

    if mutant_list:
        get_mutant_seqs(
            mutant_accessions=mutant_list,
            db_order=db_order,
            wt_seq_file=wt_seq_file,
            output_file=args.output_file,
            reference_accession=args.reference_accession,
            output_format=args.output_format,
            allow_wt=not args.no_wt
        )
    else:
        print("No mutants to process. Exiting.")


if __name__ == '__main__':
    main()