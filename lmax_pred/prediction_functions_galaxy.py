# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import pandas as pd
import tempfile
import argparse
import os
import csv  # For CSV export
import time

def process_sequence(sequence, selected_model):
    data_dir = "/home/PIA/galaxy/tools/optics/data"
    model_datasets = {
    "Whole Dataset Model": f"{data_dir}/whole_dataset.fasta",
    "Wild-Type Model": f"{data_dir}/wild_type_opsins.fasta",
    "Vertebrate Model": f"{data_dir}/vertebrate_opsins.fasta",
    "Invertebrate Model": f"{data_dir}/invertebrate_opsins.fasta",
    "Rod Model": f"{data_dir}/rod_opsins.fasta",
    }   

    model_dir = "/home/PIA/galaxy/tools/optics/models"
    model_directories = {
        "Whole Dataset Model": f"{model_dir}/wds_gbc.pkl",
        "Wild-Type Model": f"{model_dir}/wt_gbc.pkl",
        "Vertebrate Model": f"{model_dir}/vert_gbc.pkl",
        "Invertebrate Model": f"{model_dir}/invert_lgbm.pkl",
        "Rod Model": f"{model_dir}/rod_gbc.pkl",
    }

    if sequence == None:
        return ('Error: No sequence given')
    #print(sequence)
    if selected_model == None:
        return ('Error: No model selected')
    
    alignment_data = model_datasets[selected_model]
    #print(alignment_data)
    selected_model = model_directories[selected_model]

    temp_seq = "/home/PIA/galaxy/tools/optics/tmp/temp_seq.fasta"
    with open(temp_seq, "w") as temp_file:  # Key change
        if '>' in sequence:
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence) # Write your data to the file object

    with open(temp_seq, 'r') as f:
        for lines in f:
            print(lines)
        

    new_ali = '/home/PIA/galaxy/tools/optics/tmp/temp_ali.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    mafft_exe ='mafft' #change to your own directory for mafft.bat or mafft execution file
    seq_type = 'aa'

    cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
    with open(new_ali, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

    new_seq_test = read_data(new_ali, seq_type = seq_type, is_main=True, gap_threshold=0.6)
    ref_copy = read_data(alignment_data, seq_type = seq_type, is_main=True, gap_threshold=0.6)
    last_seq = ref_copy.shape[0]
    new_seq_test = new_seq_test.iloc[last_seq:].copy()
    #print(new_seq_test)

    # ... (Load the selected model and make a prediction)
    load_top_mod = load_obj(selected_model)
    prediction = load_top_mod.predict(new_seq_test)
    
    return(prediction[0])
 

def process_sequences_from_file(file,selected_model):
    if file == None:
        return ('Error: No file given')

    with open(file, 'r') as f:
        sequences = []
        names = []
        i = 0
        line_count = 0
        entry = ""
        lines = f.readlines()
        num_lines = len(lines)
        #print(num_lines)

        for line in lines:
            if '>' in line:
                if i == 1:
                    names.append(line.replace('>','').strip())
                    sequences.append(entry)
                    entry = ""
                    entry += line
                    #print(sequences)
                    line_count+=1
                else:
                    names.append(line.replace('>','').strip())
                    entry += line
                    i+=1
                    line_count+=1
            else:
                entry += line
                line_count+=1
                if line_count >= num_lines:
                     sequences.append(entry)
    #print(sequences)
    #print(names)
                     
    predictions = []
    for seq in sequences:
        print(seq)
        prediction = process_sequence(seq, selected_model)  # Process each sequence
        predictions.append(prediction)
    #print(predictions)

    return(names,predictions)

def main():
    
    model_dir = "/home/PIA/galaxy/tools/optics/models"
    model_directories = {
        "Whole Dataset Model": f"{model_dir}/wds_gbc.pkl",
        "Wild-Type Model": f"{model_dir}/models/wt_gbc.pkl",
        "Vertebrate Model": f"{model_dir}/vert_gbc.pkl",
        "Invertebrate Model": f"{model_dir}/invert_lgbm.pkl",
        "Rod Model": f"{model_dir}/rod_gbc.pkl",
    }


    parser = argparse.ArgumentParser(description="Process sequences using a selected model")
    parser.add_argument("input", help="Either a single sequence or a path to a FASTA file", type=str)
    parser.add_argument("output", help="Name for output file", type=str, default = 'output.txt')
    parser.add_argument("-m", "--model", help="Model to use for prediction", 
                        choices=list(model_directories.keys()), required=True)  

    args = parser.parse_args()

    if os.path.isfile(args.input):
        names, predictions = process_sequences_from_file(args.input, args.model)
        with open(args.output, 'w') as f:
            i = 0
            while i in range(len(names)):
                if i == 0:
                    f.write('Names\tPredictions\n')
                f.write(f"{names[i]}:\t{predictions[i]}\n")
                print(f"{names[i]}:\t{predictions[i]}\n")
                i+=1
    else:
        prediction = process_sequence(args.input, args.model)
        print(prediction) 

if __name__ == "__main__":
    main()