# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import pandas as pd
import tempfile
import os
import csv  # For CSV export

def write_sequence_to_temp(sequence, temp_seq):
    print(sequence)
    with open(temp_seq, 'w') as g:
        
        if '>' in sequence:
            g.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            g.write(sequence)


def process_sequence(sequence, selected_model, is_multi_seq = False):
    model_datasets = {
    "Whole Dataset Model": "./data/whole_dataset.fasta",
    "Vertebrate Model": "./data/vertebrate_opsins.fasta",
    "Rod Model": "./data/rod_opsins.fasta",

    }   

    model_directories = {
        "Whole Dataset Model": "./models/wds_gbc.pkl",
        "Vertebrate Model": "./models/vert_gbc.pkl",
        "Rod Model": "./models/rod_gbc.pkl",


    }


    if sequence == None:
        return ('Error: No sequence given')
    #print(sequence)
    if selected_model == None:
        return ('Error: No model selected')
    
    alignment_data = model_datasets[selected_model]
    selected_model = model_directories[selected_model]
    temp_seq = 'temp.fasta'
    write_sequence_to_temp(sequence = sequence, temp_seq = temp_seq)

    # ... (Perform alignment using MAFFT with alignment_data)
    mafft_exe = 'C:/Users/safra/mafft-win/mafft.bat' #change to your own directory for mafft.bat file
    seq_type = 'aa'
    new_ali = 'temp_ali.fasta'
    cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data , '>', new_ali]
    aligner = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out_put = aligner.communicate()[0].decode('utf8')
    
    new_seq_test = read_data(new_ali, seq_type = seq_type, is_main=True, gap_threshold=0.6)
    print(new_seq_test)
    ref_copy = read_data(alignment_data, seq_type = seq_type, is_main=True, gap_threshold=0.6)
    last_seq = ref_copy.shape[0]
    new_seq_test = new_seq_test.iloc[last_seq:].copy()
    print(new_seq_test)

    # ... (Load the selected model and make a prediction)
    load_top_mod = load_obj(selected_model)
    prediction = load_top_mod.predict(new_seq_test)
    
    return(prediction[0])
 


def process_sequences_from_file(file,selected_model):
    if file == None:
        return ('Error: No file given')
    
    temp_dir = os.path.join(os.getcwd(), 'tmp')  
    with tempfile.NamedTemporaryFile(dir=temp_dir, delete=False) as temp_file:  # Key change
        file.save(temp_file.name) 


        with open(temp_file.name, 'r') as f:
            sequences = []
            names = []
            i = 0
            entry = ""
            for lines in f:
                if '>' in lines:
                    if i == 1:
                        names.append(lines.replace('>','').strip())
                        sequences.append(entry)
                        entry = ""
                        entry += lines
                        print(sequences)
                    else:
                        names.append(lines.replace('>','').strip())
                        entry += lines
                        i+=1
                else:
                    entry += lines

        predictions = []
        i = 0
        #is_multi_seq = True
        for seq in sequences:
            print(seq)
            prediction = process_sequence(seq, selected_model)  # Process each sequence
            predictions.append(prediction)
            i+=1
        #is_multi_seq = False
        print(predictions)

        return(names,predictions)