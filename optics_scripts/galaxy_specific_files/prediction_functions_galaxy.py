# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from fpdf import FPDF
import pandas as pd
import argparse
import os
import random
from optics_scripts.blastp_align import seq_sim_report
from optics_scripts.bootstrap_predictions import calculate_ensemble_CI , plot_predictions_with_CI 

def process_sequence(sequence, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict = None, encoding_method='one_hot'):
    data_dir = "/home/PIA/galaxy/tools/optics/data"
    model_datasets = {
    "Whole Dataset Model": f"{data_dir}/fasta/vpod_1.2/wds_aligned_VPOD_1.2_het.fasta",
    "Wild-Type Model": f"{data_dir}/fasta/vpod_1.2/wt_aligned_VPOD_1.2_het.fasta",
    "Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/vert_aligned_VPOD_1.2_het.fasta",
    "Invertebrate Model": f"{data_dir}/fasta/vpod_1.2/inv_only_aligned_VPOD_1.2_het.fasta",
    "Wild-Type Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/wt_vert_aligned_VPOD_1.2_het.fasta",
    }   
    model_raw_data = {
    "Whole Dataset Model": f"{data_dir}/fasta/vpod_1.2/wds.txt",
    "Wild-Type Model": f"{data_dir}/fasta/vpod_1.2/wt.txt",
    "Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/vert.txt",
    "Invertebrate Model": f"{data_dir}/fasta/vpod_1.2/inv_only.txt",
    "Wild-Type Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/wt_vert.txt",
    }   
    model_metadata = {
    "Whole Dataset Model": f"{data_dir}/fasta/vpod_1.2/wds_meta.tsv",
    "Wild-Type Model": f"{data_dir}/fasta/vpod_1.2/wt_meta.tsv",
    "Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/vert_meta.tsv",
    "Invertebrate Model": f"{data_dir}/fasta/vpod_1.2/inv_meta.tsv",
    "Wild-Type Vertebrate Model": f"{data_dir}/fasta/vpod_1.2/wt_vert_meta.tsv",
    }   
    model_blast_db = {
    "Whole Dataset Model": f"{data_dir}/blast_dbs/vpod_1.2/wds_db",
    "Wild-Type Model": f"{data_dir}/blast_dbs/vpod_1.2/wt_db",
    "Vertebrate Model": f"{data_dir}/blast_dbs/vpod_1.2/vert_db",
    "Invertebrate Model": f"{data_dir}/blast_dbs/vpod_1.2/invert_db",
    "Wild-Type Vertebrate Model": f"{data_dir}/blast_dbs/vpod_1.2/wt_vert_db",
    }   

    model_dir = "/home/PIA/galaxy/tools/optics/models"
    #add if statements here for the different encoding methods...
    if encoding_method == 'aa_prop':
        model_directories = {
            "Whole Dataset Model": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wds_gbr.pkl",
            "Wild-Type Model": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wt_gbr.pkl",
            "Vertebrate Model": f"{model_dir}/reg_models/vpod_1.2/aa_prop/vert_gbr.pkl",
            "Invertebrate Model": f"{model_dir}/reg_models/vpod_1.2/aa_prop/invert_gbr.pkl",
            "Wild-Type Vertebrate Model": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wt_vert_gbr.pkl",
        }
        
        model_bs_dirs = {
            "Whole Dataset Model": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wds_bootstrap",
            "Wild-Type Model": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_bootstrap",
            "Vertebrate Model": f"{model_dir}/bs_models/vpod_1.2/aa_prop/vert_bootstrap",
            "Invertebrate Model": f"{model_dir}/bs_models/vpod_1.2/aa_prop/invert_bootstrap",
            "Wild-Type Vertebrate Model": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_vert_bootstrap",
        }
    else:
        model_directories = {
            "Whole Dataset Model": f"{model_dir}/reg_models/vpod_1.2/one_hot/wds_xgb.pkl",
            "Wild-Type Model": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_xgb.pkl",
            "Vertebrate Model": f"{model_dir}/reg_models/vpod_1.2/one_hot/vert_xgb.pkl",
            "Invertebrate Model": f"{model_dir}/reg_models/vpod_1.2/one_hot/invert_BayesianRidge.pkl",
            "Wild-Type Vertebrate Model": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_vert_xgb.pkl",
        }
        
        model_bs_dirs = {
            "Whole Dataset Model": f"{model_dir}/bs_models/vpod_1.2/one_hot/wds_bootstrap",
            "Wild-Type Model": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_bootstrap",
            "Vertebrate Model": f"{model_dir}/bs_models/vpod_1.2/one_hot/vert_bootstrap",
            "Invertebrate Model": f"{model_dir}/bs_models/vpod_1.2/one_hot/invert_bootstrap",
            "Wild-Type Vertebrate Model": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_vert_bootstrap",
        }


    if sequence == None:
        return ('Error: No sequence given')
    #print(sequence)
    if selected_model == None:
        return ('Error: No model selected')
    
    alignment_data = model_datasets[selected_model]
    raw_data = model_raw_data[selected_model]
    metadata = model_metadata[selected_model]
    blast_db = model_blast_db[selected_model]
    model_bs_folder = model_bs_dirs[selected_model]
    model = model_directories[selected_model]

    random_number = str(random.randint(1, 10000))
    temp_seq = f'/home/PIA/galaxy/tools/optics/tmp/temp_seq_{random_number}.fasta'
    with open(temp_seq, "w") as temp_file:  # Key change
        if '>' in sequence:
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence) # Write your data to the file object

    if blastp == 'no' or blastp == False:
        pass
    else:
        seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile)
        print('Query sequence processed via blastp')

    new_ali = f'/home/PIA/galaxy/tools/optics/tmp/temp_ali_{random_number}.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    cmd = ['mafft', '--add', temp_seq, '--keeplength', alignment_data]
    with open(new_ali, 'w') as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)

    seq_type = 'aa'
    new_seq_test = read_data(new_ali, seq_type = seq_type, is_main=True, gap_threshold=0.5)
    ref_copy = read_data(alignment_data, seq_type = seq_type, is_main=True, gap_threshold=0.5)
    last_seq = int(ref_copy.shape[0])
    new_seq_test = new_seq_test.iloc[last_seq:].copy()
    #print(new_seq_test)

    try:
        os.remove(new_ali)
        os.remove(temp_seq)
    except FileNotFoundError:
        raise Exception("File does not exist")

    if bootstrap == 'yes':
        # ... (Load the selected model and make a bootstrap prediction)
        mean_prediction, ci_lower, ci_upper, prediction_dict = calculate_ensemble_CI(model_bs_folder, new_seq_test, name, prediction_dict)
        print(f'{mean_prediction}, {ci_lower}, {ci_upper}, {prediction_dict}')
        
        # ... (Load the selected model and make a single prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return (round(float(mean_prediction),1), round(float(ci_lower),1), round(float(ci_upper),1), prediction_dict, round(float(prediction[0]),1))

    else:
        # ... (Load the selected model and make a prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return(round(float(prediction[0]),1))
 

def process_sequences_from_file(file,selected_model, identity_report, blastp, refseq, reffile, bootstrap, encoding_method):
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
                    names.append(line.replace('>','').strip().replace(' ','_'))
                    sequences.append(entry)
                    entry = ""
                    entry += line
                    #print(sequences)
                    line_count+=1
                else:
                    names.append(line.replace('>','').strip().replace(' ','_'))
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
    mean_predictions = []
    ci_lowers = []
    ci_uppers = []
    prediction_dict = {}
    i = 0
    for seq in sequences:
        if bootstrap == 'no':    
            print(seq)
            prediction = process_sequence(seq, names[i], selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method)  # Process each sequence
            predictions.append(prediction)

        else:
            if len(names) > 10:
                raise(Exception("More than TEN sequences detected for bootstrap prediction protocol.\nDue to computing limitations you will need to resubmit with only 10 sequences or less."))
            mean_prediction, ci_lower, ci_upper, prediction_dict, prediction = process_sequence(seq, names[i], selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method)  # Process each sequence
            mean_predictions.append(mean_prediction)
            predictions.append(prediction)
            ci_lowers.append(ci_lower)
            ci_uppers.append(ci_upper)

        i+=1
    #print(predictions)

    return(names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions)

def main():
    
    model_dir = "/home/PIA/galaxy/tools/optics/models/reg_models/vpod_1.2"
    model_directories = {
        "Whole Dataset Model": f"{model_dir}/wds_gbc.pkl",
        "Wild-Type Model": f"{model_dir}/wt_gbc.pkl",
        "Vertebrate Model": f"{model_dir}/vert_gbc.pkl",
        "Invertebrate Model": f"{model_dir}/invert_lgbm.pkl",
        "Wild-Type Vertebrate Model": f"{model_dir}/wt_vert_gbc.pkl",
    }
    encoding_methods=['one_hot','aa_prop']

    parser = argparse.ArgumentParser(description="Process sequences using a selected model")
    parser.add_argument("input", help="Either a single sequence or a path to a FASTA file", type=str)
    parser.add_argument("output", help="Name for output file", type=str, default = 'output.txt')
    parser.add_argument("iden_output", help="Name for the sequence identity report output file", type=str, default = 'seq_identity_report.txt')
    parser.add_argument("pdffile", help="Name for the pdf file output file for visualizing bootstrap predictions", type=str, default = 'pdf_dummy.txt')

    parser.add_argument("-m", "--model", help="Model to use for prediction", 
                    choices=list(model_directories.keys()), required=True)
    parser.add_argument("-b", "--blastp", help="Option to enable blastp analsis on query sequences", 
                    type = str or bool , required=True)
    parser.add_argument("-r", "--refseq", help="Reference sequence used for blastp analysis.", 
                type = str, required=False)
    parser.add_argument("-f", "--reffile", help="Cutom reference sequence file used for blastp analysis.", 
                type = str, required=False)
    parser.add_argument("-s", "--bootstrap", help="Option to enable bootstrap predictions on query sequences", 
                    type = str or bool , required=True)
    parser.add_argument("-e", "--encoding_method", help="Select preferred encoding method used to train model and make predictions", 
                choices=encoding_methods, type = str, required=True)

    args = parser.parse_args()

    if args.blastp == "no":
        with open(args.iden_output, 'w') as f:
            f.write('Blastp analysis did not occur.')

    if args.bootstrap == "no":
        pdf = FPDF()
        pdf.add_page()
        pdf.set_font('Arial', size=12)
        pdf.cell(40, 10, txt="Boostraping was not enabled, so this file is empty.")
        pdf.output(args.pdffile)

    if os.path.isfile(args.input):
        names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions  = process_sequences_from_file(args.input, args.model, args.iden_output, args.blastp, args.refseq, args.reffile, args.bootstrap, args.encoding_method)
        with open(args.output, 'w') as f:
            i = 0
            while i in range(len(names)):
                if args.bootstrap == 'no':
                    if i == 0:
                        f.write('Names\tPredictions\n')
                    f.write(f"{names[i]}:\t{predictions[i]}\n")
                    print(f"{names[i]}:\t{predictions[i]}\n")
                    i+=1
                else:
                    if i == 0:
                        f.write('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\n')
                        median_predictions = plot_predictions_with_CI(names, prediction_dict, mean_predictions, args.pdffile)
                    f.write(f"{names[i]}:\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\n")
                    print(f"{names[i]}:\t{predictions[i]}\t{mean_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\n")
                    i+=1
            f.write(f"\nModel Used:\t{args.model}\nEncoding Method:\t{args.encoding_method}\n")
    else:
        raise Exception("No file passed to input for predictions")

if __name__ == "__main__":
    main()