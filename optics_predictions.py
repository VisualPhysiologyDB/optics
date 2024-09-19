# prediction_functions.py
import subprocess 
from optics_scripts.deepBreaks.utils import load_obj
from optics_scripts.deepBreaks.preprocessing import read_data
from fpdf import FPDF
import pandas as pd
import argparse
import os
import random
from optics_scripts.blastp_align import seq_sim_report
from optics_scripts.bootstrap_predictions import calculate_ensemble_CI , plot_predictions_with_CI 

def process_sequence(sequence, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict = None, encoding_method='one_hot'):
    data_dir = "./data"
    model_datasets = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds_aligned_VPOD_1.2_het.fasta",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt_aligned_VPOD_1.2_het.fasta",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert_aligned_VPOD_1.2_het.fasta",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_only_aligned_VPOD_1.2_het.fasta",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert_aligned_VPOD_1.2_het.fasta",
    }   
    model_raw_data = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds.txt",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt.txt",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert.txt",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_only.txt",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert.txt",
    }   
    model_metadata = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds_meta.tsv",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt_meta.tsv",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert_meta.tsv",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_meta.tsv",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert_meta.tsv",
    }   
    model_blast_db = {
    "whole-dataset": f"{data_dir}/blast_dbs/vpod_1.2/wds_db",
    "wildtype": f"{data_dir}/blast_dbs/vpod_1.2/wt_db",
    "vertebrate": f"{data_dir}/blast_dbs/vpod_1.2/vert_db",
    "invertebrate": f"{data_dir}/blast_dbs/vpod_1.2/invert_db",
    "wildtype-vert": f"{data_dir}/blast_dbs/vpod_1.2/wt_vert_db",
    }   

    model_dir = "./models"
    #add if statements here for the different encoding methods...
    if encoding_method == 'aa_prop':
        model_directories = {
            "whole-dataset": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wds_gbr.pkl",
            "wildtype": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wt_gbr.pkl",
            "vertebrate": f"{model_dir}/reg_models/vpod_1.2/aa_prop/vert_gbr.pkl",
            "invertebrate": f"{model_dir}/reg_models/vpod_1.2/aa_prop/invert_gbr.pkl",
            "wildtype-vert": f"{model_dir}/reg_models/vpod_1.2/aa_prop/wt_vert_gbr.pkl",
        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wds_bootstrap",
            "wildtype": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_bootstrap",
            "vertebrate": f"{model_dir}/bs_models/vpod_1.2/aa_prop/vert_bootstrap",
            "invertebrate": f"{model_dir}/bs_models/vpod_1.2/aa_prop/invert_bootstrap",
            "wildtype-vert": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_vert_bootstrap",
        }
    else:
        model_directories = {
            "whole-dataset": f"{model_dir}/reg_models/vpod_1.2/one_hot/wds_xgb.pkl",
            "wildtype": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_xgb.pkl",
            "vertebrate": f"{model_dir}/reg_models/vpod_1.2/one_hot/vert_xgb.pkl",
            "invertebrate": f"{model_dir}/reg_models/vpod_1.2/one_hot/invert_BayesianRidge.pkl",
            "wildtype-vert": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_vert_xgb.pkl",
        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/vpod_1.2/one_hot/wds_bootstrap",
            "wildtype": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_bootstrap",
            "vertebrate": f"{model_dir}/bs_models/vpod_1.2/one_hot/vert_bootstrap",
            "invertebrate": f"{model_dir}/bs_models/vpod_1.2/one_hot/invert_bootstrap",
            "wildtype-vert": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_vert_bootstrap",
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
    temp_seq = f'./tmp/temp_seq_{random_number}.fasta'
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

    new_ali = f'./tmp/temp_ali_{random_number}.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    
    try:
        print('Trying Linux execution of MAFFT')
        cmd = ['mafft', '--add', temp_seq, '--keeplength', alignment_data]
        with open(new_ali, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
    except:
        print('Trying Windows execution of MAFFT')
        mafft_exe = './optics_scripts/mafft/mafft-win/mafft.bat'
        cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
        with open(new_ali, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
    finally: 
        print('Trying Mac execution of MAFFT')
        mafft_exe = './optics_scripts/mafft/mafft-mac/mafft.bat'
        cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
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
        raise Exception('Error: No file given')
    
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
        if bootstrap == 'no' or bootstrap == False:    
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
    

    models = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert']
    encoding_methods=['one_hot','aa_prop']
    ref_seq_choices = ['bovine', 'squid', 'custom']

    parser = argparse.ArgumentParser(description="Process sequences using a selected model")
    parser.add_argument("-in","--input", help="Either a single sequence or a path to a FASTA file", type=str, required = True)
    parser.add_argument("-out","--output", help="Name for output file", type=str, default = 'optics_predictions.txt', required = False)
    parser.add_argument("-ir","--iden_output", help="Name for the sequence identity report output file", type=str, default = 'blastp_report.txt', required = False)
    parser.add_argument("-bsv","--pdffile", help="Name for the pdf file output file for visualizing bootstrap predictions", type=str, default = 'bootstrap_viz.pdf', required = False)
    parser.add_argument("-m", "--model", help="Model to use for prediction", 
                    choices=models, default="whole-dataset", required=False)
    parser.add_argument("-b", "--blastp", help="Option to enable blastp analsis on query sequences", 
                    type = str or bool , default = True, required=False)
    parser.add_argument("-r", "--refseq", help="Reference sequence used for blastp analysis.", 
                type = str, choisces = ref_seq_choices, default= 'bovine', required=False)
    parser.add_argument("-f", "--reffile", help="Custom reference sequence file used for blastp analysis.", 
                type = str, default = 'not_real.txt', required=False)
    parser.add_argument("-s", "--bootstrap", help="Option to enable bootstrap predictions on query sequences", 
                    type = str or bool , deafult = False, required=False)
    parser.add_argument("-e", "--encoding_method", help="Select preferred encoding method used to train model and make predictions", 
                choices=encoding_methods, default='aa_prop', type = str, required=False)
    # python prediction_functions_galaxy.py -in -out -ir -bsv -m -b -r -f -s -e

    args = parser.parse_args()

    if args.blastp == "no" or args.blastp == False:
        with open(args.iden_output, 'w') as f:
            f.write('Blastp analysis did not occur.')

    if args.bootstrap == "no" or args.bootstrap == False:
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
                if args.bootstrap == 'no' or args.bootstrap == False:
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