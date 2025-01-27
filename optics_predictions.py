# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import pandas as pd
import json
import argparse
import os
import sys
import numpy as np
import random
import copy
import datetime
import matplotlib
from joblib import Parallel, delayed
import tempfile
from multiprocessing import Manager
from tqdm import tqdm
from optics_scripts.blastp_align import seq_sim_report
from optics_scripts.bootstrap_predictions import calculate_ensemble_CI, plot_prediction_subsets_with_CI, wavelength_to_rgb

def extract_fasta(file):
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
    return names,sequences

def process_sequence(sequence=None, name=None, selected_model=None, identity_report=None, blastp=None, refseq=None, reffile=None, bootstrap=None, prediction_dict=None, encoding_method='one_hot', wrk_dir = '', only_blast = False):
    data_dir = f"./data"
    model_datasets = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds_aligned_VPOD_1.2_het.fasta",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt_aligned_VPOD_1.2_het.fasta",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert_aligned_VPOD_1.2_het.fasta",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_only_aligned_VPOD_1.2_het.fasta",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert_aligned_VPOD_1.2_het.fasta",
    "type-one": f"{data_dir}/fasta/vpod_1.2/Karyasuyama_T1_ops_aligned.fasta"
    }   
    model_raw_data = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds.txt",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt.txt",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert.txt",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_only.txt",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert.txt",
    "type-one": f"{data_dir}/fasta/vpod_1.2/Karyasuyama_T1_ops.txt"

    }   
    model_metadata = {
    "whole-dataset": f"{data_dir}/fasta/vpod_1.2/wds_meta.tsv",
    "wildtype": f"{data_dir}/fasta/vpod_1.2/wt_meta.tsv",
    "vertebrate": f"{data_dir}/fasta/vpod_1.2/vert_meta.tsv",
    "invertebrate": f"{data_dir}/fasta/vpod_1.2/inv_meta.tsv",
    "wildtype-vert": f"{data_dir}/fasta/vpod_1.2/wt_vert_meta.tsv",
    "type-one": f"{data_dir}/fasta/vpod_1.2/Karyasuyama_T1_ops_meta.tsv"
    }   
    model_blast_db = {
    "whole-dataset": f"{data_dir}/blast_dbs/vpod_1.2/wds_db",
    "wildtype": f"{data_dir}/blast_dbs/vpod_1.2/wt_db",
    "vertebrate": f"{data_dir}/blast_dbs/vpod_1.2/vert_db",
    "invertebrate": f"{data_dir}/blast_dbs/vpod_1.2/invert_db",
    "wildtype-vert": f"{data_dir}/blast_dbs/vpod_1.2/wt_vert_db",
    "type-one": f"{data_dir}/blast_dbs/vpod_1.2/t1_db",
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
            "type-one": f"{model_dir}/reg_models/vpod_1.2/aa_prop/t1_xgb.pkl",
        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wds_bootstrap",
            "wildtype": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_bootstrap",
            "vertebrate": f"{model_dir}/bs_models/vpod_1.2/aa_prop/vert_bootstrap",
            "invertebrate": f"{model_dir}/bs_models/vpod_1.2/aa_prop/invert_bootstrap",
            "wildtype-vert": f"{model_dir}/bs_models/vpod_1.2/aa_prop/wt_vert_bootstrap",
            "type-one": f"{model_dir}/bs_models/vpod_1.2/aa_prop/t1_bootstrap",
        }
        

    else:
        model_directories = {
            "whole-dataset": f"{model_dir}/reg_models/vpod_1.2/one_hot/wds_xgb.pkl",
            "wildtype": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_xgb.pkl",
            "vertebrate": f"{model_dir}/reg_models/vpod_1.2/one_hot/vert_xgb.pkl",
            "invertebrate": f"{model_dir}/reg_models/vpod_1.2/one_hot/invert_BayesianRidge.pkl",
            "wildtype-vert": f"{model_dir}/reg_models/vpod_1.2/one_hot/wt_vert_xgb.pkl",
            "type-one": f"{model_dir}/reg_models/vpod_1.2/one_hot/t1_xgb.pkl",

        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/vpod_1.2/one_hot/wds_bootstrap",
            "wildtype": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_bootstrap",
            "vertebrate": f"{model_dir}/bs_models/vpod_1.2/one_hot/vert_bootstrap",
            "invertebrate": f"{model_dir}/bs_models/vpod_1.2/one_hot/invert_bootstrap",
            "wildtype-vert": f"{model_dir}/bs_models/vpod_1.2/one_hot/wt_vert_bootstrap",
            "type-one": f"{model_dir}/bs_models/vpod_1.2/one_hot/t1_bootstrap",
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
    
    #wrk_dir = os.getcwd().replace('\\','/')
    with tempfile.NamedTemporaryFile(mode="w", dir=f"{wrk_dir}/tmp", suffix=".fasta", delete=False) as temp_seq_file:
        temp_seq = temp_seq_file.name  # Get the unique filename
        if '>' in sequence:
            temp_seq_file.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            temp_seq_file.write(sequence)
    
    # This is a special case for when a cached sequence is detected
    if only_blast == True:
        if blastp == 'no' or blastp == False or blastp == 'False':
            percent_iden = '-'
        else:
            percent_iden = seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile)
        return percent_iden

    if blastp == 'no' or blastp == False or blastp == 'False':
        percent_iden = '-'
    else:
        percent_iden = seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile)
        #print('Query sequence processed via blastp')

    with tempfile.NamedTemporaryFile(mode="w", dir=f"{wrk_dir}/tmp", suffix=".fasta", delete=False) as temp_ali_file:
        new_ali = temp_ali_file.name 
    try:
        #print('Trying Linux execution of MAFFT')
        cmd = ['mafft', '--add', temp_seq, '--keeplength', alignment_data]
        with open(new_ali, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
    except:
        #print('Trying Windows execution of MAFFT')
        try:
            mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-win/mafft.bat'
            cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
            with open(new_ali, 'w') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
        except: 
            #print('Trying Mac execution of MAFFT')
            try:
                mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-mac/mafft.bat'
                cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
                with open(new_ali, 'w') as f:
                    subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                raise Exception(f'MAFFT alignment failed for all three options (Linux, Windows, and Max).\nCheck your FASTA file to make sure it is formatted correctly as that could be a source of error.\n{e.stderr.decode()}')
            
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

    if bootstrap == True:
        # ... (Load the selected model and make a bootstrap prediction)
        mean_prediction, ci_lower, ci_upper, prediction_dict, median_prediction, std_dev, predictions_all = calculate_ensemble_CI(model_bs_folder, new_seq_test, name, prediction_dict)
        #print(f'{mean_prediction}, {ci_lower}, {ci_upper}, {prediction_dict}')
        
        # ... (Load the selected model and make a single prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return (round(float(mean_prediction),1), round(float(ci_lower),1), round(float(ci_upper),1), prediction_dict, round(float(prediction[0]),1), round(float(median_prediction),1), str(percent_iden), round(float(std_dev),1), predictions_all.tolist())

    else:
        # ... (Load the selected model and make a prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return(round(float(prediction[0]),1), str(percent_iden))
 

def process_sequences_from_file(file, selected_model, identity_report, blastp, refseq, reffile, bootstrap, encoding_method):
    wrk_dir = os.getcwd().replace('\\','/')
    cache_dir = f"./data/cached_predictions/"
    
    if (bootstrap == True or bootstrap == 'True' or bootstrap == 'true' or bootstrap == 'yes'):
        bootstrap = True
        model_type = 'bs_models'
    else:
        bootstrap = False
        model_type = 'reg_models'
    model_cache_dirs = {
        "whole-dataset": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/wds_pred_dict.json",
        "wildtype": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/wt_pred_dict.json",
        "vertebrate": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/vert_pred_dict.json",
        "invertebrate": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/invert_pred_dict.json",
        "wildtype-vert": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/wt_pred_dict.json",
        "type-one": f"{cache_dir}/{model_type}/vpod_1.2/{encoding_method}/t1_pred_dict.json",
    }

    cache_file = model_cache_dirs[selected_model]
    
    if file == None:
        raise Exception('Error: No file given')
    
    names,sequences = extract_fasta(file)       
    predictions = []
    mean_predictions = []
    median_predictions = []
    std_dev_list = []
    ci_lowers = []
    ci_uppers = []
    per_iden_list = []
    seq_lens = []
    # Load the existing prediction dictionary to pull from if a sequence has already been predicted by the chosen model
    # Check to see if an existing prediiction dictionary already exists to save time.
    # In this dictionary the sequence will be our keys, and the len(seq), prediction, percent_iden, mean_prediction, ci_lower, ci_upper, median_prediction, std_dev will be the corresponding values...
    if os.path.isfile(cache_file):
        try:
            with open(cache_file, 'r') as f:
                cached_pred_dict = json.load(f)
        except FileNotFoundError:
            cached_pred_dict = {}
    else:
        cached_pred_dict = {}
    
    manager = Manager()
    prediction_dict = manager.dict()  # Use a shared dictionary for prediction outputs
    mp_cached_pred_dict = manager.dict(cached_pred_dict)  # Initialize a multiprocess available cached dict with cached data
    
    def process_sequence_wrapper(seq, name):  # Helper function
        just_seq = seq.split('\n')[1]
        if bootstrap == False:
            if just_seq not in list(mp_cached_pred_dict.keys()):
                prediction, percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method, wrk_dir)
                mp_cached_pred_dict[just_seq] = {'len': len(just_seq), 'single_prediction': prediction}
                return len(just_seq), prediction, percent_iden, None, None, None, None, None  # Consistent return values
            else:
                percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method, wrk_dir, only_blast=True)
                return mp_cached_pred_dict[just_seq]['len'], mp_cached_pred_dict[just_seq]['single_prediction'], percent_iden, None, None, None, None, None  # Consistent return values

        else:
            if just_seq not in list(mp_cached_pred_dict.keys()):
                mean_prediction, ci_lower, ci_upper, updated_prediction_dict, prediction, median_prediction, percent_iden, std_dev, predictions_all = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method, wrk_dir)
                prediction_dict.update(updated_prediction_dict) # Update the shared dictionary with the returned dictionary
                mp_cached_pred_dict[just_seq] = {'len': len(just_seq),
                                            'single_prediction': prediction,
                                            'mean_prediction' : mean_prediction,
                                            'ci_lower' : ci_lower,
                                            'ci_upper' : ci_upper,
                                            'median_prediction' : median_prediction,
                                            'std_deviation' : std_dev,
                                            'all_bs_predictions' : predictions_all
                                            }
                # If the seq is in our cache, I need to be able to update the prediction_dict with that info...
                return len(just_seq), prediction, percent_iden, mean_prediction, ci_lower, ci_upper, median_prediction, std_dev
            else:
                percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method, wrk_dir, only_blast=True)
                prediction_dict[name] = np.array(mp_cached_pred_dict[just_seq]['all_bs_predictions'])
                return mp_cached_pred_dict[just_seq]['len'], mp_cached_pred_dict[just_seq]['single_prediction'], percent_iden, mp_cached_pred_dict[just_seq]['mean_prediction'], mp_cached_pred_dict[just_seq]['ci_lower'], mp_cached_pred_dict[just_seq]['ci_upper'], mp_cached_pred_dict[just_seq]['median_prediction'], mp_cached_pred_dict[just_seq]['std_deviation']


    try:
        with tqdm_joblib(tqdm(total=len(sequences), desc="Processing Sequences", ascii = True, bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]", dynamic_ncols=True)) as pbar:  # Use tqdm for progress bar
            results = Parallel(n_jobs=-1)(delayed(process_sequence_wrapper)(seq, names[i]) for i, seq in enumerate(sequences))
    except Exception as e:
        updated_cached_pred_dict = dict(mp_cached_pred_dict)        
        # Save the taxon dictionary if it doesn't yet exist or if it has been updated since being loaded 
        if list(updated_cached_pred_dict.keys()) > list(cached_pred_dict.keys()):  
            #print('Saving Updated Dictionary') 
            try:
                with open(cache_file, 'w') as f:
                    json.dump(updated_cached_pred_dict, f, indent=4)  # indent for pretty formatting
                    raise Exception(f"Error Occured During the Prediction Process:\n{e}\n")

            except FileNotFoundError:
                raise Exception(f"Error: Cached prediction file can't be saved...\n")

    # Extract results
    for result in results:
        seq_len, prediction, percent_iden, mean_prediction, ci_lower, ci_upper, median_prediction, std_dev = result
        seq_lens.append(seq_len)
        predictions.append(prediction)
        per_iden_list.append(percent_iden)
        if mean_prediction is not None:  # Handle bootstrap case
            mean_predictions.append(mean_prediction)
            ci_lowers.append(ci_lower)
            ci_uppers.append(ci_upper)
            median_predictions.append(median_prediction)
            std_dev_list.append(std_dev)
            
    updated_cached_pred_dict = dict(mp_cached_pred_dict)        
    # Save the taxon dictionary if it doesn't yet exist or if it has been updated since being loaded 
    if list(updated_cached_pred_dict.keys()) > list(cached_pred_dict.keys()):  
        #print('Saving Updated Dictionary') 
        try:
            with open(cache_file, 'w') as f:
                json.dump(updated_cached_pred_dict, f, indent=4)  # indent for pretty formatting
        except FileNotFoundError:
            print(f"Error: Cached prediction file can't be saved...\n")

    return(names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens)

def main():
    
    models = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert','type-one']
    encoding_methods=['one_hot','aa_prop']
    ref_seq_choices = ['bovine', 'squid', 'microbe','custom']
    bool_choices = ['true', 'True', 'false','False', True, False]
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    
    parser = argparse.ArgumentParser(description="Process sequences using a selected model")
    parser.add_argument("-in","--input", help="Either a single sequence or a path to a FASTA file", type=str, required = True)
    parser.add_argument("-rd","--report_dir", help="Name of folder directory to create", type=str, required = False, default = f'optics_on_unamed_{dt_label}')
    parser.add_argument("-out","--output", help="Name for output file", type=str, default = 'optics_predictions.txt', required = False)
    parser.add_argument("-m", "--model", help="Model to use for prediction", 
                    choices=models, default="whole-dataset", required=False)
    parser.add_argument("-e", "--encoding_method", help="Select preferred encoding method used to train model and make predictions", 
                    choices=encoding_methods, default='aa_prop', type = str, required=False)
    parser.add_argument("-b", "--blastp", help="Option to enable blastp analsis on query sequences", 
                    type = str or bool, choices=bool_choices, default = True, required=False)
    parser.add_argument("-ir","--iden_report", help="Name for the blastp report output file", type=str, default = 'blastp_report.txt', required = False)
    parser.add_argument("-r", "--refseq", help="Reference sequence used for blastp analysis.", 
                type = str, choices = ref_seq_choices, default= 'bovine', required=False)
    parser.add_argument("-f", "--reffile", help="Custom reference sequence file used for blastp analysis.", 
                type = str, default = '', required=False)
    parser.add_argument("-s", "--bootstrap", help="Option to enable bootstrap predictions on query sequences", 
                    type = str or bool , choices=bool_choices, default = True, required=False)
    parser.add_argument("-viz", "--visualize_bootstrap", help="Option to enable/disable visualization of bootstrap predictions on query sequences", 
                type = str or bool, choices=bool_choices, default = True, required=False)
    parser.add_argument("-bsv","--bootstrap_viz_file", help="Name for the pdf file output file for visualizing bootstrap predictions", type=str, default = 'bootstrap_viz', required = False)

    # python optics_predictions.py -in ./examples/msp_erg_raw.txt -rd msp_test_of_optics -out msp_predictions.tsv -m whole-dataset -e aa_prop -b True -ir msp_blastp_report.tsv -r squid -s False
    # python optics_predictions.py -in ./examples/msp_erg_raw.txt -rd msp_test_of_optics -out msp_predictions.tsv -m whole-dataset -e aa_prop -b True -ir msp_blastp_report.tsv -r squid -s True -bsv msp_bs_viz

    args = parser.parse_args()

    if os.path.isdir('./tmp'):
        pass
    else:    
        os.makedirs(f'./tmp')
        
    if os.path.isdir(f'./prediction_outputs'):
        pass
    else:    
        os.makedirs(f'./prediction_outputs')   

    if 'optics_on_unamed' in args.report_dir:
        report_dir = args.report_dir
    else:
        report_dir = f'./prediction_outputs/optics_on_{args.report_dir}_{dt_label}'
    os.makedirs(report_dir)

    if '.txt' in args.iden_report or 'tsv in args.iden_report':
        blastp_file = f'{report_dir}/{args.iden_report}'
    else:
        blastp_file = f'{report_dir}/{args.iden_report}.txt'

    bootstrap_file = f'{report_dir}/{args.bootstrap_viz_file}'
   
    log_file = f'{report_dir}/arg_log.txt'
    
    if os.path.isfile(args.input):
        names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens_list = process_sequences_from_file(args.input, args.model, blastp_file, args.blastp, args.refseq, args.reffile, args.bootstrap, args.encoding_method)
        if '.tsv' in args.output or '.txt' in args.output:
            output = f'{report_dir}/{args.output}'
            sub_output = args.output.replace('.tsv','')
            excel_output = f'{report_dir}/{sub_output}_for_excel.xlsx'
        else:
            output = f'{report_dir}/{args.output}.tsv'
            excel_output = f'{report_dir}/{args.output}_for_excel.xlsx'

        with open(output, 'w') as f:
            i = 0
            while i in range(len(names)):
                if args.bootstrap == False or args.bootstrap == 'False' or args.bootstrap == 'false':
                    if i == 0:
                        f.write('Names\tPredictions\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\n')
                        colors = [wavelength_to_rgb(pred) for pred in predictions]
                        hex_color_list = [matplotlib.colors.to_hex(color) for color in colors]
                        write_to_excel(names, predictions, per_iden_list, excel_output, hex_color_list=hex_color_list, seq_lens_list=seq_lens_list)
                    f.write(f"{names[i]}\t{predictions[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")
                    print(f"{names[i]}\t{predictions[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")
                    i+=1

                else:
                    if i == 0:
                        f.write('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\tStd_Deviation\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\tLmax_Hex_Color\n')
                        print('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\tStd_Deviation\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\n')
                        
                        # colors for hex_color_list generated from the mean prediction of the bootstraped predictions during the visulization steps    
                        hex_color_list = plot_prediction_subsets_with_CI(names, prediction_dict, mean_predictions, 
                                                                         bootstrap_file, args.visualize_bootstrap)
                        
                        write_to_excel(names, predictions, per_iden_list, excel_output, 
                                        mean_predictions, median_predictions, ci_lowers, 
                                        ci_uppers, std_dev_list, hex_color_list, seq_lens_list)
                        
                    f.write(f"{names[i]}\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\t{std_dev_list[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\t{hex_color_list[i]}\n")
                    print(f"{names[i]}\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\t{std_dev_list[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")
                    i+=1
        with open(log_file, 'w') as f:
            command_line_input = ' '.join(sys.argv)
            f.write(f"Command executed:\t{command_line_input}\n")
            f.write(f"Model Used:\t{args.model}\nEncoding Method:\t{args.encoding_method}\n")
            print(f"\nModel Used:\t{args.model}\nEncoding Method:\t{args.encoding_method}\n")
                            
                        
        with open(f'{report_dir}/fig_tree_color_annotation.txt', 'w') as g:
            g.write("Name\t!color\n")  # Header row
            for name, hex_color in zip(names, hex_color_list):
                g.write(f"{name}\t{hex_color}\n") 
        with open(f'{report_dir}/itol_color_annotation.txt', 'w') as g:
            g.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
            for name, hex_color in zip(names, hex_color_list):
                g.write(f"{name}\tlabel_background\t{hex_color}\n") 
            
        print('Predictions Complete!')
        #os.remove('./tmp')

    else:
        raise Exception("No file passed to input for predictions")
    
from openpyxl import Workbook
from openpyxl.styles import PatternFill

def write_to_excel(names, predictions, per_iden_list, output_filename="output.xlsx", 
                  mean_predictions=None, median_predictions=None, ci_lowers=None, 
                  ci_uppers=None, std_dev_list=None, hex_color_list=None, seq_lens_list=None):
    """
    Writes data to an Excel sheet, including bootstrap statistics and
    hexadecimal color codes, and colors the cells based on the hex codes.

    Args:
        names: List of names.
        predictions: List of predictions.
        per_iden_list: List of percentage identities.
        output_filename: Name of the output Excel file.
        mean_predictions: List of mean predictions (optional, for bootstrap).
        median_predictions: List of median predictions (optional, for bootstrap).
        ci_lowers: List of lower confidence intervals (optional, for bootstrap).
        ci_uppers: List of upper confidence intervals (optional, for bootstrap).
        std_dev_list: List of standard deviations (optional, for bootstrap).
        hex_color_list: List of hexadecimal color codes.
        seq_lens_list: List of sequence lengths
        """

    wb = Workbook()
    ws = wb.active

    
    if mean_predictions == None:
        ws.append(['Names', 'Single_Prediction', '%Identity_Nearest_VPOD_Sequence', 'Sequence_Length','Lmax_Hex_Color'])
        for i in range(len(names)):
            # Because openpyxel is picky about hex-codes we need to remove the '#' symbol for it to accept it as a fill color.
            hex_color = hex_color_list[i].replace('#','') 
            ws.append([names[i], predictions[i], per_iden_list[i], seq_lens_list[i], hex_color_list[i]])
            ws.cell(row=i+2, column=5).fill = PatternFill(start_color=hex_color, 
                                                        end_color=hex_color, 
                                                        fill_type="solid")
    else:
        ws.append(['Names', 'Single_Prediction', 'Prediction_Means', 'Prediction_Medians',
                    'Prediction_Lower_Bounds', 'Prediction_Upper_Bounds', 'Std_Deviation', 
                    '%Identity_Nearest_VPOD_Sequence', 'Sequence_Lengths','Lmax_Hex_Color'])
        for i in range(len(names)):
            # Because openpyxel is picky about hex-codes we need to remove the '#' symbol for it to accept it as a fill color.
            hex_color = hex_color_list[i].replace('#','') 
            ws.append([names[i], predictions[i], mean_predictions[i], median_predictions[i],
                        ci_lowers[i], ci_uppers[i], std_dev_list[i], per_iden_list[i], seq_lens_list[i], hex_color_list[i]])
            ws.cell(row=i+2, column=10).fill = PatternFill(start_color=hex_color, 
                                                        end_color=hex_color, 
                                                        fill_type="solid")
    wb.save(output_filename)
    
    
import contextlib  
import joblib  
@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()
        
        
if __name__ == "__main__":
    main()