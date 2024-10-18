# prediction_functions.py
import subprocess 
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import pandas as pd
import argparse
import os
import sys
import random
import datetime
from progress.bar import ShadyBar
from optics_scripts.blastp_align import seq_sim_report
from optics_scripts.bootstrap_predictions import calculate_ensemble_CI, plot_prediction_subsets_with_CI

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
    
    wrk_dir = os.getcwd().replace('\\','/')
    random_number = str(random.randint(1, 10000))
    temp_seq = f'{wrk_dir}/tmp/temp_seq_{random_number}.fasta'
    with open(temp_seq, "w") as temp_file:  # Key change
        if '>' in sequence:
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence)
        else:
            sequence = ">placeholder_name\n" + sequence
            #print(f"here is the sequence: {sequence}")
            temp_file.write(sequence) # Write your data to the file object

    if blastp == 'no' or blastp == False or blastp == 'False':
        percent_iden = '-'
    else:
        percent_iden = seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile)
        #print('Query sequence processed via blastp')

    new_ali = f'{wrk_dir}/tmp/temp_ali_{random_number}.fasta'  
    # ... (Perform alignment using MAFFT with alignment_data)
    
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

    if bootstrap == 'yes' or bootstrap == True or bootstrap == 'True'  or bootstrap == 'true':
        # ... (Load the selected model and make a bootstrap prediction)
        mean_prediction, ci_lower, ci_upper, prediction_dict, median_prediction, std_dev = calculate_ensemble_CI(model_bs_folder, new_seq_test, name, prediction_dict)
        #print(f'{mean_prediction}, {ci_lower}, {ci_upper}, {prediction_dict}')
        
        # ... (Load the selected model and make a single prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return (round(float(mean_prediction),1), round(float(ci_lower),1), round(float(ci_upper),1), prediction_dict, round(float(prediction[0]),1), round(float(median_prediction),1), str(percent_iden), round(float(std_dev),1))

    else:
        # ... (Load the selected model and make a prediction)
        loaded_mod = load_obj(model)
        prediction = loaded_mod.predict(new_seq_test)
        
        return(round(float(prediction[0]),1), str(percent_iden))
 

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
    median_predictions = []
    std_dev_list = []
    ci_lowers = []
    ci_uppers = []
    prediction_dict = {}
    per_iden_list = []
    seq_lens = []
    i = 0

    bar = ShadyBar('Processing Sequences', max=len(names), charset='ascii')
    for seq in sequences:
        seq_lens.append(len(seq))
        if bootstrap == 'no' or bootstrap == False or bootstrap == 'False' or bootstrap == 'false':    
            #print(seq)
            prediction, percent_iden = process_sequence(seq, names[i], selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method)  # Process each sequence
            predictions.append(prediction)
            per_iden_list.append(percent_iden)
        else:
            mean_prediction, ci_lower, ci_upper, prediction_dict, prediction, median_prediction, percent_iden, std_dev = process_sequence(seq, names[i], selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, encoding_method)  # Process each sequence
            mean_predictions.append(mean_prediction)
            predictions.append(prediction)
            ci_lowers.append(ci_lower)
            ci_uppers.append(ci_upper)
            median_predictions.append(median_prediction)
            per_iden_list.append(percent_iden)
            std_dev_list.append(std_dev)
        bar.next()
        i+=1
    #print(predictions)
    bar.finish()
    return(names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens)

def main():
    
    models = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert']
    encoding_methods=['one_hot','aa_prop']
    ref_seq_choices = ['bovine', 'squid', 'custom']
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
                    type = str or bool , default = True, required=False)
    parser.add_argument("-ir","--iden_report", help="Name for the blastp report output file", type=str, default = 'blastp_report.txt', required = False)
    parser.add_argument("-r", "--refseq", help="Reference sequence used for blastp analysis.", 
                type = str, choices = ref_seq_choices, default= 'bovine', required=False)
    parser.add_argument("-f", "--reffile", help="Custom reference sequence file used for blastp analysis.", 
                type = str, default = 'not_real.txt', required=False)
    parser.add_argument("-s", "--bootstrap", help="Option to enable bootstrap predictions on query sequences", 
                    type = str or bool , default = True, required=False)
    parser.add_argument("-viz", "--visualize_bootstrap", help="Option to enable/disable visualization of bootstrap predictions on query sequences", 
                type = str or bool , default = True, required=False)
    parser.add_argument("-bsv","--bootstrap_viz_file", help="Name for the pdf file output file for visualizing bootstrap predictions", type=str, default = 'bootstrap_viz', required = False)

    # python optics_predictions.py -in ./examples/msp_erg_raw.txt -rd msp_test_of_optics -out msp_predictions.tsv -m whole-dataset -e aa_prop -b True -ir msp_blastp_report.tsv -r squid -s False
    # python optics_predictions.py -in ./examples/msp_erg_raw.txt -rd msp_test_of_optics -out msp_predictions.tsv -m whole-dataset -e aa_prop -b True -ir msp_blastp_report.tsv -r squid -s True -bsv msp_bs_viz

    args = parser.parse_args()
    
    if os.path.isdir('./tmp'):
        pass
    else:    
        os.makedirs(f'./tmp')
    
    if os.path.isdir('./prediction_outputs'):
        pass
    else:    
        os.makedirs('./prediction_outputs')
        
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
                if args.bootstrap == 'no' or args.bootstrap == False or args.bootstrap == 'False' or args.bootstrap == 'false':
                    if i == 0:
                        f.write('Names\tPredictions\t%Identity_Nearest_VPOD_Sequence\n')
                    f.write(f"{names[i]}\t{predictions[i]}\t{per_iden_list[i]}\n")
                    print(f"{names[i]}\t{predictions[i]}\t{per_iden_list[i]}\n")
                    i+=1
                else:
                    if i == 0:
                        f.write('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\tStd_Deviation\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\tLmax_Hex_Color\n')
                        print('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\tStd_Deviation\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\n')
                        
                        # colors for hex_color_list generated from the mean prediction of the bootstraped predictions during the visulization steps    
                        hex_color_list = plot_prediction_subsets_with_CI(names, prediction_dict, mean_predictions, bootstrap_file, args.visualize_bootstrap)

                    f.write(f"{names[i]}\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\t{std_dev_list[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\t{hex_color_list[i]}\n")
                    print(f"{names[i]}\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\t{std_dev_list[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")
                    i+=1
        with open(log_file, 'w') as f:
            command_line_input = ' '.join(sys.argv)
            f.write(f"Command executed:\t{command_line_input}\n")
            f.write(f"Model Used:\t{args.model}\nEncoding Method:\t{args.encoding_method}\n")
            print(f"\nModel Used:\t{args.model}\nEncoding Method:\t{args.encoding_method}\n")
                            
        if hex_color_list:
            write_to_excel(names, predictions, per_iden_list, excel_output, 
                mean_predictions, median_predictions, ci_lowers, 
                ci_uppers, std_dev_list, hex_color_list, seq_lens_list)
                        
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
    
if __name__ == "__main__":
    main()
