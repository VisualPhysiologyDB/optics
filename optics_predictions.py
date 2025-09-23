import subprocess 
import os
import sys
import json
import argparse
import datetime
import tempfile
import pathlib
from multiprocessing import Manager
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use 'Agg' to prevent Mac crash when using GUI

# The VPOD/OPTICS special sauce ~
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from optics_scripts.utils import extract_fasta_entries, write_to_excel
from optics_scripts.blastp_align import seq_sim_report
from optics_scripts.bootstrap_predictions import calculate_ensemble_CI, plot_prediction_subsets_with_CI, wavelength_to_rgb

import contextlib
from tqdm import tqdm
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

def process_sequence(sequence=None, name=None, selected_model=None, identity_report=None, blastp=None, refseq=None, reffile=None, 
                     bootstrap=None, prediction_dict=None, wrk_dir='', only_blast=False, model_version='vpod_1.3', 
                     loaded_mod=None, loaded_bs_models=None, preload_to_memory=False, bs_model_folder_path=None, bootstrap_num=100):
    
    """Processes a single opsin amino acid sequence for lmax prediction.

    This helper function orchestrates the processing of a single input opsin amino acid sequence.
    It determines file paths based on the selected model and encoding method,
    optionally performs a BLASTp search against a curated blastp database, 
    aligns the sequence to a reference alignment using MAFFT, and predicts 
    lmax using a pre-trained regression model. It can also perform
    bootstrap analysis to estimate prediction confidence intervals.

    It handles temporary file creation for sequence input and alignment output,
    and attempts to run MAFFT using system-installed, Windows, or macOS
    executables bundled with the application.

    Args:
        sequence (str, optional): The input amino acid sequence. Can be a
            raw sequence string or a string containing a full FASTA entry
            (header starting with '>'). Defaults to None.
        name (str, optional): The identifier or name for the sequence, used
            primarily for tracking in reports and dictionaries. Defaults to None.
        selected_model (str, optional): A key indicating which pre-trained model
            and associated dataset files (alignment, raw data, metadata, BLAST DB,
            model files) to use. Valid keys are defined within the function's
            internal dictionaries (e.g., "whole-dataset", "wildtype", etc.).
            Defaults to None.
        identity_report (str, optional): Path to write a detailed BLASTp identity
            report file by the `seq_sim_report` function. Defaults to None.
        blastp (bool | str, optional): If True or not 'no'/'False', runs BLASTp
            via `seq_sim_report` to find the percent identity to the best hit
            in the selected model's BLAST database. Defaults to None.
        refseq (str, optional): Reference sequence identifier used by
            `seq_sim_report` for specific BLAST comparisons if needed.
            Defaults to None.
        reffile (str, optional): Path to a file containing reference sequences,
            used by `seq_sim_report` if needed for BLAST comparison.
            Defaults to None.
        bootstrap (bool, optional): If True, performs bootstrap analysis using
            pre-generated bootstrap models to calculate prediction mean,
            confidence intervals, median, and standard deviation. Defaults to None.
        prediction_dict (dict, optional): A dictionary, typically a
            `multiprocessing.Manager().dict()` passed from a calling function,
            used to aggregate raw bootstrap predictions across multiple sequence
            processes when `bootstrap` is True. This function *receives* the
            updated dict from `calculate_ensemble_CI` and *returns* it.
            Defaults to None.
        encoding_method (str, optional): Specifies the sequence encoding method
            used for model training ('one_hot' or 'aa_prop'). This determines
            which set of model files and bootstrap directories are used.
            Defaults to 'one_hot'.
        wrk_dir (str, optional): The working directory for the analysis. If '',
            it defaults to the directory containing the script. Used for locating
            data, models, and temporary files. Defaults to ''.
        only_blast (bool, optional): If True, performs only the BLASTp step
            (if `blastp` is enabled) and returns the percent identity. Skips
            alignment and prediction. Useful for quickly getting identity for
            cached sequences. Defaults to False.
        loaded_mod (object, optional): A pre-loaded model object. If provided,
            the function will skip loading the model from disk, improving performance
            and memory usage in parallel operations. Defaults to None.
        loaded_bs_models (dict, optional): A dictionary of pre-loaded bootstrap
            model objects

    Returns:
        tuple | str:
            - If `bootstrap` is True: A tuple containing:
                (mean_prediction (float), ci_lower (float), ci_upper (float),
                 prediction_dict (dict), single_best_prediction (float),
                 median_prediction (float), percent_identity (str),
                 std_deviation (float), all_bootstrap_predictions (list))
            - If `bootstrap` is False: A tuple containing:
                (single_best_prediction (float), percent_identity (str))
            - If `only_blast` is True: The percent identity (str).
            - If `sequence` or `selected_model` is None: An error message (str).

    Raises:
        Exception: If MAFFT alignment fails across all attempted execution
            methods (Linux, Windows, macOS).
        FileNotFoundError: If temporary files created during the process cannot
            be deleted (indicates potential filesystem issues).
        Exception: Can propagate exceptions from called functions like
            `seq_sim_report`, `read_data`, `load_obj`, `calculate_ensemble_CI`.
            
    Requires:
        - MAFFT executable accessible via system PATH, or bundled versions at
          '{wrk_dir}/optics_scripts/mafft/mafft-win/mafft.bat' (Windows) or
          '{wrk_dir}/optics_scripts/mafft/mafft-mac/mafft.bat' (macOS).
        - Pre-generated data files (FASTA alignments, raw sequences, metadata,
          BLAST databases) and model files (.pkl regression models, bootstrap
          model directories) organized as expected relative to `wrk_dir`.
        - External Python libraries: pathlib, tempfile, subprocess, os, numpy.
        - Helper functions/modules: `read_data`, `load_obj`, `seq_sim_report`,
          `calculate_ensemble_CI` (assumed to be defined elsewhere).
    """
    
    if wrk_dir == '':
        script_path = pathlib.Path(__file__).resolve()  # Get absolute path
        wrk_dir = str(script_path.parent).replace('\\', '/')

    data_dir = f"{wrk_dir}/data"
    model_datasets = {
        "whole-dataset": f"{data_dir}/fasta/{model_version}/wds_aligned_VPOD_1.2_het.fasta",
        "wildtype": f"{data_dir}/fasta/{model_version}/wt_aligned_VPOD_1.2_het.fasta",
        "vertebrate": f"{data_dir}/fasta/{model_version}/vert_aligned_VPOD_1.2_het.fasta",
        "invertebrate": f"{data_dir}/fasta/{model_version}/inv_only_aligned_VPOD_1.2_het.fasta",
        "wildtype-vert": f"{data_dir}/fasta/{model_version}/wt_vert_aligned_VPOD_1.2_het.fasta",
        "type-one": f"{data_dir}/fasta/{model_version}/Karyasuyama_T1_ops_aligned.fasta",
        "whole-dataset-mnm": f"{data_dir}/fasta/{model_version}/wds_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-mnm": f"{data_dir}/fasta/{model_version}/wt_mnm_aligned_VPOD_1.2_het.fasta",
        "vertebrate-mnm": f"{data_dir}/fasta/{model_version}/vert_mnm_aligned_VPOD_1.2_het.fasta",
        "invertebrate-mnm": f"{data_dir}/fasta/{model_version}/inv_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-vert-mnm": f"{data_dir}/fasta/{model_version}/wt_vert_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-mut": f"{data_dir}/fasta/{model_version}/wt_mut_added_aligned_VPOD_1.2_het.fasta",

    }   
    model_raw_data = {
        "whole-dataset": f"{data_dir}/fasta/{model_version}/wds.txt",
        "wildtype": f"{data_dir}/fasta/{model_version}/wt.txt",
        "vertebrate": f"{data_dir}/fasta/{model_version}/vert.txt",
        "invertebrate": f"{data_dir}/fasta/{model_version}/inv_only.txt",
        "wildtype-vert": f"{data_dir}/fasta/{model_version}/wt_vert.txt",
        "type-one": f"{data_dir}/fasta/{model_version}/Karyasuyama_T1_ops.txt",
        "whole-dataset-mnm": f"{data_dir}/fasta/{model_version}/wds_mnm.txt",
        "wildtype-mnm": f"{data_dir}/fasta/{model_version}/wt_mnm.txt",
        "vertebrate-mnm": f"{data_dir}/fasta/{model_version}/vert_mnm.txt",
        "invertebrate-mnm": f"{data_dir}/fasta/{model_version}/inv_mnm.txt",
        "wildtype-vert-mnm": f"{data_dir}/fasta/{model_version}/wt_vert_mnm.txt",
        "wildtype-mut": f"{data_dir}/fasta/{model_version}/wt_mut_added.txt",
    }   
    model_metadata = {
        "whole-dataset": f"{data_dir}/fasta/{model_version}/wds_meta.tsv",
        "wildtype": f"{data_dir}/fasta/{model_version}/wt_meta.tsv",
        "vertebrate": f"{data_dir}/fasta/{model_version}/vert_meta.tsv",
        "invertebrate": f"{data_dir}/fasta/{model_version}/inv_meta.tsv",
        "wildtype-vert": f"{data_dir}/fasta/{model_version}/wt_vert_meta.tsv",
        "type-one": f"{data_dir}/fasta/{model_version}/Karyasuyama_T1_ops_meta.tsv",
        "whole-dataset-mnm": f"{data_dir}/fasta/{model_version}/wds_mnm_meta.csv",
        "wildtype-mnm": f"{data_dir}/fasta/{model_version}/wt_mnm_meta.csv",
        "vertebrate-mnm": f"{data_dir}/fasta/{model_version}/vert_mnm_meta.csv",
        "invertebrate-mnm": f"{data_dir}/fasta/{model_version}/inv_mnm_meta.csv",
        "wildtype-vert-mnm": f"{data_dir}/fasta/{model_version}/wt_vert_mnm_meta.csv",
        "wildtype-mut": f"{data_dir}/fasta/{model_version}/wt_meta.tsv",
    }   
    model_blast_db = {
        "whole-dataset": f"{data_dir}/blast_dbs/{model_version}/wds_db",
        "wildtype": f"{data_dir}/blast_dbs/{model_version}/wt_db",
        "vertebrate": f"{data_dir}/blast_dbs/{model_version}/vert_db",
        "invertebrate": f"{data_dir}/blast_dbs/{model_version}/invert_db",
        "wildtype-vert": f"{data_dir}/blast_dbs/{model_version}/wt_vert_db",
        "type-one": f"{data_dir}/blast_dbs/{model_version}/t1_db",
        "whole-dataset-mnm": f"{data_dir}/blast_dbs/{model_version}/wds_mnm_db",
        "wildtype-mnm": f"{data_dir}/blast_dbs/{model_version}/wt_mnm_db",
        "vertebrate-mnm": f"{data_dir}/blast_dbs/{model_version}/vert_mnm_db",
        "invertebrate-mnm": f"{data_dir}/blast_dbs/{model_version}/inv_mnm_db",
        "wildtype-vert-mnm": f"{data_dir}/blast_dbs/{model_version}/wt_vert_mnm_db",
        "wildtype-mut": f"{data_dir}/blast_dbs/{model_version}/wt_db",
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
    
    with tempfile.NamedTemporaryFile(mode="w", dir=f"{wrk_dir}/tmp", suffix=".fasta", delete=False) as temp_seq_file:
        temp_seq = temp_seq_file.name
        if '>' in sequence:
            temp_seq_file.write(sequence)
        else:
            temp_seq_file.write(f">placeholder_name\n{sequence}")
    
    if only_blast:
        percent_iden = '-'
        if blastp not in ['no', False, 'False']:
            percent_iden = seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile, wrk_dir)
        os.remove(temp_seq)
        return percent_iden

    percent_iden = '-'
    if blastp not in ['no', False, 'False']:
        percent_iden = seq_sim_report(temp_seq, name, refseq, blast_db, raw_data, metadata, identity_report, reffile, wrk_dir)

    
    # perform temporary alignment to training data for preprocessing
    with tempfile.NamedTemporaryFile(mode="w", dir=f"{wrk_dir}/tmp", suffix=".fasta", delete=False) as temp_ali_file:
        new_ali = temp_ali_file.name 
    try:
        cmd = ['mafft', '--add', temp_seq, '--keeplength', alignment_data]
        with open(new_ali, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        try:
            mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-win/mafft.bat'
            cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
            with open(new_ali, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
        except (subprocess.CalledProcessError, FileNotFoundError): 
            try:
                mafft_exe = f'{wrk_dir}/optics_scripts/mafft/mafft-mac/mafft.bat'
                cmd = [mafft_exe, '--add', temp_seq, '--keeplength', alignment_data]
                with open(new_ali, 'w') as f:
                    subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                raise Exception(f'MAFFT alignment failed. Check FASTA formatting.\n{e.stderr.decode()}')
            
    seq_type = 'aa'
    prediction = None
    for gap_thresh in [0.5, 0.501, 0.505, 0.51]:
        try:
            new_seq_test = read_data(new_ali, seq_type=seq_type, is_main=True, gap_threshold=gap_thresh)
            ref_copy = read_data(alignment_data, seq_type=seq_type, is_main=True, gap_threshold=gap_thresh)
            last_seq = int(ref_copy.shape[0])
            new_seq_test = new_seq_test.iloc[last_seq:].copy()
            
            # Use the pre-loaded model object
            if loaded_mod:
                prediction = loaded_mod.predict(new_seq_test)
                break # Exit loop on success

        except Exception:
            continue # Try next threshold
    
    if prediction is None:
        raise Exception("Failed to process sequence with all gap thresholds.")

    os.remove(new_ali)
    os.remove(temp_seq)

    if bootstrap and preload_to_memory:
        # Use the pre-loaded bootstrap models
        if not loaded_bs_models:
             raise ValueError("Bootstrap is enabled, but no pre-loaded bootstrap models were provided.")
        mean_prediction, ci_lower, ci_upper, prediction_dict, median_prediction, std_dev, predictions_all = calculate_ensemble_CI(prediction, loaded_bs_models, new_seq_test, name, prediction_dict, bootstrap_num=bootstrap_num)
        return (round(float(mean_prediction),1), round(float(ci_lower),1), round(float(ci_upper),1), prediction_dict, round(float(prediction[0]),1), round(float(median_prediction),1), str(percent_iden), round(float(std_dev),1), predictions_all.tolist())
    elif bootstrap:
        mean_prediction, ci_lower, ci_upper, prediction_dict, median_prediction, std_dev, predictions_all = calculate_ensemble_CI(prediction, loaded_bs_models, new_seq_test, name, prediction_dict, bs_model_folder_path, bootstrap_num)
        return (round(float(mean_prediction),1), round(float(ci_lower),1), round(float(ci_upper),1), prediction_dict, round(float(prediction[0]),1), round(float(median_prediction),1), str(percent_iden), round(float(std_dev),1), predictions_all.tolist())
    else:
        return (round(float(prediction[0]),1), str(percent_iden))
 

def process_sequences_from_file(file, selected_model, identity_report, blastp, refseq, reffile, 
                                bootstrap, bootstrap_num, encoding_method, wrk_dir, model_version, preload_to_memory, n_jobs):
    # Extract sequences and their names from the input fasta file
    if file == None:
        raise Exception('Error: No file given')
    names,sequences = extract_fasta_entries(file)       
    
    model_dir = f"{wrk_dir}/models"
    if encoding_method == 'aa_prop':
        model_directories = {
            "whole-dataset": f"{model_dir}/reg_models/{model_version}/aa_prop/wds_xgb.pkl",
            "wildtype": f"{model_dir}/reg_models/{model_version}/aa_prop/wt_gbr.pkl",
            "vertebrate": f"{model_dir}/reg_models/{model_version}/aa_prop/vert_xgb.pkl",
            "invertebrate": f"{model_dir}/reg_models/{model_version}/aa_prop/inv_gbr.pkl",
            "wildtype-vert": f"{model_dir}/reg_models/{model_version}/aa_prop/wt_vert_gbr.pkl",
            "type-one": f"{model_dir}/reg_models/{model_version}/aa_prop/t1_xgb.pkl",
            "whole-dataset-mnm": f"{model_dir}/reg_models/{model_version}/aa_prop/wds_mnm_xgb.pkl",
            "wildtype-mnm": f"{model_dir}/reg_models/{model_version}/aa_prop/wt_mnm_gbr.pkl",
            "vertebrate-mnm": f"{model_dir}/reg_models/{model_version}/aa_prop/vert_mnm_xgb.pkl",
            "invertebrate-mnm": f"{model_dir}/reg_models/{model_version}/aa_prop/inv_mnm_gbr.pkl",
            "wildtype-vert-mnm": f"{model_dir}/reg_models/{model_version}/aa_prop/wt_vert_mnm_xgb.pkl",
            "wildtype-mut": f"{model_dir}/reg_models/{model_version}/aa_prop/wt_mut_gbr.pkl",
        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wds_H1_H2_H3_P2_V_SCT_PKA_bootstrap_100_2025-03-23_22-05-50",
            "wildtype": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_H2_H3_P1_NCI_PKA_bootstrap_100_2025-03-24_10-19-45",
            "vertebrate": f"{model_dir}/bs_models/{model_version}/{encoding_method}/vert_H2_H3_NCI_SCT_PKB_bootstrap_100_2025-03-21_17-47-40",
            "invertebrate": f"{model_dir}/bs_models/{model_version}/{encoding_method}/inv_H1_H3_bootstrap_100_2025-03-21_17-40-08",
            "wildtype-vert": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_vert_H2_P2_V_MASS_bootstrap_100_2025-03-21_17-25-47",
            "type-one": f"{model_dir}/bs_models/{model_version}/{encoding_method}/t1_H3_P1_PKB_bootstrap_100_2025-03-24_10-31-03",
            "whole-dataset-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wds_mnm_H2_H3_NCI_MASS_bootstrap_100_2025-04-15_09-35-46",
            "wildtype-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_mnm_H1_H2_H3_NCI_MASS_PKA_bootstrap_100_2025-04-15_10-17-06",
            "vertebrate-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/vert_mnm_H3_P2_SCT_PKA_PKB_bootstrap_100_2025-04-15_11-00-24",
            "invertebrate-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/inv_mnm_H1_P1_SCT_bootstrap_100_2025-04-15_10-52-20",
            "wildtype-vert-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_vert_mnm_H2_H3_PKB_bootstrap_100_2025-04-15_10-40-35",
            "wildtype-mut": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_mut_H1_H2_H3_SCT_bootstrap_100_2025-04-25_17-14-36",
        }
        
    else:
        model_directories = {
            "whole-dataset": f"{model_dir}/reg_models/{model_version}/one_hot/wds_xgb.pkl",
            "wildtype": f"{model_dir}/reg_models/{model_version}/one_hot/wt_xgb.pkl",
            "vertebrate": f"{model_dir}/reg_models/{model_version}/one_hot/vert_xgb.pkl",
            "invertebrate": f"{model_dir}/reg_models/{model_version}/one_hot/invert_BayesianRidge.pkl",
            "wildtype-vert": f"{model_dir}/reg_models/{model_version}/one_hot/wt_vert_xgb.pkl",
            "type-one": f"{model_dir}/reg_models/{model_version}/one_hot/t1_xgb.pkl",
            "whole-dataset-mnm": f"{model_dir}/reg_models/{model_version}/one_hot/wds_mnm_xgb.pkl",
            "wildtype-mnm": f"{model_dir}/reg_models/{model_version}/one_hot/wt_mnm_gbr.pkl",
            "vertebrate-mnm": f"{model_dir}/reg_models/{model_version}/one_hot/vert_mnm_xgb.pkl",
            "invertebrate-mnm": f"{model_dir}/reg_models/{model_version}/one_hot/invert_mnm_gbr.pkl",
            "wildtype-vert-mnm": f"{model_dir}/reg_models/{model_version}/one_hot/wt_vert_mnm_xgb.pkl",
        }
        
        model_bs_dirs = {
            "whole-dataset": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wds_bootstrap",
            "wildtype": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_bootstrap",
            "vertebrate": f"{model_dir}/bs_models/{model_version}/{encoding_method}/vert_bootstrap",
            "invertebrate": f"{model_dir}/bs_models/{model_version}/{encoding_method}/invert_bootstrap",
            "wildtype-vert": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_vert_bootstrap",
            "type-one": f"{model_dir}/bs_models/{model_version}/{encoding_method}/t1_bootstrap",
            "whole-dataset-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wds_mnm_bootstrap",
            "wildtype-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_mnm_bootstrap",
            "vertebrate-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/vert_mnm_bootstrap",
            "invertebrate-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/invert_mnm_bootstrap",
            "wildtype-vert-mnm": f"{model_dir}/bs_models/{model_version}/{encoding_method}/wt_vert_mnm_bootstrap",
        }

    cache_dir = f"{wrk_dir}/data/cached_predictions"
    model_type = 'bs_models' if bootstrap else 'reg_models'
    cache_file = f"{cache_dir}/{model_type}/{model_version}/{encoding_method}/{selected_model}_pred_dict.json"

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
        

    # Load the base model once before starting parallel jobs to reduce memory overhead.
    model_path = model_directories[selected_model]
    main_model_preloaded = load_obj(model_path)
    
    # Load all bootstrap models once, if requested, to reduce memory overhead
    bs_models_preloaded = []
    if bootstrap and preload_to_memory:
        bs_model_folder_path = model_bs_dirs[selected_model]
        if os.path.isdir(bs_model_folder_path):
            # Sort to ensure deterministic order
            model_files = sorted([f for f in os.listdir(bs_model_folder_path) if f.endswith('.pkl')])
            for i in range(bootstrap_num):
                filename = model_files[i]
                full_path = os.path.join(bs_model_folder_path, filename)
                bs_models_preloaded.append(load_obj(full_path))
        if not bs_models_preloaded:
            print(f"Warning: Bootstrap is enabled but no models found in {bs_model_folder_path}")


    manager = Manager()
    prediction_dict = manager.dict()  # Use a shared dictionary for prediction outputs
    mp_cached_pred_dict = manager.dict(cached_pred_dict)  # Initialize a multiprocess available cached dict with cached data
    
    def process_sequence_wrapper(seq, name):  # Helper function
        just_seq = seq.split('\n', 1)[1]
        just_seq = just_seq.replace('\n', '')
        if bootstrap == False:
            if just_seq not in list(mp_cached_pred_dict.keys()):
                prediction, percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, wrk_dir, model_version=model_version, loaded_mod=main_model_preloaded)
                mp_cached_pred_dict[just_seq] = {'len': len(just_seq), 'single_prediction': prediction}
                return len(just_seq), prediction, percent_iden, None, None, None, None, None  # Consistent return values
            else:
                percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, wrk_dir, only_blast=True, model_version=model_version)
                return mp_cached_pred_dict[just_seq]['len'], mp_cached_pred_dict[just_seq]['single_prediction'], percent_iden, None, None, None, None, None  # Consistent return values

        else:
            if just_seq not in list(mp_cached_pred_dict.keys()):
                mean_prediction, ci_lower, ci_upper, updated_prediction_dict, prediction, median_prediction, percent_iden, std_dev, predictions_all = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, wrk_dir, model_version=model_version, loaded_mod=main_model_preloaded, loaded_bs_models=bs_models_preloaded, bs_model_folder_path=bs_model_folder_path, bootstrap_num=bootstrap_num)
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
                percent_iden = process_sequence(seq, name, selected_model, identity_report, blastp, refseq, reffile, bootstrap, prediction_dict, wrk_dir, only_blast=True, model_version=model_version)
                prediction_dict[name] = np.array(mp_cached_pred_dict[just_seq]['all_bs_predictions'])
                return mp_cached_pred_dict[just_seq]['len'], mp_cached_pred_dict[just_seq]['single_prediction'], percent_iden, mp_cached_pred_dict[just_seq]['mean_prediction'], mp_cached_pred_dict[just_seq]['ci_lower'], mp_cached_pred_dict[just_seq]['ci_upper'], mp_cached_pred_dict[just_seq]['median_prediction'], mp_cached_pred_dict[just_seq]['std_deviation']


    try:
        with tqdm_joblib(tqdm(total=len(sequences), desc="Processing Sequences", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]", 
                              dynamic_ncols=True, colour="#CF9FFF",
                              unit ='seqs',ascii="░▒▓")) as pbar:  # Use tqdm for progress bar
            results = Parallel(n_jobs=n_jobs)(delayed(process_sequence_wrapper)(seq, names[i]) for i, seq in enumerate(sequences))
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
        
        raise Exception(f"Error: {e}...\n")

    seq_lens, predictions, per_iden_list, mean_predictions, ci_lowers, ci_uppers, median_predictions, std_dev_list = [], [], [], [], [], [], [], []
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

def run_optics_predictions(input_sequence, pred_dir=None, output='optics_predictions',
                           model="whole-dataset", encoding_method='aa_prop', blastp=True,
                           iden_report='blastp_report.txt', refseq='bovine', reffile=None,
                           bootstrap=True, bootstrap_num = 100, visualize_bootstrap=True, bootstrap_viz_file='bootstrap_viz', save_as='svg', full_spectrum_xaxis=False,
                           model_version='vpod_1.3', preload_to_memory=False, n_jobs=-1):
    """
    Processes sequences using a selected model and generates prediction outputs.  This function
    encapsulates the logic from the original `main` function, making it callable from
    other scripts or notebooks.

    Args:
        input_sequence (str): Either a single sequence or a path to a FASTA file.
        report_dir (str, optional):  Name of folder directory to create.  Defaults to a
            timestamped directory if None.  If a directory name is provided, a timestamp
            will be appended to avoid overwrites, and all results put in a 'prediction_outputs' subdirectory.
        output (str, optional): Name for output file. Defaults to 'optics_predictions.txt'.
        model (str, optional): Model to use for prediction. Defaults to "whole-dataset".
        encoding_method (str, optional): Encoding method. Defaults to 'aa_prop'.
        blastp (bool, optional): Enable blastp analysis. Defaults to True.
        iden_report (str, optional): Blastp report output file name. Defaults to 'blastp_report.txt'.
        refseq (str, optional): Reference sequence for blastp. Defaults to 'bovine'.
        reffile (str, optional): Custom reference sequence file. Defaults to None.
        bootstrap (bool, optional): Enable bootstrap predictions. Defaults to True.
        visualize_bootstrap (bool, optional): Enable visualization of bootstrap predictions. Defaults to True.
        bootstrap_viz_file (str, optional): Output file name for bootstrap visualization. Defaults to 'bootstrap_viz'.

    Returns:
        tuple: A tuple containing the following lists (in order):
            - names: List of sequence names.
            - mean_predictions: List of mean predictions (if bootstrapping).
            - ci_lowers: List of lower confidence interval bounds (if bootstrapping).
            - ci_uppers: List of upper confidence interval bounds (if bootstrapping).
            - prediction_dict: Dictionary of all bootstrap predictions (if bootstrapping).
            - predictions: List of single predictions (or first prediction from bootstrap).
            - median_predictions: List of median predictions (if bootstrapping).
            - per_iden_list: List of percent identities from BLASTp.
            - std_dev_list: List of standard deviations (if bootstrapping).
            - seq_lens_list: list of sequence lengths

        Also creates several output files within the specified report directory:
            - Main results file (TSV or TXT)
            - Excel file with additional formatting
            - BLASTp report (if blastp is enabled)
            - Bootstrap visualization (PDF, if bootstrapping and visualization are enabled)
            - Argument log file
            - Color annotation files for FigTree and iTOL

    Raises:
        Exception: If `input_sequence` is not a file or a valid sequence.  (Improved error handling.)
    """
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    script_path = pathlib.Path(__file__).resolve()  # Get absolute path
    wrk_dir = str(script_path.parent).replace('\\', '/')
    #print(f"Script directory (pathlib): {wrk_dir}")

    # Argument validation (mimicking argparse behavior, but for function inputs)
    models = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert', 'type-one', 'whole-dataset-mnm', 'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 'wildtype-vert-mnm', 'wildtype-mut']
    encoding_methods = ['one_hot', 'aa_prop']
    ref_seq_choices = ['bovine', 'squid', 'microbe', 'custom']
    if bootstrap_num > 100:
        bootstrap_num = 100
    
    if model not in models:
        raise ValueError(f"Invalid model choice.  Must be one of {models}")
    if encoding_method not in encoding_methods:
        raise ValueError(f"Invalid encoding method choice.  Must be one of {encoding_methods}")
    if refseq not in ref_seq_choices:
        raise ValueError(f"Invalid refseq choice. Must be one of {ref_seq_choices}")
    # No need to check bool_choices.  We'll convert directly to boolean.

    # Directory setup (with added handling for pre-existing directories)
    if not os.path.isdir(f'{wrk_dir}/tmp'):
        os.makedirs(f'{wrk_dir}/tmp')
    
    if not pred_dir is None:
        if not os.path.isdir(pred_dir):
            os.makedirs(pred_dir)
    else:
        if not os.path.isdir('./prediction_outputs'):
            os.makedirs('./prediction_outputs') 
    
    output = output.replace('.tsv', '').replace('.txt', '').replace('.csv', '')
    if pred_dir is None:
        report_dir = f'./prediction_outputs/optics_on_unamed_{dt_label}'
    else:
        report_dir = f'{pred_dir}/optics_on_{output}_{dt_label}'
    os.makedirs(report_dir, exist_ok=True)  # exist_ok=True prevents errors if dir exists



    blastp_file = f'{report_dir}/{iden_report}'
    if not (blastp_file.endswith(('.txt','.tsv','.csv'))):
           blastp_file += '.txt'

    bootstrap_file = f'{report_dir}/{bootstrap_viz_file}'
    log_file = f'{report_dir}/arg_log.txt'

    # Input handling (file or sequence string)
    if os.path.isfile(input_sequence):
        names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens_list = process_sequences_from_file(input_sequence, model, blastp_file, blastp, refseq, reffile, bootstrap, bootstrap_num, encoding_method, wrk_dir, model_version, preload_to_memory, n_jobs)
        
        # Output file handling (TSV or TXT, Excel)
        if 'predictions' not in {output}:
            output += '_predictions'
        output_path = f'{report_dir}/{output}.tsv'
        excel_output = f'{report_dir}/{output}_for_excel.xlsx'

    else:  # Assume it's a single sequence
        #  create a temporary file.
        temp_input_file = os.path.join(f'{wrk_dir}/tmp', f'temp_input_{dt_label}.fasta')
        with open(temp_input_file, "w") as f:
            f.write(f">temp_seq\n{input_sequence}\n") #write temp file to be consistant with process_sequence_from_file function
        
        names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens_list = process_sequences_from_file(temp_input_file, model, blastp_file, blastp, refseq, reffile, bootstrap, bootstrap_num, encoding_method, wrk_dir, model_version, preload_to_memory, n_jobs)
        output_path = f'{report_dir}/{output}'
        if not output_path.endswith(('.tsv', '.txt')):
            output_path += '.tsv'
        excel_output = output_path.replace('.tsv', '_for_excel.xlsx').replace('.txt', '_for_excel.xlsx')
        os.remove(temp_input_file) #clean up temp file


    # Write main results file
    with open(output_path, 'w') as f:
        if not bootstrap: #cast to bool to simplify logic
            colors = [wavelength_to_rgb(pred) for pred in predictions]
            hex_color_list = [matplotlib.colors.to_hex(color) for color in colors]
            write_to_excel(names, predictions, per_iden_list, excel_output, hex_color_list=hex_color_list, seq_lens_list=seq_lens_list)
            # Make Prediction Dataframe
            pred_df = pd.DataFrame({
            'Names': names,
            'Single_Prediction': predictions,
            '%Identity_Nearest_VPOD_Sequence': per_iden_list,
            'Sequence_Length': seq_lens_list,
            'Lmax_Hex_Color': hex_color_list
            })
            # Write to file
            pred_df.to_csv(f, sep='\t', index=False)
            for i in range(len(names)):
                print(f"{names[i]}\t{predictions[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")
        else:
            hex_color_list = plot_prediction_subsets_with_CI(names, prediction_dict, mean_predictions,
                                                            bootstrap_file, visualize_bootstrap, save_as=save_as, full_spectrum_xaxis=full_spectrum_xaxis)
            write_to_excel(names, predictions, per_iden_list, excel_output,
                            mean_predictions, median_predictions, ci_lowers,
                            ci_uppers, std_dev_list, hex_color_list, seq_lens_list)
            # Make Prediction Dataframe
            pred_df = pd.DataFrame({
            'Names': names,
            'Single_Prediction': predictions,
            'Prediction_Means': mean_predictions,
            'Prediction_Medians': median_predictions,
            'Prediction_Lower_Bounds': ci_lowers,
            'Prediction_Upper_Bounds': ci_uppers,
            'Std_Deviation': std_dev_list,
            '%Identity_Nearest_VPOD_Sequence': per_iden_list,
            'Sequence_Length': seq_lens_list,
            'Lmax_Hex_Color': hex_color_list
            })
            # Write to file
            pred_df.to_csv(f, sep='\t', index=False)
            
            print('Names\tSingle_Prediction\tPrediction_Means\tPrediction_Medians\tPrediction_Lower_Bounds\tPrediction_Upper_Bounds\tStd_Deviation\t%Identity_Nearest_VPOD_Sequence\tSequence_Length\n')
            for i in range(len(names)):
                print(f"{names[i]}\t{predictions[i]}\t{mean_predictions[i]}\t{median_predictions[i]}\t{ci_lowers[i]}\t{ci_uppers[i]}\t{std_dev_list[i]}\t{per_iden_list[i]}\t{seq_lens_list[i]}\n")



    # Write log file
    with open(log_file, 'w') as f:
        #  We don't have sys.argv, so we reconstruct the equivalent.
        arg_string = (f"input_fasta: {input_sequence}\nreport_dir: {report_dir}\n"
                      f"output_file: {output}\nmodel: {model}\nencoding_method: {encoding_method}\n"
                      f"blastp: {blastp}\nblastp_report: {iden_report}\nrefseq: {refseq}\n"
                      f"custom_ref_file: {reffile}\nbootstrap: {bootstrap}\n"
                      f"\tvisualize_bootstrap: {visualize_bootstrap}\n\t\tbootstrap_viz_file: {bootstrap_viz_file}\n\t\tsave_as: {save_as}\n\t\tfull_spectrum_xaxis: {full_spectrum_xaxis}")
        
        if blastp == True:
            bp_cmd = f'--blastp --blastp_report {iden_report} --refseq {refseq} '
        else:
            bp_cmd = ''
        if reffile is not None and refseq == 'custom':
            rf_cmd = f"--custom_ref_file {reffile} "
        else:
            rf_cmd = ''
            
        if bootstrap == True:
            bs_cmd = '--bootstrap '
        else:
            bs_cmd = ''
        if visualize_bootstrap == True:
            vb_cmd = f'--visualize_bootstrap --bootstrap_viz_file {bootstrap_viz_file} --save_viz_as {save_as} --full_spectrum_xaxis {full_spectrum_xaxis}'
        else:
            vb_cmd = ''
    
        exec_cmd =  (f"python optics_predictions.py " 
                      f"-i {input_sequence} " 
                      f"-o {pred_dir} " 
                      f"-p {output} " 
                      f"-v {model_version} "
                      f"-m {model} " 
                      f"-e {encoding_method} " 
                      f"{bp_cmd}" 
                      f"{rf_cmd}" 
                      f"{bs_cmd}" 
                      f"{vb_cmd}\n")
        f.write(f"Selected Options...\n{arg_string}\n")  # More informative
        f.write(f"Command executed (reconstructed): {exec_cmd}\n")
        #f.write(f"Model Used:\t{model}\nEncoding Method:\t{encoding_method}\n")
        print(f"\nModel Used:\t{model}\nEncoding Method:\t{encoding_method}\n")

    # Write color annotation files
    with open(f'{report_dir}/fig_tree_color_annotation.txt', 'w') as g:
        g.write("Name\t!color\n")
        for name, hex_color in zip(names, hex_color_list):
            g.write(f"{name}\t{hex_color}\n")

    with open(f'{report_dir}/itol_color_annotation.txt', 'w') as g:
        g.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
        for name, hex_color in zip(names, hex_color_list):
            g.write(f"{name}\tlabel_background\t{hex_color}\n")

    print('Predictions Complete!')

    return pred_df, output_path
        
if __name__ == '__main__':
    # This part will only execute when the script is run directly (not imported)
    # It provides a simple way to test the module functionality from the command line
        
    parser = argparse.ArgumentParser(description="Predict protein properties using OPTICS.")

    # Input sequence or FASTA file
    parser.add_argument("-i", "--input", 
                        help="Either a single sequence or a path to a FASTA file", 
                        type=str, 
                        required=True)

    # Output directory for all results
    parser.add_argument("-o", "--output_dir", 
                        help="Desired directory to save output folder/files (optional).", 
                        type=str, 
                        default=None,
                        required=False)

    # Base filename for prediction folder and all subsequent results
    parser.add_argument("-p", "--prediction_prefix", 
                        help="Base filename for prediction outputs (optional).", 
                        type=str, 
                        default="unnamed", 
                        required=False)

    # Prediction model
    parser.add_argument("-v", "--model_version",
                        help="Version of models to use (optional).\nBased on the version of VPOD used to train models.", 
                        type=str, 
                        default="vpod_1.3",
                        choices=['vpod_1.3'],
                        required=False)
    parser.add_argument("-m", "--model", 
                        help="Prediction model to use (optional).", 
                        type=str, 
                        default="whole-dataset", 
                        choices=['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert', 'type-one', 'whole-dataset-mnm', 'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 'wildtype-vert-mnm', 'wildtype-mut'],
                        required=False)

    # Encoding method
    parser.add_argument("-e", "--encoding", 
                        help="Encoding method to use (optional).", 
                        type=str, 
                        default="aa_prop",
                        choices=['one_hot', 'aa_prop'],
                        required=False)
    parser.add_argument("--n_jobs",
                        help="Number of parallel processes to run.\n-1 is the default, utilizing all avaiable processors.", 
                        type=int,
                        default=-1)

    # BLASTp options
    blastp_group = parser.add_argument_group("BLASTp analysis (optional)") # Group related args

    blastp_group.add_argument("--blastp", 
                            help="Enable BLASTp analysis.", 
                            action="store_true")
    blastp_group.add_argument("--blastp_report", 
                            help="Filename for BLASTp report.", 
                            type=str, 
                            default="blastp_report.txt")
    blastp_group.add_argument("--refseq",
                            help="Reference sequence used for blastp analysis.", 
                            type=str, 
                            default="bovine",
                            choices=['bovine', 'squid', 'microbe', 'custom'])
    blastp_group.add_argument("--custom_ref_file",
                            help="Path to a custom reference sequence file for BLASTp.", 
                            type=str)  # No default, as it's optional

    # Bootstrap options
    bootstrap_group = parser.add_argument_group("Bootstrap analysis (optional)")

    bootstrap_group.add_argument("--bootstrap", 
                                help="Enable bootstrap predictions.", 
                                action="store_true")
    bootstrap_group.add_argument("--bootstrap_num", 
                                help="Number of bootstrap models to load for prediction replicates. Default and max is 100", 
                                type=int,
                                default=100)
    bootstrap_group.add_argument("--preload_bootstrap_models", 
                                help="Enable preloading of bootstrap models to memory.\nCan be quite cumbersome, but will theoretically make predictions faster.", 
                                action="store_true")
    bootstrap_group.add_argument("--visualize_bootstrap", 
                                help="Enable visualization of bootstrap predictions.", 
                                action="store_true")
    bootstrap_group.add_argument("--bootstrap_viz_file", 
                                help="Filename prefix for bootstrap visualizations.", 
                                type=str, 
                                default="bootstrap_viz")
    bootstrap_group.add_argument("--save_viz_as", 
                                help="File type for bootstrap visualizations (SVG, PNG, or PDF).", 
                                type=str, 
                                default="svg",
                                choices=['svg','png','pdf'])
    bootstrap_group.add_argument("--full_spectrum_xaxis", 
                            help="Enable visualization of predictions on a full spectrum x-axis (300-650nm).\nOtherwise, x-axis is scaled with predictions.", 
                            action="store_true")
    
    args = parser.parse_args()

    blastp_report_output = f"{args.blastp_report}" if args.blastp else None
    bootstrap_viz_output = f"{args.bootstrap_viz_file}" if args.visualize_bootstrap else None
    
    run_optics_predictions(args.input, args.output_dir,
                        args.prediction_prefix, args.model, args.encoding,
                        args.blastp, blastp_report_output, args.refseq, args.custom_ref_file,
                        args.bootstrap, args.bootstrap_num, args.visualize_bootstrap, bootstrap_viz_output, args.save_viz_as, 
                        args.full_spectrum_xaxis, args.model_version, args.preload_bootstrap_models, args.n_jobs)
