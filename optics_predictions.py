import subprocess 
import os
import json
import argparse
import time
import datetime
import tempfile
import pathlib
from multiprocessing import Manager
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import matplotlib
import warnings

matplotlib.use('Agg')  # Use 'Agg' to prevent Mac crash when using GUI

# The VPOD/OPTICS special sauce ~
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from optics_scripts.utils import extract_fasta_entries, write_to_excel
from optics_scripts.blastp_analysis import run_blastp_analysis
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

# Set of the 20 standard amino acid one-letter codes
STANDARD_AMINO_ACIDS = {
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
}
# Global flag to ensure the warning message is printed only once
_non_standard_amino_acid_warning_shown = False

def filter_non_standard_aa(sequence: str) -> str:
    """
    Filters a protein sequence to keep only the 20 standard amino acids.

    Removes any non-standard amino acids (e.g., U, O, X, B, Z) and issues a
    one-time warning to the console if any are found.

    Args:
        sequence: The input protein sequence string.

    Returns:
        A new string containing only the standard amino acids.
    """
    global _non_standard_amino_acid_warning_shown
    
    # Using a list comprehension for a concise and performant filter
    filtered_sequence_list = [aa for aa in sequence.upper() if aa in STANDARD_AMINO_ACIDS]
    filtered_sequence = "".join(filtered_sequence_list)
    
    # Check if any non-standard amino acids were removed
    if len(filtered_sequence) != len(sequence) and not _non_standard_amino_acid_warning_shown:
        warnings.warn("Warning: Non-standard amino acids were detected in your input. These have been removed for the prediction. Only the 20 standard amino acids are supported.")
        _non_standard_amino_acid_warning_shown = True
        
    return filtered_sequence

### REFACTOR NOTE: This new worker function is designed to be self-contained.
### It loads its own models, preventing the memory duplication that caused crashes.
### It no longer uses Manager.dict(), instead returning a simple dictionary.
def _worker_predict_sequence(name, sequence, selected_model, bootstrap, wrk_dir, model_version, model_path, bs_model_folder_path, bootstrap_num, preload_to_memory):
    """
    Worker function executed by each parallel process.
    It loads its own models to prevent memory duplication in the main process.
    """
    time.sleep(0.1)
    ### REFACTOR NOTE: All model loading now happens HERE, inside the worker.
    # Each process gets its own copy without overwhelming the system.
    main_model = load_obj(model_path)
    
    # Initialize bootstrap models list.
    bootstrap_models = []
    if bootstrap and preload_to_memory:
        if os.path.isdir(bs_model_folder_path):
            model_files = sorted([f for f in os.listdir(bs_model_folder_path) if f.endswith('.pkl')])
            # Only load the number of models requested
            for i in range(min(bootstrap_num, len(model_files))):
                full_path = os.path.join(bs_model_folder_path, model_files[i])
                bootstrap_models.append(load_obj(full_path))
        if not bootstrap_models:
             print(f"Warning: Bootstrap enabled but no models found in {bs_model_folder_path} for worker process.")


    alignment_data_paths = {
        "whole-dataset": f"{wrk_dir}/data/fasta/{model_version}/wds_aligned_VPOD_1.2_het.fasta",
        "wildtype": f"{wrk_dir}/data/fasta/{model_version}/wt_aligned_VPOD_1.2_het.fasta",
        "vertebrate": f"{wrk_dir}/data/fasta/{model_version}/vert_aligned_VPOD_1.2_het.fasta",
        "invertebrate": f"{wrk_dir}/data/fasta/{model_version}/inv_only_aligned_VPOD_1.2_het.fasta",
        "wildtype-vert": f"{wrk_dir}/data/fasta/{model_version}/wt_vert_aligned_VPOD_1.2_het.fasta",
        "type-one": f"{wrk_dir}/data/fasta/{model_version}/Karyasuyama_T1_ops_aligned.fasta",
        "whole-dataset-mnm": f"{wrk_dir}/data/fasta/{model_version}/wds_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-mnm": f"{wrk_dir}/data/fasta/{model_version}/wt_mnm_aligned_VPOD_1.2_het.fasta",
        "vertebrate-mnm": f"{wrk_dir}/data/fasta/{model_version}/vert_mnm_aligned_VPOD_1.2_het.fasta",
        "invertebrate-mnm": f"{wrk_dir}/data/fasta/{model_version}/inv_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-vert-mnm": f"{wrk_dir}/data/fasta/{model_version}/wt_vert_mnm_aligned_VPOD_1.2_het.fasta",
        "wildtype-mut": f"{wrk_dir}/data/fasta/{model_version}/wt_mut_added_aligned_VPOD_1.2_het.fasta",
    }
    alignment_data = alignment_data_paths[selected_model]

    process_id = os.getpid()
    safe_name = "".join(c for c in name if c.isalnum()) # Clean the name for use in a filename
    temp_dir = os.path.join(wrk_dir, 'tmp')
    temp_seq_path = os.path.join(temp_dir, f"{process_id}_{safe_name}_seq.fasta")
    temp_ali_path = os.path.join(temp_dir, f"{process_id}_{safe_name}_ali.fasta")
    
    try:
        with open(temp_seq_path, "w") as f:
            f.write(f">{name}\n{sequence}")


        mafft_executables = [
            'mafft',
            str(pathlib.Path(wrk_dir) / 'optics_scripts/mafft/mafft-win/mafft.bat'),
            str(pathlib.Path(wrk_dir) / 'optics_scripts/mafft/mafft-mac/mafft.bat')
        ]
        
        alignment_successful = False
        last_error = ""
        for exe in mafft_executables:
            try:
                mafft_cmd = [exe, '--add', temp_seq_path, '--keeplength', alignment_data]
                with open(temp_ali_path, 'w') as f_out:
                    subprocess.run(mafft_cmd, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)
                alignment_successful = True
                break
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                last_error = str(e)
                continue

        if not alignment_successful:
            warnings.warn(f'ERROR: MAFFT alignment failed for {name}. Last error:\n{last_error}')
            return name, None # Return name to identify failure

        seq_type = 'aa'
        prediction = None
        new_seq_test = None
        ref_copy = read_data(alignment_data, seq_type=seq_type, is_main=True, gap_threshold=0.50)
        last_seq = int(ref_copy.shape[0])
        for gap_thresh in [0.5, 0.501, 0.505, 0.51, 0.515, 0.520]:
            try:
                new_seq_test = read_data(temp_ali_path, seq_type=seq_type, is_main=True, gap_threshold=gap_thresh)
                new_seq_test = new_seq_test.iloc[last_seq:].copy()
                                
                ### CHANGE NOTE: Added a defensive check here.
                # Predicting on an empty DataFrame can cause C-level crashes.
                if new_seq_test.empty:
                    # This can happen with very gappy sequences. It's not an error, just can't predict.
                    warnings.warn(f'You might have a gappy sequence causing an issue: {name}')
                    time.sleep(0.1)
                    continue
                
                prediction = main_model.predict(new_seq_test)
                break 
            except Exception:
                warnings.warn(f'Trying additional distance thresholds for sequence {name}')
                time.sleep(0.1)
                continue

        if prediction is None:
            warnings.warn(f"Failed to process sequence with all gap thresholds: {name}")
            return name, None
        
        # This dictionary will be returned and used to update the main cache
        result_dict = {
            'len': len(sequence),
            'single_prediction': round(float(prediction[0]), 1)
        }

        if bootstrap:
            ### REFACTOR NOTE: We pass the loaded bootstrap_models list here.
            # No Manager.dict() is needed for prediction_dict; it's temporary for this call.
            mean_pred, ci_low, ci_up, median_pred, std_dev, all_preds = calculate_ensemble_CI(
                prediction, bootstrap_models, new_seq_test, name, bootstrap_num, bs_model_folder_path
            )
            result_dict.update({
                'mean_prediction': round(float(mean_pred), 1),
                'ci_lower': round(float(ci_low), 1),
                'ci_upper': round(float(ci_up), 1),
                'median_prediction': round(float(median_pred), 1),
                'std_deviation': round(float(std_dev), 1),
                'all_bs_predictions': all_preds.tolist()
            })
        
        return sequence, result_dict

    finally:
        # Manual cleanup of the files created by this specific worker
        if os.path.exists(temp_seq_path):
            os.remove(temp_seq_path)
        if os.path.exists(temp_ali_path):
            os.remove(temp_ali_path)

def process_sequences_from_file(file, selected_model, identity_report, blastp, refseq, reffile, 
                                bootstrap, bootstrap_num, encoding_method, wrk_dir, model_version, preload_to_memory, n_jobs, tolerate_non_standard_aa=False):
    if file is None:
        raise ValueError('Error: No input file was provided.')
        
    names_unfiltered, sequences_unfiltered = extract_fasta_entries(file)
    
    # --- REFACTOR 1: Process all sequences first, tracking valid entries ---
    # This list will store all sequences that pass initial checks, preserving duplicates and order.
    all_valid_entries = []
    removed_sequences = []
    
    for name, seq_entry in zip(names_unfiltered, sequences_unfiltered):
        seq_body = seq_entry.split('\n', 1)[1].replace('\n', '')
        clean_seq_body = filter_non_standard_aa(seq_body)
        
        # Condition 1: Check for non-standard amino acids
        if (clean_seq_body != seq_body) and not tolerate_non_standard_aa:
            print(f'WARNING: Sequence {name} contained non-standard amino acids and will be skipped.')
            removed_sequences.append(name)
            continue # Skip to the next sequence

        # Condition 2: Check for valid length
        if not (300 <= len(clean_seq_body) <= 600):
            print(f'WARNING: Sequence {name} (length {len(clean_seq_body)}) is outside the 300-600 aa range and will be skipped.')
            if seq_body != clean_seq_body:
                print(f'NOTE: This sequence was originally {len(seq_body)} aa but was cleaned to {len(clean_seq_body)} aa.')
            removed_sequences.append(name)
            continue # Skip to the next sequence
            
        # If all checks pass, add it to our list of valid entries
        all_valid_entries.append({'name': name, 'sequence': clean_seq_body})

    if tolerate_non_standard_aa:
        print(f'\n{len(removed_sequences)} sequences were removed due to length constraints.')
    else:
        print(f'\n{len(removed_sequences)} sequences were removed due to length constraints and/or non-standard amino acids.')

    # Get the unique set of sequences that actually need prediction
    unique_seq_to_name_map = {}
    for entry in all_valid_entries:
        sequence = entry['sequence']
        if sequence not in unique_seq_to_name_map:
            unique_seq_to_name_map[sequence] = entry['name']
    print(f'Found {len(all_valid_entries)} valid sequences to process, corresponding to {len(unique_seq_to_name_map)} unique sequences.')

    data_dir = f"{wrk_dir}/data"
    # Dictionaries for paths remain the same...
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
    
    raw_data = model_raw_data[selected_model]
    metadata = model_metadata[selected_model]
    blast_db = model_blast_db[selected_model]
    
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
    
    model_path = model_directories[selected_model]
    bs_model_folder_path = model_bs_dirs.get(selected_model, '') # Use .get for safety

    # --- Caching Logic ---
    cache_dir = f"{wrk_dir}/data/cached_predictions"
    os.makedirs(cache_dir, exist_ok=True)
    model_type = 'bs_models' if bootstrap else 'reg_models'
    cache_file = f"{cache_dir}/{model_type}/{model_version}/{encoding_method}/{selected_model}_pred_dict.json"

    try:
        with open(cache_file, 'r') as f:
            cached_pred_dict = json.load(f)
        print('\nPrediction cache file successfully loaded.\n')
    except (json.JSONDecodeError, FileNotFoundError):
        cached_pred_dict = {}
        print('\nCache file not found or invalid. A new cache will be created.\n')
    
    prediction_results = {}
    # The list for multiprocessing will now store tuples of (name, sequence)
    sequences_for_mp = []

    # Iterate through the unique sequences using the new map
    for seq, name in unique_seq_to_name_map.items():
        if seq in cached_pred_dict:
            prediction_results[seq] = cached_pred_dict[seq]
        else:
            # Add the real name and sequence to the list for the worker
            sequences_for_mp.append((name, seq))
            
    print(f"{len(prediction_results)} unique sequences found in cache. Predicting {len(sequences_for_mp)} new unique sequences.")

    if sequences_for_mp:
        try:
            with tqdm_joblib(tqdm(total=len(sequences_for_mp), desc="Processing New Sequences", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]")) as pbar:
                # The worker now gets the real name and the sequence.
                # IMPORTANT: The worker MUST return the sequence as the first element in the tuple,
                # as the sequence is our unique key for mapping results back.
                # Expected return format: (sequence, result_data)
                mp_results = Parallel(n_jobs=n_jobs, backend='loky')(
                    delayed(_worker_predict_sequence)(
                        name, # The actual name of the first-seen sequence
                        seq,  # The sequence string
                        selected_model, bootstrap, wrk_dir, model_version, model_path, 
                        bs_model_folder_path, bootstrap_num, preload_to_memory
                    ) for name, seq in sequences_for_mp
                )
            
            newly_predicted_count = 0
            for seq_key, result_data in mp_results:
                if result_data:
                    prediction_results[seq_key] = result_data # Store result
                    cached_pred_dict[seq_key] = result_data   # Update cache object
                    newly_predicted_count += 1
                else:
                    # This case is tricky; if a sequence fails, we need to remove all entries with it.
                    print(f"Failed to get prediction for a sequence, it will be excluded from final results.")
                    # Mark this sequence as failed in our results dictionary
                    prediction_results[seq_key] = None
            
            if newly_predicted_count > 0:
                print(f"Saving {newly_predicted_count} new predictions to cache file...")
                try:
                    with open(cache_file, 'w') as f:
                        json.dump(cached_pred_dict, f, indent=4)
                except Exception as e:
                    print(f"Error: Could not save cache file: {e}")

        except Exception as e:
            print(f"A critical error occurred during multiprocessing: {e}")
            # Attempt to save any results we managed to get before the error
            if any(prediction_results):
                with open(cache_file, 'w') as f:
                    json.dump(cached_pred_dict, f, indent=4)
            raise

    # --- REFACTOR 3: Assemble final results by iterating through the original valid entries ---
    names, predictions, mean_predictions, ci_lowers, ci_uppers, median_predictions, std_dev_list, seq_lens = [], [], [], [], [], [], [], []
    prediction_dict = {}
    final_blast_entries = [] # For running BLASTp only on the final set

    for entry in all_valid_entries:
        name = entry['name']
        sequence = entry['sequence']
        
        # Look up the result for this sequence
        result = prediction_results.get(sequence)
        
        if result:
            names.append(name)
            final_blast_entries.append({'name': name, 'sequence': sequence})
            
            predictions.append(result.get('single_prediction'))
            seq_lens.append(result.get('len'))
            if bootstrap:
                mean_predictions.append(result.get('mean_prediction'))
                ci_lowers.append(result.get('ci_lower'))
                ci_uppers.append(result.get('ci_upper'))
                median_predictions.append(result.get('median_prediction'))
                std_dev_list.append(result.get('std_deviation'))
                if 'all_bs_predictions' in result:
                    prediction_dict[name] = np.array(result['all_bs_predictions'])
        else:
            # This sequence failed prediction or was not processed; add its name to the removed list.
            if name not in removed_sequences:
                removed_sequences.append(name)

    per_iden_list = ['N/A'] * len(names)
    if blastp and names:
        # Extract unique sequences and their first corresponding names for BLAST
        unique_blast_seqs = {item['sequence']: item['name'] for item in reversed(final_blast_entries)}
        blast_sequences = list(unique_blast_seqs.keys())
        blast_names = list(unique_blast_seqs.values())

        blastp_results_df = run_blastp_analysis(blast_sequences, blast_names, blast_db, raw_data, metadata, 
                                                output_file=identity_report, 
                                                ref_seq_id=refseq, reffile=reffile, wrk_dir=wrk_dir, n_jobs=n_jobs)
        
        if not blastp_results_df.empty:
            # Map results back to all sequences
            blast_map = blastp_results_df.set_index('query_id')['percent_identity'].to_dict()
            name_to_seq_map = {item['name']: item['sequence'] for item in final_blast_entries}
            seq_to_blast_result = {seq: blast_map.get(name) for seq, name in unique_blast_seqs.items()}
            
            per_iden_list = [seq_to_blast_result.get(name_to_seq_map.get(name), 'N/A') for name in names]


    return names, mean_predictions, ci_lowers, ci_uppers, prediction_dict, predictions, median_predictions, per_iden_list, std_dev_list, seq_lens, removed_sequences


def run_optics_predictions(input_sequence, pred_dir=None, output='optics_predictions',
                           model="whole-dataset", encoding_method='aa_prop', blastp=True,
                           iden_report='blastp_report.txt', refseq='bovine', reffile=None,
                           bootstrap=True, bootstrap_num = 100, visualize_bootstrap=True, bootstrap_viz_file='bootstrap_viz', save_as='svg', full_spectrum_xaxis=False,
                           model_version='vpod_1.3', preload_to_memory=False, n_jobs=-1, tolerate_non_standard_aa=False):
    ### REFACTOR NOTE: The 'preload_to_memory' argument is now ignored as it's an anti-pattern.
    #if preload_to_memory:
    #    print("Warning: --preload_bootstrap_models is deprecated and has no effect. Models are now loaded within each worker process to ensure stability.")

    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    script_path = pathlib.Path(__file__).resolve()
    wrk_dir = str(script_path.parent).replace('\\', '/')
    
    models = ['whole-dataset', 'wildtype', 'vertebrate', 'invertebrate', 'wildtype-vert', 'type-one', 'whole-dataset-mnm', 'wildtype-mnm', 'vertebrate-mnm', 'invertebrate-mnm', 'wildtype-vert-mnm', 'wildtype-mut']
    if model not in models:
        raise ValueError(f"Invalid model choice. Must be one of {models}")
    if bootstrap_num > 100:
        print("Warning: bootstrap_num cannot exceed 100. Setting to 100.")
        bootstrap_num = 100

    # Directory setup
    os.makedirs(f'{wrk_dir}/tmp', exist_ok=True)
    os.makedirs(f'{wrk_dir}/data/cached_blastp_analysis', exist_ok=True)
    
    if pred_dir:
        os.makedirs(pred_dir, exist_ok=True)
    else:
        pred_dir = './prediction_outputs'
        os.makedirs(pred_dir, exist_ok=True)
    
    output_prefix = output.replace('.tsv', '').replace('.txt', '').replace('.csv', '')
    report_dir = f'{pred_dir}/optics_on_{output_prefix}_{dt_label}'
    os.makedirs(report_dir, exist_ok=True)

    blastp_file = f'{report_dir}/{iden_report}'.replace('.tsv', '').replace('.txt', '') + '.csv'
    bootstrap_file_path = f'{report_dir}/{bootstrap_viz_file}'
    log_file = f'{report_dir}/arg_log.txt'

    temp_input_file_path = None
    if not os.path.isfile(input_sequence):
        # Treat as a single sequence string, write to a temp file
        temp_input_file = tempfile.NamedTemporaryFile(mode="w", dir=f'{wrk_dir}/tmp', suffix=".fasta", delete=False)
        temp_input_file.write(f">input_sequence\n{input_sequence}\n")
        temp_input_file.close()
        input_sequence_path = temp_input_file.name
        temp_input_file_path = input_sequence_path
    else:
        input_sequence_path = input_sequence

    try:
        names, mean_preds, ci_lows, ci_ups, pred_dict, preds, median_preds, iden_list, std_devs, seq_lens_list, removed = process_sequences_from_file(
            input_sequence_path, model, blastp_file, blastp, refseq, reffile, bootstrap, 
            bootstrap_num, encoding_method, wrk_dir, model_version, preload_to_memory, n_jobs, tolerate_non_standard_aa
        )

        output_path = f'{report_dir}/{output_prefix}_predictions.tsv'
        excel_output = f'{report_dir}/{output_prefix}_predictions_for_excel.xlsx'

        if not names:
            print("No sequences were successfully processed. No output files will be generated.")
            return None, None

        # Write main results file
        with open(output_path, 'w') as f:
            if not bootstrap:
                colors = [wavelength_to_rgb(p) for p in preds if p is not None]
                hex_color_list = [matplotlib.colors.to_hex(c) for c in colors]
                write_to_excel(names, preds, iden_list, excel_output, hex_color_list=hex_color_list, seq_lens_list=seq_lens_list)
                pred_df = pd.DataFrame({
                    'Names': names,
                    'Single_Prediction': preds,
                    '%Identity_Nearest_VPOD_Sequence': iden_list,
                    'Sequence_Length': seq_lens_list,
                    'Lmax_Hex_Color': hex_color_list
                })
                pred_df.to_csv(f, sep='\t', index=False)
            else:
                hex_color_list = plot_prediction_subsets_with_CI(names, pred_dict, mean_preds,
                                                                bootstrap_file_path, visualize_bootstrap, save_as=save_as, full_spectrum_xaxis=full_spectrum_xaxis)
                write_to_excel(names, preds, iden_list, excel_output,
                                mean_preds, median_preds, ci_lows,
                                ci_ups, std_devs, hex_color_list, seq_lens_list)
                pred_df = pd.DataFrame({
                    'Names': names,
                    'Single_Prediction': preds,
                    'Prediction_Means': mean_preds,
                    'Prediction_Medians': median_preds,
                    'Prediction_Lower_Bounds': ci_lows,
                    'Prediction_Upper_Bounds': ci_ups,
                    'Std_Deviation': std_devs,
                    '%Identity_Nearest_VPOD_Sequence': iden_list,
                    'Sequence_Length': seq_lens_list,
                    'Lmax_Hex_Color': hex_color_list
                })
                pred_df.to_csv(f, sep='\t', index=False)
        
        # ... (rest of the file writing, logging, and color annotations remain the same) ...
        # Write a text file for sequences/ids removed from the optics analysis
        if len(removed) > 0:
            print(f'Saving text file with accessions/ids {len(removed)} removed sequence(s)')
            with open(f'{report_dir}/removed_sequences.txt', 'w') as f:
                for seq_name in removed:
                    f.write(f"{seq_name}\n")

        # Write log file and color annotations...
        with open(log_file, 'w') as f:
            arg_string = (f"input_sequence: {input_sequence}\nreport_dir: {report_dir}\n"
                        f"output_file: {output}\nmodel: {model}\nencoding_method: {encoding_method}\n"
                        f"blastp: {blastp}\nblastp_report: {iden_report}\nrefseq: {refseq}\n"
                        f"custom_ref_file: {reffile}\nbootstrap: {bootstrap}\n"
                        f"\tvisualize_bootstrap: {visualize_bootstrap}\n\t\tbootstrap_viz_file: {bootstrap_viz_file}\n\t\tsave_as: {save_as}\n\t\tfull_spectrum_xaxis: {full_spectrum_xaxis}")
            f.write(f"Selected Options...\n{arg_string}\n")
            print(f"\nModel Used:\t{model}\nEncoding Method:\t{encoding_method}\n")
        
        if 'hex_color_list' in locals() and hex_color_list:
            with open(f'{report_dir}/fig_tree_color_annotation.txt', 'w') as g:
                g.write("Name\t!color\n")
                for name, hex_color in zip(names, hex_color_list):
                    g.write(f"{name}\t{hex_color}\n")

            with open(f'{report_dir}/itol_color_annotation.txt', 'w') as g:
                g.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")
                for name, hex_color in zip(names, hex_color_list):
                    g.write(f"{name}\tlabel_background\t{hex_color}\n")
        
        print(f"Predictions Complete! Results are in: {report_dir}")
        return pred_df, output_path

    finally:
        # Clean up the temporary file if it was created
        if temp_input_file_path and os.path.exists(temp_input_file_path):
            os.remove(temp_input_file_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Predict protein properties using OPTICS.")
    # Arguments remain the same...
    parser.add_argument("-i", "--input", 
                        help="Either a single sequence or a path to a FASTA file", 
                        type=str, 
                        required=True)
    parser.add_argument("-o", "--output_dir", 
                        help="Desired directory to save output folder/files (optional).", 
                        type=str, 
                        default=None,
                        required=False)
    parser.add_argument("-p", "--prediction_prefix", 
                        help="Base filename for prediction outputs (optional).", 
                        type=str, 
                        default="unnamed", 
                        required=False)
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
    blastp_group = parser.add_argument_group("BLASTp analysis (optional)")
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
                            type=str)
    bootstrap_group = parser.add_argument_group("Bootstrap analysis (optional)")
    bootstrap_group.add_argument("--bootstrap", 
                                help="Enable bootstrap predictions.", 
                                action="store_true")
    bootstrap_group.add_argument("--bootstrap_num", 
                                help="Number of bootstrap models to load for prediction replicates. Default and max is 100", 
                                type=int,
                                default=100)
    bootstrap_group.add_argument("--preload_bootstrap_models", 
                                help="DEPRECATED: This flag has no effect. Models are always loaded in worker processes for stability.", 
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

    run_optics_predictions(args.input, args.output_dir,
                        args.prediction_prefix, args.model, args.encoding,
                        args.blastp, args.blastp_report, args.refseq, args.custom_ref_file,
                        args.bootstrap, args.bootstrap_num, args.visualize_bootstrap, args.bootstrap_viz_file, args.save_viz_as, 
                        args.full_spectrum_xaxis, args.model_version, args.preload_bootstrap_models, args.n_jobs)
