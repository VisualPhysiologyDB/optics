import os
import sys
import subprocess
import argparse
import datetime
import tempfile
import pathlib
import itertools
import json
import time
import contextlib
import warnings

import joblib
from joblib import Parallel, delayed
import shap
import xgboost
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from tqdm import tqdm

from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from optics_scripts.utils import extract_fasta_entries

# --- Progress Bar Context Manager (Same as optics_predictions) ---
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

# --- Worker Function for Multiprocessing ---
def _worker_process_shap_sequence(name, sequence, model_path, alignment_data, wrk_dir, cached_prediction=None):
    """
    Worker function to process a single sequence:
    1. Runs MAFFT alignment.
    2. Encodes the sequence using the model's preprocessor.
    3. Generates a prediction (only if cached_prediction is None).
    """
    # Load model locally in worker to avoid serialization issues/memory bloat
    loaded_mod = load_obj(model_path)
    
    process_id = os.getpid()
    safe_name = "".join(c for c in name if c.isalnum())
    temp_dir = os.path.join(wrk_dir, 'tmp')
    temp_seq_path = os.path.join(temp_dir, f"{process_id}_{safe_name}_seq.fasta")
    temp_ali_path = os.path.join(temp_dir, f"{process_id}_{safe_name}_ali.fasta")
    
    try:
        # Write temp input file
        with open(temp_seq_path, "w") as f:
            f.write(f">{name}\n{sequence}")

        # Try MAFFT execution
        mafft_executables = [
            'mafft',
            str(pathlib.Path(wrk_dir) / 'optics_scripts/mafft/mafft-win/mafft.bat'),
            str(pathlib.Path(wrk_dir) / 'optics_scripts/mafft/mafft-mac/mafft.bat')
        ]
        
        alignment_successful = False
        for exe in mafft_executables:
            try:
                cmd = [exe, '--add', temp_seq_path, '--keeplength', alignment_data]
                with open(temp_ali_path, 'w') as f_out:
                    subprocess.run(cmd, stdout=f_out, stderr=subprocess.PIPE, check=True, text=True)
                alignment_successful = True
                break
            except (FileNotFoundError, subprocess.CalledProcessError):
                continue

        if not alignment_successful:
            return sequence, None # Signal failure

        # Process Alignment and Predict
        seq_type = 'aa'
        prediction = None
        encoded_test = None
        
        # Load reference just to get the shape for slicing (optimization: could pass length, but this is safe)
        ref_copy = read_data(alignment_data, seq_type=seq_type, is_main=True, gap_threshold=0.50)
        last_seq = int(ref_copy.shape[0])

        for gap_thresh in [0.5, 0.501, 0.505, 0.51, 0.515, 0.520]:
            try:
                new_seq_test = read_data(temp_ali_path, seq_type=seq_type, is_main=True, gap_threshold=gap_thresh)
                new_seq_test = new_seq_test.iloc[last_seq:].copy()

                if new_seq_test.empty:
                    continue

                # Encode the sequence (Required for SHAP)
                encoded_test = loaded_mod[0].transform(new_seq_test)
                
                # Check cache for prediction
                if cached_prediction is not None:
                    # Use the cached value, wrap in list/array to match model output format
                    prediction = [cached_prediction]
                else:
                    # Run the prediction step
                    prediction = loaded_mod.predict(new_seq_test)
                    
                break 
            except Exception:
                continue

        if prediction is None:
            return sequence, None

        # Return format: (sequence_key, data_dict)
        # We include 'encoded_seq' here which is NOT in the standard cache, 
        # but is needed for SHAP. We will separate them later.
        return sequence, {
            'name': name,
            'prediction': round(float(prediction[0]), 1),
            'encoded_seq': encoded_test,
            'len': len(sequence)
        }

    except Exception as e:
        return sequence, None
    finally:
        # Cleanup
        if os.path.exists(temp_seq_path):
            os.remove(temp_seq_path)
        if os.path.exists(temp_ali_path):
            os.remove(temp_ali_path)

def generate_shap_explanation(
    input_file, pred_dir=None, output='unnamed_shap_comparison', save_as='svg',
    model="whole-dataset", encoding_method='aa_prop',
    model_version='vpod_1.3', cmd_line=None, 
    mode='both', use_reference_sites=True, n_jobs=-1
):
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    clean_output_name = output.replace('.tsv', '').replace('.txt', '').replace('.csv', '')
    
    if pred_dir is None:
        report_dir = f'./prediction_outputs/optics_shap_{clean_output_name}_{dt_label}'
    else:
        report_dir = f'{pred_dir}/optics_shap_on_{clean_output_name}_{dt_label}'
    
    os.makedirs(report_dir, exist_ok=True)

    # Logging
    if cmd_line and report_dir:
        log_file_path = os.path.join(report_dir, 'shap_script_log.txt')
        with open(log_file_path, 'a') as f:
            f.write(f"Command executed:\n{cmd_line}\n\n")

    # Path setup
    script_path = pathlib.Path(__file__).resolve()
    wrk_dir = str(script_path.parent).replace('\\', '/')
    if not os.path.isdir(f'{wrk_dir}/tmp'):
        os.makedirs(f'{wrk_dir}/tmp')
        
    data_dir = f"{wrk_dir}/data"
    
    # Model Dictionaries
    model_datasets = {
        "whole-dataset": f"{data_dir}/fasta/{model_version}/wds_aligned_VPOD_1.3_het.fasta",
        "wildtype": f"{data_dir}/fasta/{model_version}/wt_aligned_VPOD_1.3_het.fasta",
        "vertebrate": f"{data_dir}/fasta/{model_version}/vert_aligned_VPOD_1.3_het.fasta",
        "invertebrate": f"{data_dir}/fasta/{model_version}/inv_aligned_VPOD_1.3_het.fasta",
        "wildtype-vert": f"{data_dir}/fasta/{model_version}/wt_vert_aligned_VPOD_1.3_het.fasta",
        "type-one": f"{data_dir}/fasta/{model_version}/Karyasuyama_T1_ops_aligned.fasta",
        "whole-dataset-mnm": f"{data_dir}/fasta/{model_version}/wds_mnm_aligned_VPOD_1.3_het.fasta",
        "wildtype-mnm": f"{data_dir}/fasta/{model_version}/wt_mnm_aligned_VPOD_1.3_het.fasta",
        "vertebrate-mnm": f"{data_dir}/fasta/{model_version}/vert_mnm_aligned_VPOD_1.3_het.fasta",
        "invertebrate-mnm": f"{data_dir}/fasta/{model_version}/inv_mnm_aligned_VPOD_1.3_het.fasta",
        "wildtype-vert-mnm": f"{data_dir}/fasta/{model_version}/wt_vert_mnm_aligned_VPOD_1.3_het.fasta",
        "wildtype-mut": f"{data_dir}/fasta/{model_version}/wt_mut_added_aligned_VPOD_1.3_het.fasta",
    }
    
    imp_report_directory = f"{data_dir}/importance_reports/{model_version}/{encoding_method}/{model}_translated_importance_report.csv"
    if model == "invertebrate" or model == "invertebrate-mnm":
        ref_seq_name = "Squid"
    elif model == "type-one":
        ref_seq_name = "Microbe"
    else:
        ref_seq_name = "Bovine"
    
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

    if input_file is None:
        raise Exception('Error: No file given')

    # Load sequences
    names, sequences = extract_fasta_entries(input_file)
    if len(names) == 0:
        raise Exception("No sequences found in input file.")

    # --- Caching Setup ---
    # We use 'reg_models' type for SHAP usually, matching the model dict
    model_type = 'reg_models' 
    cache_dir = f"{wrk_dir}/data/cached_predictions/{model_type}/{model_version}/{encoding_method}"
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = f"{cache_dir}/{model}_pred_dict.json"

    try:
        with open(cache_file, 'r') as f:
            cached_pred_dict = json.load(f)
        print('\nPrediction cache file successfully loaded.\n')
    except (json.JSONDecodeError, FileNotFoundError):
        cached_pred_dict = {}
        print('\nCache file not found or invalid. A new cache will be created.\n')

    # --- Parallel Processing ---
    alignment_data = model_datasets[model]
    model_path = model_directories[model]
    
     # Prepare list of sequences to process
    # We create a tuple: (name, sequence, cached_prediction_value_or_None)
    sequences_to_process = []
    
    # Track how many we found in cache for the user's info
    cached_count = 0
    
    for name, seq in zip(names, sequences):
        cached_val = None
        if seq in cached_pred_dict:
            # Retrieve the 'single_prediction' from the cache dict
            # Expected format: {'len': 500, 'single_prediction': 495.2, ...}
            cached_val = cached_pred_dict[seq].get('single_prediction')
            cached_count += 1
        
        sequences_to_process.append((name, seq, cached_val))

    print(f"Found {cached_count} sequences in cache. Predictions for these will be retrieved.\n")
    print(f"Processing {len(sequences_to_process)} sequences (alignment & encoding required for all)...\n")
        
    with tqdm_joblib(tqdm(total=len(sequences_to_process), desc="Processing Sequences", unit="seq", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]")) as pbar:
        # We pass the cached_val to the worker
        mp_results = Parallel(n_jobs=n_jobs, backend='loky')(
            delayed(_worker_process_shap_sequence)(
                name, seq, model_path, alignment_data, wrk_dir, cached_val
            ) for name, seq, cached_val in sequences_to_process
        )

    # Process results from workers
    newly_predicted_count = 0
    valid_processed_data = []

    for seq_key, result_data in mp_results:
        if result_data:
            valid_processed_data.append(result_data)
            
            # Update cache only if it wasn't there before
            if seq_key not in cached_pred_dict:
                cached_pred_dict[seq_key] = {
                    'len': result_data['len'],
                    'single_prediction': result_data['prediction']
                }
                newly_predicted_count += 1
        else:
            print(f"Failed to process sequence: {seq_key[:20]}...")

    # Update Cache File
    if newly_predicted_count > 0:
        print(f"Saving {newly_predicted_count} new predictions to cache file...")
        try:
            with open(cache_file, 'w') as f:
                json.dump(cached_pred_dict, f, indent=4)
        except Exception as e:
            print(f"Error: Could not save cache file: {e}")

    if not valid_processed_data:
        raise Exception("No sequences were successfully processed.")

    # --- Setup SHAP Explainer (Once) ---
    print("\nInitializing SHAP Explainer (Loading reference dataset)...")
    # Load the model in main thread for SHAP
    loaded_mod = load_obj(model_path)
    
    # We must load the background data for the TreeExplainer
    # We do this once here instead of inside the loop
    ref_seq_type = 'aa'
    ref_df = read_data(alignment_data, seq_type=ref_seq_type, is_main=True, gap_threshold=0.50)
    encoded_refs = loaded_mod[0].transform(ref_df)
    
    explainer = shap.TreeExplainer(loaded_mod[-1], encoded_refs)
    base_value = explainer.expected_value
    print(f"SHAP Base Value (Average Predicted λmax): {float(base_value):.1f} nm")

    # --- Calculate SHAP Values ---
    print(f"Calculating SHAP values for {len(valid_processed_data)} sequences...")
    
    # We iterate and calculate SHAP.
    final_data_for_plotting = []
    
    for item in tqdm(valid_processed_data, desc="Computing SHAP Values", unit="seq", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]"):
        encoded_seq = item['encoded_seq']
        shap_values = explainer.shap_values(encoded_seq)
        
        final_data_for_plotting.append({
            'name': item['name'],
            'prediction': item['prediction'],
            'encoded_seq': encoded_seq,
            'shap_values': shap_values
        })

    # --- Visualization Modes ---
    if use_reference_sites == True and model!="type-one":
        site_translation_df = pd.read_csv(imp_report_directory)
        site_translation_dict = dict(zip(site_translation_df['feature'],site_translation_df['true_position']))
    else:
        site_translation_dict = None

    # MODE: SINGLE
    if mode in ['single', 'both']:
        print("\n--- Generating Individual SHAP Explanations ---")
        for data in final_data_for_plotting:
            save_single_plot(
                data, base_value, report_dir, save_as, encoding_method, site_translation_dict, ref_seq_name
            )

    # MODE: COMPARISON
    if mode in ['comparison', 'both']:
        if len(final_data_for_plotting) < 2:
            print("\n[WARNING] Comparison mode requires at least 2 sequences. Skipping comparisons.")
        else:
            pairs = list(itertools.combinations(final_data_for_plotting, 2))
            print(f"\n--- Generating Pairwise Comparisons ({len(pairs)} pairs) ---")
            
            for seq1_data, seq2_data in pairs:
                save_comparison_plot(
                    seq1_data, seq2_data, report_dir, save_as, encoding_method, 
                    cmd_line, log_file_path if cmd_line else None, site_translation_dict, ref_seq_name
                )

    print(f"\nAnalysis complete. Results saved to: {report_dir}")


def save_single_plot(data, base_value, report_dir, save_as, encoding_method, site_translation_dict, ref_seq_name):
    """Generates and saves a bar plot explaining the prediction of a single sequence."""
    
    name = data['name']
    prediction = data['prediction']
    shap_vals = data['shap_values'][0]
    features = data['encoded_seq']
    
    if isinstance(site_translation_dict, dict):
        if encoding_method == "one_hot":
            trans_site_list = []
            for feature in features.columns:
                try:
                    trans_site_list.append(f"{site_translation_dict[int(feature.split('_')[0][1:])]:.0f}_{feature.split('_')[1]}")
                except:
                    #print(f'Exception was triggered for {feature}')
                    trans_site_list.append(feature)
        else:
            trans_site_list = [f"{site_translation_dict[feature[1:]]:.0f}_{feature.split('_')[1]}" for feature in features.columns]
        
        df = pd.DataFrame({
            'feature': features.columns,
            'state': features.values[0],
            'shap_effect': shap_vals,
            'true_position': trans_site_list
        })
    else:
        df = pd.DataFrame({
            'feature': features.columns,
            'state': features.values[0],
            'shap_effect': shap_vals,
            'true_position': features.columns # Fallback
        })
    
    df['abs_shap'] = df['shap_effect'].abs()
    df = df.sort_values(by='abs_shap', ascending=False).head(15) 
    df = df.sort_values(by='shap_effect')

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = ["#30c898d9" if x < 0 else "#f4a365ec" for x in df['shap_effect']]
    
    if isinstance(site_translation_dict, dict):
        bars = ax.barh(df['true_position'], df['shap_effect'], color=colors)
    else:
        bars = ax.barh(df['feature'], df['shap_effect'], color=colors)

    for i, bar in enumerate(bars):
        state = df['state'].iloc[i]
        val_str = f"{state:.0f}" if encoding_method == 'one_hot' else f"{state:.2f}"
        ax.text(0, bar.get_y() + bar.get_height()/2, f' {val_str} ',
                va='center', ha='right' if bar.get_width() < 0 else 'left',
                color='black', fontweight='bold', fontsize=9)

    ax.axvline(0, color='grey', linewidth=2, linestyle='--')
    ax.set_xlabel('Impact on Predicted λmax (nm)', fontsize=12)
    
    if isinstance(site_translation_dict, dict):
        ax.set_ylabel(f'Amino Acid Position ({ref_seq_name} Reference)', fontsize=12)
    else: 
        ax.set_ylabel('Feature', fontsize=12)

    ax.set_title(f'SHAP Explanation: {name}\nPrediction: {prediction:.2f}nm | Base Shape Value: {float(base_value):.2f}nm')
    ax.tick_params(axis='y', labelrotation=45)

    
    if encoding_method == 'one_hot':
        legend_text = "1 = Present // 0 = Absent"

        anchored_text = AnchoredText(legend_text, loc='lower right',
                                    prop=dict(size=10, ha='right'), frameon=True,
                                    pad=0.5, borderpad=0.5)
        anchored_text.patch.set_boxstyle("round,pad=0.3")
        anchored_text.patch.set_facecolor("whitesmoke")
        anchored_text.patch.set_edgecolor("grey")
        anchored_text.patch.set_linewidth(0.5)
        anchored_text.patch.set_alpha(0.8)
        ax.add_artist(anchored_text)
    
    safe_name = name.replace('|', '_').replace('/', '_')
    plt.tight_layout()
    plt.savefig(f'{report_dir}/{safe_name}_individual_shap.{save_as}', dpi=300, format=save_as)
    plt.close()


def save_comparison_plot(seq1, seq2, report_dir, save_as, encoding_method, cmd_line=None, log_path=None, site_translation_dict=None, ref_seq_name="Bovine"):
    """Generates and saves a comparative SHAP plot (Sequence 1 vs Sequence 2)."""    
    
    name1 = seq1['name']
    name2 = seq2['name']
    
    shap_diff = seq1['shap_values'][0] - seq2['shap_values'][0]
    
    feat1_vals = seq1['encoded_seq'].values[0]
    feat2_vals = seq2['encoded_seq'].values[0]
    cols = seq1['encoded_seq'].columns
    
    mask = (feat1_vals != feat2_vals)
    
    if isinstance(site_translation_dict, dict):
        if encoding_method == "one_hot":
            trans_site_list = []
            for feature in np.array(cols)[mask]:
                try:
                    trans_site_list.append(f"{site_translation_dict[int(feature.split('_')[0][1:])]:.0f}_{feature.split('_')[1]}")
                except:
                    print(f'Exception was triggered for {feature}')
                    trans_site_list.append(feature)          
        else:
            trans_site_list = [f"{site_translation_dict[feature[1:]]:.0f}_{feature.split('_')[1]}" for feature in np.array(cols)[mask]]
        
        comparison_df = pd.DataFrame({
            'feature': np.array(cols)[mask],
            f'{name1}_states': feat1_vals[mask],
            f'{name2}_states': feat2_vals[mask],
            f'{name1}_shap': seq1['shap_values'][0][mask],
            f'{name2}_shap': seq2['shap_values'][0][mask],
            'shap_difference': shap_diff[mask], 
            'true_position': trans_site_list
        })
    else:
        comparison_df = pd.DataFrame({
            'feature': np.array(cols)[mask],
            f'{name1}_states': feat1_vals[mask],
            f'{name2}_states': feat2_vals[mask],
            f'{name1}_shap': seq1['shap_values'][0][mask],
            f'{name2}_shap': seq2['shap_values'][0][mask],
            'shap_difference': shap_diff[mask]
        })

    if encoding_method == 'one_hot':
        comparison_df[f'{name1}_states'] = comparison_df[f'{name1}_states'].astype(int)
        comparison_df[f'{name2}_states'] = comparison_df[f'{name2}_states'].astype(int)

    comparison_df['abs_shap_diff'] = comparison_df['shap_difference'].abs()
    
    safe_n1 = name1.replace('|', '_')
    safe_n2 = name2.replace('|', '_')
    csv_name = f'{safe_n1}_vs_{safe_n2}_shap_data.csv'
    comparison_df.sort_values(by='abs_shap_diff', ascending=False).to_csv(f'{report_dir}/{csv_name}', index=False)
    
    pred_diff = seq1['prediction'] - seq2['prediction']
    log_msg = (f"\nSHAP Comparison: {name1} vs {name2}\n"
               f"\tSum SHAP diff: {shap_diff.sum():.2f}nm\n"
               f"\tActual Pred diff: {pred_diff:.2f}nm")
    print(log_msg)
    if log_path:
        with open(log_path, 'a') as f:
            f.write(log_msg)

    comparison_df = comparison_df.sort_values(by='abs_shap_diff', ascending=False).head(10)
    comparison_df = comparison_df.sort_values(by='shap_difference')

    if comparison_df.empty:
        return
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ["#30c898d9" if x < 0 else "#f4a365ec" for x in comparison_df['shap_difference']]
    
    if isinstance(site_translation_dict, dict):
        bars = ax.barh(comparison_df['true_position'], comparison_df['shap_difference'], color=colors)
    else:
        bars = ax.barh(comparison_df['feature'], comparison_df['shap_difference'], color=colors)


    for i, bar in enumerate(bars):
        aa_a = comparison_df[f'{name1}_states'].iloc[i]
        aa_b = comparison_df[f'{name2}_states'].iloc[i]
        
        val_fmt = "{:.0f}" if encoding_method == 'one_hot' else "{:.2f}"
        label = f" {val_fmt.format(aa_a)} → {val_fmt.format(aa_b)} "
        
        ax.text(0, bar.get_y() + bar.get_height()/2, label,
                va='center', ha='right' if bar.get_width() < 0 else 'left',
                color='black', fontweight='bold', fontsize=10)

    ax.axvline(0, color='grey', linewidth=2, linestyle='--')
    ax.set_xlabel(f'Contribution to Difference in Predicted λmax (nm)', fontsize=12)
    if isinstance(site_translation_dict, dict):
        ax.set_ylabel(f'Amino Acid Position ({ref_seq_name} Reference)', fontsize=12)
    else:
        ax.set_ylabel('Feature', fontsize=12)
    
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    ax.set_ylim(ymin, ymax + 0.75)
    ax.text((xmin + xmax)/2, ymax + 0.5,
            f"{name1} → {name2}\n({pred_diff:.1f}nm Difference in Predicted λmax)",
            ha='center', va='top', fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', fc="whitesmoke", ec='grey', lw=1, alpha=0.9))

    ax.tick_params(axis='y', labelrotation=45)

    if encoding_method == 'one_hot':
        legend_text = "1 = Present // 0 = Absent"

        anchored_text = AnchoredText(legend_text, loc='lower right',
                                    prop=dict(size=10, ha='right'), frameon=True,
                                    pad=0.5, borderpad=0.5)
        anchored_text.patch.set_boxstyle("round,pad=0.3")
        anchored_text.patch.set_facecolor("whitesmoke")
        anchored_text.patch.set_edgecolor("grey")
        anchored_text.patch.set_linewidth(0.5)
        anchored_text.patch.set_alpha(0.8)
        ax.add_artist(anchored_text)

    plt.tight_layout()
    plt.savefig(f'{report_dir}/{safe_n1}_vs_{safe_n2}_viz.{save_as}', dpi=350, bbox_inches='tight', format=save_as)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Predict protein properties using OPTICS.")

    parser.add_argument("-i", "--input", help="Path to FASTA file", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save output", type=str, default=None)
    parser.add_argument("-p", "--prediction_prefix", help="Prefix for filenames", type=str, default="unnamed")
    parser.add_argument("--mode", 
                        help="Analysis mode: 'comparison' (pairwise), 'single' (individual), or 'both'.", 
                        type=str, default="both", choices=['comparison', 'single', 'both'])

    parser.add_argument("-v", "--model_version", help="Model version", type=str, default="vpod_1.3", choices=['vpod_1.3'])
    parser.add_argument("-m", "--model", help="Prediction model", type=str, default="whole-dataset")
    parser.add_argument("-e", "--encoding", help="Encoding method", type=str, default="aa_prop", choices=['one_hot', 'aa_prop'])
    parser.add_argument("--save_viz_as", help="File type", type=str, default="svg", choices=['svg','png','pdf'])
    parser.add_argument("--use_reference_sites", help="Use reference site numbering, instead of feature names", action="store_true")
    
    # New argument for multiprocessing
    parser.add_argument("--n_jobs",
                        help="Number of parallel processes to run. -1 (default) uses all available.", 
                        type=int,
                        default=-1)

    args = parser.parse_args()
    command_run = " ".join(sys.argv)

    generate_shap_explanation(
        args.input, args.output_dir, args.prediction_prefix, args.save_viz_as,
        args.model, args.encoding, args.model_version,
        cmd_line=command_run, mode=args.mode,
        use_reference_sites=args.use_reference_sites,
        n_jobs=args.n_jobs
    )