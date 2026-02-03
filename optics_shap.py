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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.offsetbox import AnchoredText
from tqdm import tqdm
from openpyxl import Workbook
from openpyxl.styles import PatternFill



from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from optics_scripts.utils import extract_fasta_entries

# --- Load physiochemical properties data ---
def load_amino_acid_properties():
    """Load both scaled and non-scaled amino acid properties from Excel file."""
    script_path = pathlib.Path(__file__).resolve()
    wrk_dir = str(script_path.parent).replace('\\', '/')
    prop_file = os.path.join(wrk_dir, 'data/aa_property_index/aa_property_values.xlsx')
    
    try:
        # Load both sheets
        scaled_df = pd.read_excel(prop_file, sheet_name='scaled')
        non_scaled_df = pd.read_excel(prop_file, sheet_name='non_scaled')
        
        # Create mapping dictionaries
        scaled_dict = {}
        non_scaled_dict = {}
        
        # Get property names from columns (excluding the first column which is amino acid)
        prop_names = scaled_df.columns[1:].tolist()  # Skip the first column (amino acid)
        
        for _, row in scaled_df.iterrows():
            aa = row.iloc[0]  # Amino acid code
            scaled_dict[aa] = {}
            for prop in prop_names:
                scaled_dict[aa][prop] = row[prop]
        
        for _, row in non_scaled_df.iterrows():
            aa = row.iloc[0]  # Amino acid code
            non_scaled_dict[aa] = {}
            for prop in prop_names:
                non_scaled_dict[aa][prop] = row[prop]
        
        return scaled_dict, non_scaled_dict, prop_names
        
    except Exception as e:
        print(f"Warning: Could not load amino acid properties file: {e}")
        print("Using scaled values for display as fallback.")
        return {}, {}, []

# Load amino acid properties globally
SCALED_AA_PROPS, NON_SCALED_AA_PROPS, PROPERTY_NAMES = load_amino_acid_properties()

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
        new_seq_test = None
        
        ali_names, ali_seqs = extract_fasta_entries(temp_ali_path)
        if ali_seqs:
            aligned_query_string = str(ali_seqs[-1]).split('\n')[1]
            aligned_ref_string = str(ali_seqs[0]).split('\n')[1]

            # Find the equivalent position in the aligned sequence
            ref_to_query_mapping = {}
            for aligned_idx, char in enumerate(aligned_ref_string):
                if (char != '-' and aligned_query_string[aligned_idx] != '-'):
                    ref_gaps = aligned_ref_string[:aligned_idx+1].count('-')                        
                    query_gaps = aligned_query_string[:aligned_idx+1].count('-')
                    ref_to_query_mapping.update({(aligned_idx+1-ref_gaps):(aligned_idx+1-query_gaps)})
                        
        
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
        return sequence, {
            'name': name,
            'prediction': round(float(prediction[0]), 1),
            'encoded_seq': encoded_test,
            'aligned_seq_df': new_seq_test,  # Raw aligned sequence (dataframe) for AA extraction
            'ref_to_query_map': ref_to_query_mapping,
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

def get_obs_aa(feature, aligned_df, encoding_method):
    """Helper to extract observed Amino Acid from the aligned dataframe based on feature name."""
    try:
        parts = feature.split('_')
        col_idx = parts[0]
        
        # Retrieve from the single row in aligned_df
        obs_aa = aligned_df.iloc[0][col_idx]
        if str(obs_aa) == 'nan':
            obs_aa = '-'
        return obs_aa

    except Exception as e:
        return "?"

def get_property_display_value(feature_name, encoded_value, observed_aa, encoding_method):
    """
    Get the display value for a property feature.
    For aa_prop encoding, returns non-scaled value if available, otherwise scaled.
    For one_hot encoding, returns the binary value.
    """
    if encoding_method != 'aa_prop':
        # For one_hot encoding, just return the binary value
        return encoded_value
    
    try:
        # Extract property name from feature name
        # Feature format: "P123_H1" or similar
        parts = feature_name.split('_')
        if len(parts) >= 2:
            property_name = parts[1]
            
            # Get non-scaled value if available
            if observed_aa in NON_SCALED_AA_PROPS and property_name in NON_SCALED_AA_PROPS[observed_aa]:
                return NON_SCALED_AA_PROPS[observed_aa][property_name]
            
            # Fall back to scaled value
            return encoded_value
        else:
            return encoded_value
    except Exception:
        return encoded_value

def generate_shap_explanation(
    input_file, pred_dir=None, output='unnamed_shap_comparison', save_as='svg',
    model="whole-dataset", encoding_method='aa_prop',
    model_version='vpod_1.3', cmd_line=None, 
    mode='both', n_positions=10, use_reference_sites=True, n_jobs=-1
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
    
    sequences_to_process = []
    cached_count = 0
    
    for name, seq in zip(names, sequences):
        cached_val = None
        if seq in cached_pred_dict:
            cached_val = cached_pred_dict[seq].get('single_prediction')
            cached_count += 1
        
        sequences_to_process.append((name, seq, cached_val))

    if cached_count >= 1:
        print(f"Found {cached_count} sequences in cache. Predictions for these will be retrieved.\n")
    else:
        print(f"Found 0 sequences in cache. Predictions for all sequences will be processed.\n")

    print(f"Processing {len(sequences_to_process)} sequences (alignment & encoding required for all)...\n")
        
    with tqdm_joblib(tqdm(total=len(sequences_to_process), desc="Processing Sequences", unit="seq", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]")) as pbar:
        mp_results = Parallel(n_jobs=n_jobs, backend='loky')(
            delayed(_worker_process_shap_sequence)(
                name, seq, model_path, alignment_data, wrk_dir, cached_val
            ) for name, seq, cached_val in sequences_to_process
        )    
    
    newly_predicted_count = 0
    valid_processed_data = []
    names = []
    predictions = []
    seq_lens_list = []
        
    for seq_key, result_data in mp_results:
        if result_data:
            valid_processed_data.append(result_data)
            
            names.append(result_data['name'])
            predictions.append(result_data['prediction'])
            seq_lens_list.append(result_data['len'])
            
            if seq_key not in cached_pred_dict:
                cached_pred_dict[seq_key] = {
                    'len': result_data['len'],
                    'single_prediction': result_data['prediction']
                }
                newly_predicted_count += 1
        else:
            print(f"Failed to process sequence: {seq_key[:20]}...")
    
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
    loaded_mod = load_obj(model_path)
    
    ref_seq_type = 'aa'
    ref_df = read_data(alignment_data, seq_type=ref_seq_type, is_main=True, gap_threshold=0.50)
    encoded_refs = loaded_mod[0].transform(ref_df)
    
    explainer = shap.TreeExplainer(loaded_mod[-1], encoded_refs)
    base_value = explainer.expected_value
    print(f"SHAP Base Value (Average Predicted λmax): {float(base_value):.1f} nm")

    # --- Load Translation Dictionary (Mandatory for True Positions) ---
    # We load this even if use_reference_sites is False for plotting, 
    # because we need it for the output CSV column.
    try:
        site_translation_df = pd.read_csv(imp_report_directory)
        site_translation_dict = dict(zip(site_translation_df['feature'], site_translation_df['true_position']))
        print("Site translation dictionary loaded for 'reference_position' mapping.")
    except Exception as e:
        print(f"Warning: Could not load site translation file ({imp_report_directory}). 'true_position' column will default to feature names.")
        site_translation_dict = None

    # --- Calculate SHAP Values ---
    print(f"Calculating SHAP values for {len(valid_processed_data)} sequences...")
    
    final_data_for_plotting = []
    shap_value_pred_list = []
    
    for item in tqdm(valid_processed_data, desc="Computing SHAP Values", unit="seq", bar_format="{l_bar}{bar:25}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}{postfix}]"):
        encoded_seq = item['encoded_seq']
        shap_values = explainer.shap_values(encoded_seq)
        shap_value_pred_list.append(shap_values[0].sum())
        final_data_for_plotting.append({
            'name': item['name'],
            'prediction': item['prediction'],
            'encoded_seq': encoded_seq,
            'aligned_seq_df': item['aligned_seq_df'],
            'ref_to_query_map': item['ref_to_query_map'],
            'shap_values': shap_values
        })

    just_prediction_data = pd.DataFrame({
        'Names': names,
        'Single_Prediction': predictions,
        'SHAP_Prediction': shap_value_pred_list,
        'Sequence_Length': seq_lens_list,
    })
    just_prediction_data.to_csv(f'{report_dir}/{clean_output_name}_predictions.csv', sep=',', index=False)
    
    # --- Prediction Matrix Generation ---
    if len(final_data_for_plotting) > 1:
        # (Same matrix logic as original)
        names_list = [d['name'] for d in final_data_for_plotting]
        preds_list = [d['prediction'] for d in final_data_for_plotting]
        n_seqs = len(names_list)
        diff_matrix = np.full((n_seqs, n_seqs), np.nan)
        for i in range(n_seqs):
            for j in range(n_seqs):
                if i <= j:
                    diff_matrix[i, j] = preds_list[i] - preds_list[j]
        diff_df = pd.DataFrame(diff_matrix, index=names_list, columns=names_list)
        diff_df.to_csv(f'{report_dir}/{clean_output_name}_prediction_diff_matrix.csv')

    # --- Visualization & Data Saving ---
    # MODE: SINGLE
    if mode in ['single', 'both']:
        print("\n--- Generating Individual SHAP Explanations ---")
        for data in final_data_for_plotting:
            save_single_plot(
                data, base_value, report_dir, save_as, encoding_method, 
                site_translation_dict, ref_seq_name, n_positions, use_reference_sites
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
                    cmd_line, log_file_path if cmd_line else None, 
                    site_translation_dict, ref_seq_name, n_positions, use_reference_sites
                )

    print(f"\nAnalysis complete. Results saved to: {report_dir}")


def save_single_plot(data, base_value, report_dir, save_as, encoding_method, site_translation_dict, ref_seq_name, n_positions, use_reference_sites_for_plot):
    """Generates and saves a bar plot explaining the prediction of a single sequence."""
    
    name = data['name']
    prediction = data['prediction']
    shap_vals = data['shap_values'][0]
    features = data['encoded_seq']
    aligned_df = data['aligned_seq_df']
    ref_to_query_dict = data['ref_to_query_map']
    
    # 1. Generate the 'reference_position' list (Reference Numbering)
    # This is critical for structure mapping, so we do it even if not plotting it.
    reference_position_col = []
    query_position_col = []
    
    if isinstance(site_translation_dict, dict):
        if encoding_method == "one_hot":
            for feature in features.columns:
                try:
                    # Extract numeric part P123_H -> 123
                    pos_idx = int(feature[1:].split('_')[0])
                    mapped_pos = site_translation_dict.get(pos_idx, feature)
                    query_pos = ref_to_query_dict.get(mapped_pos,'NA')
                    query_position_col.append(query_pos)
                    
                    # Format: 123_H
                    if isinstance(mapped_pos, (int, float)):
                         reference_position_col.append(f"{mapped_pos:.0f}_{feature.split('_')[1]}")
                    else:
                         reference_position_col.append(str(mapped_pos))
                except:
                    reference_position_col.append(feature)
                    query_position_col.append('NA')

        else:
            # AA Prop: P123_Hydrophobicity
            for feature in features.columns:
                try:
                    pos_idx = feature[1:] # P123 -> 123
                    mapped_pos = site_translation_dict.get(pos_idx, feature)
                    #print(f'this is the mapped poistion...{mapped_pos}')
                    query_pos = ref_to_query_dict.get(int(mapped_pos),'NA')
                    #print(f'this is the query position...{query_pos}')
                    query_position_col.append(query_pos)
                    
                     # Format: 123_Hydrophobicity
                    if isinstance(mapped_pos, (int, float)):
                        reference_position_col.append(f"{mapped_pos:.0f}_{feature.split('_')[1]}")
                    else:
                        reference_position_col.append(str(mapped_pos))
                except:
                    reference_position_col.append(feature)
                    query_position_col.append('NA')

    else:
        reference_position_col = features.columns.tolist()
        query_position_col.append('NA' for feats in reference_position_col)
        print('Site translation not working...')

    # Create DataFrame
    df = pd.DataFrame({
        'feature': features.columns,
        'state': features.values[0],
        'shap_effect': shap_vals,
        'reference_position': reference_position_col,
        'query_position': query_position_col
    })
    
    # Add Observed Amino Acid Column
    obs_aas = []
    for feat in df['feature']:
        obs_aas.append(get_obs_aa(feat, aligned_df, encoding_method))
    df['observed_aa'] = obs_aas
    
    # Add display values
    display_values = []
    for idx, row in df.iterrows():
        display_val = get_property_display_value(row['feature'], row['state'], row['observed_aa'], encoding_method)
        display_values.append(display_val)
    df['display_value'] = display_values

    df['abs_shap'] = df['shap_effect'].abs()
    
    # Save full data to CSV with observed AA and TRUE POSITIONS
    safe_name = name.replace('|', '_').replace('/', '_')
    df.sort_values(by='abs_shap', ascending=False).to_csv(f'{report_dir}/{safe_name}_shap_analysis.csv', index=False)

    # Filter for Plotting
    df = df.sort_values(by='abs_shap', ascending=False).head(n_positions) 
    df = df.sort_values(by='shap_effect')

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = ["#30c898d9" if x < 0 else "#f4a365ec" for x in df['shap_effect']]
    
    # Select Y-axis labels based on user preference
    if use_reference_sites_for_plot and isinstance(site_translation_dict, dict):
        y_labels = df['reference_position']
        y_axis_title = f'Amino Acid Position ({ref_seq_name} Reference)'
    else:
        y_labels = df['feature']
        y_axis_title = 'Feature'

    bars = ax.barh(y_labels, df['shap_effect'], color=colors)

    for i, bar in enumerate(bars):
        display_val = df['display_value'].iloc[i]
        obs_aa = df['observed_aa'].iloc[i]
        
        if encoding_method == 'aa_prop' and 'SCT' not in df['feature'].iloc[i]:
            try:
                float_val = float(display_val)
                val_str = f"{float_val:.2f}"
            except:
                val_str = str(display_val)
        else:
            val_str = f"{int(display_val)}"
        
        label_str = f" {val_str} ({obs_aa}) " if obs_aa != "?" else f" {val_str} "

        ax.text(0, bar.get_y() + bar.get_height()/2, label_str,
                va='center', ha='right' if bar.get_width() < 0 else 'left',
                color='black', fontweight='bold', fontsize=10)

    ax.axvline(0, color='grey', linewidth=2, linestyle='--')
    ax.set_xlabel(r'Impact on Predicted $λ_{max}$ (nm)', fontsize=13)
    ax.set_ylabel(y_axis_title, fontsize=13)

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    ax.set_ylim(ymin, ymax + 0.80)
    ax.text((xmin + xmax)/2, ymax + 0.55,
            f'SHAP Explanation: {name}\nPrediction: {prediction:.2f}nm | Base Shap Value: {float(base_value):.2f}nm',
            ha='center', va='top', fontsize=13,
            bbox=dict(boxstyle='round,pad=0.3', fc="whitesmoke", ec='grey', lw=1, alpha=0.9))
    
    ax.tick_params(axis='y', labelrotation=45, labelsize=12)
    ax.tick_params(axis='x', labelsize=12)
    
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
    plt.savefig(f'{report_dir}/{safe_name}_individual_shap.{save_as}', dpi=300, format=save_as)
    plt.close()


def save_comparison_plot(seq1, seq2, report_dir, save_as, encoding_method, cmd_line=None, log_path=None, site_translation_dict=None, ref_seq_name="Bovine", n_positions=10, use_reference_sites_for_plot=False):
    """Generates and saves a comparative SHAP plot."""    
    
    name1 = seq1['name']
    name2 = seq2['name']
    shap_diff = seq1['shap_values'][0] - seq2['shap_values'][0]
    
    feat1_vals = seq1['encoded_seq'].values[0]
    feat2_vals = seq2['encoded_seq'].values[0]
    cols = seq1['encoded_seq'].columns
    mask = (feat1_vals != feat2_vals)
    
    # Generate True Positions for masking
    masked_true_pos = []
    
    if isinstance(site_translation_dict, dict):
        masked_cols = np.array(cols)[mask]
        if encoding_method == "one_hot":
            for feature in masked_cols:
                try:
                    pos_idx = int(feature[1:].split('_')[0])
                    mapped_pos = site_translation_dict.get(pos_idx, feature)
                    if isinstance(mapped_pos, (int, float)):
                        masked_true_pos.append(f"{mapped_pos:.0f}_{feature.split('_')[1]}")
                    else:
                        masked_true_pos.append(str(mapped_pos))
                except:
                    masked_true_pos.append(feature)          
        else:
            for feature in masked_cols:
                try:
                    pos_idx = feature[1:]
                    mapped_pos = site_translation_dict.get(pos_idx, feature)
                    if isinstance(mapped_pos, (int, float)):
                         masked_true_pos.append(f"{mapped_pos:.0f}_{feature.split('_')[1]}")
                    else:
                         masked_true_pos.append(str(mapped_pos))
                except:
                    masked_true_pos.append(feature)
        
        comparison_df = pd.DataFrame({
            'feature': np.array(cols)[mask],
            f'{name1}_states': feat1_vals[mask],
            f'{name2}_states': feat2_vals[mask],
            f'{name1}_shap': seq1['shap_values'][0][mask],
            f'{name2}_shap': seq2['shap_values'][0][mask],
            'shap_difference': shap_diff[mask], 
            'reference_position': masked_true_pos
        })
    else:
        comparison_df = pd.DataFrame({
            'feature': np.array(cols)[mask],
            f'{name1}_states': feat1_vals[mask],
            f'{name2}_states': feat2_vals[mask],
            f'{name1}_shap': seq1['shap_values'][0][mask],
            f'{name2}_shap': seq2['shap_values'][0][mask],
            'shap_difference': shap_diff[mask],
            'reference_position': np.array(cols)[mask] # Fallback
        })

    if encoding_method == 'one_hot':
        comparison_df[f'{name1}_states'] = comparison_df[f'{name1}_states'].astype(int)
        comparison_df[f'{name2}_states'] = comparison_df[f'{name2}_states'].astype(int)

    # Get Observed Amino Acids
    obs_aa_1 = []
    obs_aa_2 = []
    for feat in comparison_df['feature']:
        obs_aa_1.append(get_obs_aa(feat, seq1['aligned_seq_df'], encoding_method))
        obs_aa_2.append(get_obs_aa(feat, seq2['aligned_seq_df'], encoding_method))
        
    comparison_df[f'{name1}_obs_aa'] = obs_aa_1
    comparison_df[f'{name2}_obs_aa'] = obs_aa_2
    
    # Display values
    display_vals_1 = []
    display_vals_2 = []
    for idx, row in comparison_df.iterrows():
        display_vals_1.append(get_property_display_value(row['feature'], row[f'{name1}_states'], row[f'{name1}_obs_aa'], encoding_method))
        display_vals_2.append(get_property_display_value(row['feature'], row[f'{name2}_states'], row[f'{name2}_obs_aa'], encoding_method))
    
    comparison_df[f'{name1}_display'] = display_vals_1
    comparison_df[f'{name2}_display'] = display_vals_2
    comparison_df['abs_shap_diff'] = comparison_df['shap_difference'].abs()
    
    safe_n1 = name1.replace('|', '_')
    safe_n2 = name2.replace('|', '_')
    csv_name = f'{safe_n1}_vs_{safe_n2}_shap_data.csv'
    comparison_df.sort_values(by='abs_shap_diff', ascending=False).to_csv(f'{report_dir}/{csv_name}', index=False)
    
    comparison_df = comparison_df.sort_values(by='abs_shap_diff', ascending=False).head(n_positions)
    comparison_df = comparison_df.sort_values(by='shap_difference')

    if comparison_df.empty:
        return
    
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']
    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ["#30c898d9" if x < 0 else "#f4a365ec" for x in comparison_df['shap_difference']]
    
    # Plotting Axis selection
    if use_reference_sites_for_plot and isinstance(site_translation_dict, dict):
        bars = ax.barh(comparison_df['reference_position'], -comparison_df['shap_difference'], color=colors)
        y_label = f'Amino Acid Position ({ref_seq_name} Reference)'
    else:
        bars = ax.barh(comparison_df['feature'], -comparison_df['shap_difference'], color=colors)
        y_label = 'Feature'

    for i, bar in enumerate(bars):
        display_a = comparison_df[f'{name1}_display'].iloc[i]
        display_b = comparison_df[f'{name2}_display'].iloc[i]
        obs_a = comparison_df[f'{name1}_obs_aa'].iloc[i]
        obs_b = comparison_df[f'{name2}_obs_aa'].iloc[i]

        if encoding_method == 'aa_prop' and 'SCT' not in comparison_df['feature'].iloc[i]:
            try:
                val_str_a = f"{float(display_a):.2f}"
            except: val_str_a = str(display_a)
            try:
                val_str_b = f"{float(display_b):.2f}"
            except: val_str_b = str(display_b)
        else:
            val_str_a = f"{int(display_a):.0f}"
            val_str_b = f"{int(display_b):.0f}"
        
        label = f" {val_str_a} ({obs_a}) → {val_str_b} ({obs_b}) "
        ax.text(0, bar.get_y() + bar.get_height()/2, label,
                va='center', ha='right' if bar.get_width() < 0 else 'left',
                color='black', fontweight='bold', fontsize=10)

    ax.axvline(0, color='grey', linewidth=2, linestyle='--')
    ax.set_xlabel(r'Contribution to Difference in Predicted $λ_{max}$ (nm)', fontsize=13)
    ax.set_ylabel(y_label, fontsize=12)
    
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    pred_diff = seq1['prediction'] - seq2['prediction']
    
    ax.set_ylim(ymin, ymax + 0.80)
    ax.text((xmin + xmax)/2, ymax + 0.55,
            f"{name1} ({seq1['prediction']}) → {name2} ({seq2['prediction']})\n({pred_diff:.1f}nm Difference in Predicted " + r"$λ_{max}$)",
            ha='center', va='top', fontsize=13,
            bbox=dict(boxstyle='round,pad=0.3', fc="whitesmoke", ec='grey', lw=1, alpha=0.9))

    ax.tick_params(axis='y', labelrotation=45, labelsize=12)
    ax.tick_params(axis='x', labelsize=12)

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
    parser.add_argument("--mode", help="Analysis mode: 'comparison', 'single', 'both'.", type=str, default="both", choices=['comparison', 'single', 'both'])
    parser.add_argument("-v", "--model_version", help="Model version", type=str, default="vpod_1.3", choices=['vpod_1.3'])
    parser.add_argument("-m", "--model", help="Prediction model", type=str, default="whole-dataset")
    parser.add_argument("-e", "--encoding", help="Encoding method", type=str, default="aa_prop", choices=['one_hot', 'aa_prop'])
    parser.add_argument("--n_positions", help="Number of positions to show on SHAP explanation graphs.", type=int, default=10)
    parser.add_argument("--save_viz_as", help="File type", type=str, default="svg", choices=['svg','png','pdf'])
    parser.add_argument("--use_reference_sites", help="Use reference site numbering on plots", action="store_true")
    parser.add_argument("--n_jobs", help="Number of parallel processes.", type=int, default=-1)

    args = parser.parse_args()
    command_run = " ".join(sys.argv)

    generate_shap_explanation(
        args.input, args.output_dir, args.prediction_prefix, args.save_viz_as,
        args.model, args.encoding, args.model_version,
        cmd_line=command_run, mode=args.mode,
        n_positions=args.n_positions, use_reference_sites=args.use_reference_sites,
        n_jobs=args.n_jobs
    )