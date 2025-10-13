# optics_shap.py
import os
import sys
import subprocess
import argparse
import datetime
import tempfile
import pathlib

import joblib
import shap
import xgboost
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
from optics_scripts.utils import extract_fasta_entries



def generate_shap_explanation(
    input_file, pred_dir=None, output='unnamed_shap_comparison', save_as='svg',
                           model="whole-dataset", encoding_method='aa_prop',
                           model_version='vpod_1.3', refseq='bovine', reffile=None,
                           cmd_line=None
):
    """
    Generates and saves a SHAP plot to explain the difference in model predictions between two sequences.

    This function can generate two types of plots:
    1. A force plot for a single sequence, explaining its prediction.
    2. A comparative bar plot for two sequences, explaining the difference
       in their predictions.

    """
    dt_label = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    output = output.replace('.tsv', '').replace('.txt', '').replace('.csv', '')
    if pred_dir is None:
        report_dir = f'./prediction_outputs/optics_on_unamed_{dt_label}'
    else:
        report_dir = f'{pred_dir}/optics_shap_on_{output}_{dt_label}'
    os.makedirs(report_dir, exist_ok=True)  # exist_ok=True prevents errors if dir exists

    # If a command is passed and an output directory is specified, log the command for reproducibility.
    if cmd_line and report_dir:
        log_file_path = os.path.join(report_dir, 'shap_script_log.txt')
        with open(log_file_path, 'a') as f:
            f.write(f"Command executed:\n{cmd_line}\n\n")


    script_path = pathlib.Path(__file__).resolve()  # Get absolute path
    wrk_dir = str(script_path.parent).replace('\\', '/')
    if not os.path.isdir(f'{wrk_dir}/tmp'):
        os.makedirs(f'{wrk_dir}/tmp')
        
    data_dir = f"{wrk_dir}/data"
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
    # Extract sequences and their names from the input fasta file
    if input_file == None:
        raise Exception('Error: No file given')
    names,sequences = extract_fasta_entries(input_file)
    if len(names) > 2 or len(names) == 1:
        raise Exception('OPTICS SHAP Comparison only works for comparing two sequences at a time (for now).\n')
    alignment_data = model_datasets[model]
    metadata = model_metadata[model]
    model_path = model_directories[model]
    loaded_mod = load_obj(model_path)

    predictions = []
    encoded_sequences = []
    for name, sequence in zip(names, sequences):
        prediction, encoded_seq, encoded_refs = process_sequence(sequence, name, alignment_data, loaded_mod, wrk_dir)
        predictions.append(prediction)
        encoded_sequences.append(encoded_seq)

    # The base value is the average prediction over the entire training set
    explainer = shap.TreeExplainer(loaded_mod[-1], encoded_refs)
    base_value = explainer.expected_value
    print(f"\nSHAP Base Value (Average Predicted λmax): {float(base_value):.2f} nm")
    if cmd_line and report_dir:
        with open(log_file_path, 'a') as f:
            f.write(f"SHAP Base Value (Average Predicted lambda-max): {float(base_value):.2f} nm")

    seq1_shap = explainer.shap_values(encoded_sequences[0])
    seq2_shap = explainer.shap_values(encoded_sequences[1])
    shap_diff = seq1_shap[0] - seq2_shap[0]

    # Identify which features are different *after* transformation
    changed_features_mask = (encoded_sequences[0].values[0] != encoded_sequences[1].values[0])
    comparison_df = pd.DataFrame({
    'feature': np.array(encoded_sequences[0].columns)[changed_features_mask],
    f'{names[0]}_states': encoded_sequences[0].values[0][changed_features_mask],
    f'{names[1]}_states': encoded_sequences[1].values[0][changed_features_mask],
    f'{names[0]}_shap_effects': seq1_shap[0][changed_features_mask],
    f'{names[1]}_shap_effects': seq2_shap[0][changed_features_mask],
    'shap_difference': shap_diff[changed_features_mask]
    })

    if encoding_method == 'one_hot':
        comparison_df[f'{names[0]}_states']=comparison_df[f'{names[0]}_states'].astype(int)
        comparison_df[f'{names[0]}_states']=comparison_df[f'{names[0]}_states'].astype(int)
        comparison_df[f'{names[1]}_states']=comparison_df[f'{names[1]}_states'].astype(int)
        comparison_df[f'{names[1]}_states']=comparison_df[f'{names[1]}_states'].astype(int)

    comparison_df['abs_shap_diff'] = comparison_df['shap_difference'].abs()
    
    # Save the full DataFrame of SHAP values for all differing features to a TSV file.
    # This provides the detailed data for users interested in all feature contributions.
    full_shap_output_path = f'{report_dir}/{output}_all_shap_values.csv'
    # Sort by absolute SHAP difference for better interpretability in the output file.
    comparison_df_to_save = comparison_df.sort_values(by='abs_shap_diff', ascending=False)
    comparison_df_to_save.to_csv(full_shap_output_path, index=False)
    
    comparison_df_copy = comparison_df.copy()
    comparison_df_copy=comparison_df_copy.reset_index(drop=True).set_index('feature')

    print(f"\nSum of SHAP differences: {shap_diff.sum():.2f}nm")
    print(f"\nActual prediction difference: {abs(predictions[0] - predictions[1]):.2f}nm")
    print(f"\n\tOptics Prediction for {names[0]} = {predictions[0]:.2f}nm")
    print(f"\n\tOptics Prediction for {names[1]} = {predictions[1]:.2f}nm")
    if cmd_line and report_dir:
        with open(log_file_path, 'a') as f:
            f.write(f"\nSum of SHAP differences: {shap_diff.sum():.2f}nm")
            f.write(f"\nActual prediction difference: {abs(predictions[0] - predictions[1]):.2f}nm")
            f.write(f"\n\tOptics Prediction for {names[0]} = {predictions[0]:.2f}nm")
            f.write(f"\n\tOptics Prediction for {names[1]} = {predictions[1]:.2f}nm")

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']

    comparison_df = comparison_df.sort_values(by='abs_shap_diff', ascending=False).head(10)
    comparison_df = comparison_df.sort_values(by='shap_difference')

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ["#f5954bed" if x < 0 else "#30c8a4d9" for x in comparison_df['shap_difference']]
    bars = ax.barh(comparison_df['feature'], comparison_df['shap_difference'], color=colors)

    for i, bar in enumerate(bars):
        aa_a = comparison_df[f'{names[0]}_states'].iloc[i]
        aa_b = comparison_df[f'{names[1]}_states'].iloc[i]
        if encoding_method == 'one_hot':
            ax.text(0, bar.get_y() + bar.get_height()/2, f' {aa_a:.0f} → {aa_b:.0f} ',
                    va='center', ha='right' if bar.get_width() < 0 else 'left',
                    color='black', fontweight='bold', fontsize=10)
        else:
            ax.text(0, bar.get_y() + bar.get_height()/2, f' {aa_a:.2f} → {aa_b:.2f} ',
                    va='center', ha='right' if bar.get_width() < 0 else 'left',
                    color='black', fontweight='bold', fontsize=10)

    ax.axvline(0, color='grey', linewidth=2)
    ax.set_xlabel(f'Contribution to Difference in Predicted λmax (nm)', fontsize=12)
    ax.set_ylabel('Amino Acid Position/Feature', fontsize=12)
    plt.yticks(rotation=45)

    # Get the current y-axis limits before changing them.
    ymin, ymax = ax.get_ylim()

    # Increase the top y-axis limit to make space for the text.
    ax.set_ylim(ymin, ymax + 0.5)

    # Add the base value annotation at the top center (x=0).
    ax.text(0, ymax + 0.2, # x=0 for center, y is just above the original top bar
            f'{abs(predictions[0] - predictions[1]):.2f}nm Difference in Predicted λmax',
            ha='center', # Horizontal alignment
            va='top', # Vertical alignment
            fontsize=12,
            bbox=dict(boxstyle='round,pad=0.3', fc="#bdece3ac", ec='grey', lw=1, alpha=0.9))

    # Define the legend text based on the encoding method.
    if encoding_method == 'one_hot':
        legend_text = f"{names[0]} → {names[1]}\n1 = Present\n0 = Absent"
    else:
        legend_text = f"{names[0]} → {names[1]}"

    # Create an AnchoredText object for robust, automatic placement.
    # loc='lower right' places the box in the bottom right of the axes.
    # Other options include 'upper left', 'center', etc.
    # The 'prop' dictionary styles the text, and 'frameon=True' draws the box.
    anchored_text = AnchoredText(legend_text, loc='lower right',
                                 prop=dict(size=10, ha='right'), frameon=True,
                                 pad=0.5, borderpad=0.5)

    # Style the box to match your original aesthetics.
    anchored_text.patch.set_boxstyle("round,pad=0.3")
    anchored_text.patch.set_facecolor("whitesmoke")
    anchored_text.patch.set_edgecolor("grey")
    anchored_text.patch.set_linewidth(0.5)
    anchored_text.patch.set_alpha(0.8)

    # Add the automatically positioned text box to the plot.
    ax.add_artist(anchored_text)

    plt.tight_layout()
    plt.savefig(f'{report_dir}/{output}_viz_shap_comparison.{save_as}', dpi=350, bbox_inches='tight', format=save_as)
    plt.close()

    print(f"\nComparison explanation plot and SHAP values for all differing features saved to: {report_dir}")

def process_sequence(sequence, name, alignment_data, loaded_mod, wrk_dir):
    with tempfile.NamedTemporaryFile(mode="w", dir=f"{wrk_dir}/tmp", suffix=".fasta", delete=False) as temp_seq_file:
        temp_seq = temp_seq_file.name
        if '>' in sequence:
            temp_seq_file.write(sequence)
        else:
            temp_seq_file.write(f">placeholder_name\n{sequence}")

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

    return round(float(prediction[0]),1), loaded_mod[0].transform(new_seq_test), loaded_mod[0].transform(ref_copy)



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

    parser.add_argument("--refseq",
                            help="Reference sequence used for blastp analysis.",
                            type=str,
                            default="bovine",
                            choices=['bovine', 'squid', 'microbe', 'custom'])
    parser.add_argument("--custom_ref_file",
                            help="Path to a custom reference sequence file for BLASTp.",
                            type=str)  # No default, as it's optional

    parser.add_argument("--save_viz_as",
                                help="File type for bootstrap visualizations (SVG, PNG, or PDF).",
                                type=str,
                                default="svg",
                                choices=['svg','png','pdf'])

    args = parser.parse_args()

    # Capture the full command-line instruction for logging purposes.
    command_run = " ".join(sys.argv)

    generate_shap_explanation(args.input, args.output_dir, args.prediction_prefix, args.save_viz_as,
                           args.model, args.encoding, args.model_version,
                           args.refseq, args.custom_ref_file, cmd_line=command_run)
