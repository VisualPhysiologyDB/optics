# optics_shap.py
import shap
import xgboost
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from deepBreaks.utils import load_obj
from deepBreaks.preprocessing import read_data
import os

def generate_shap_explanation(
    model,
    processed_seq_A,
    processed_seq_B=None,
    background_data=None,
    output_path="shap_explanation.png",
    sequence_name_A="Sequence A",
    sequence_name_B="Sequence B"
):
    """
    Generates and saves a SHAP plot to explain model predictions.

    This function can generate two types of plots:
    1. A force plot for a single sequence, explaining its prediction.
    2. A comparative bar plot for two sequences, explaining the difference
       in their predictions.

    Args:
        model (ML model): The trained and loaded ML model object.
        processed_seq_A (pd.DataFrame): The preprocessed feature vector for the first sequence.
        processed_seq_B (pd.DataFrame, optional): The preprocessed feature vector for the second sequence.
                                                   If provided, a comparison plot is generated. Defaults to None.
        background_data (pd.DataFrame, optional): A subset of the training data used to initialize the
                                                  SHAP explainer for a more robust explanation. Defaults to None.
        output_path (str, optional): Path to save the output plot image. Defaults to "shap_explanation.png".
        sequence_name_A (str, optional): Display name for the first sequence. Defaults to "Sequence A".
        sequence_name_B (str, optional): Display name for the second sequence. Defaults to "Sequence B".
    """
    print("Initializing SHAP explainer...")
    # Using shap.TreeExplainer is highly efficient for tree-based models like XGBoost.
    # Providing background data gives a better reference for expected values.
    explainer = shap.TreeExplainer(model, background_data)
    
    # --- Case 1: Explanation for a single sequence ---
    if processed_seq_B is None:
        print(f"Generating SHAP force plot for {sequence_name_A}...")
        shap_values = explainer.shap_values(processed_seq_A)
        prediction = model.predict(processed_seq_A)[0]

        # Create the force plot
        plt.figure()
        shap.force_plot(
            explainer.expected_value,
            shap_values,
            processed_seq_A.iloc[0],
            matplotlib=True,
            show=False
        )
        plt.title(f"SHAP Explanation for {sequence_name_A} (Predicted λmax: {prediction:.2f} nm)", fontsize=12)
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Single sequence explanation plot saved to: {output_path}")

    # --- Case 2: Explanation for the difference between two sequences ---
    else:
        print(f"Generating SHAP comparison plot for {sequence_name_A} vs {sequence_name_B}...")
        
        # Calculate SHAP values for both sequences
        shap_values_A = explainer.shap_values(processed_seq_A)
        shap_values_B = explainer.shap_values(processed_seq_B)
        
        # The difference in SHAP values explains the difference in predictions
        shap_diff = shap_values_B[0] - shap_values_A[0]
        
        # Get model predictions
        pred_A = model.predict(processed_seq_A)[0]
        pred_B = model.predict(processed_seq_B)[0]
        pred_diff = pred_B - pred_A

        # Identify which features (amino acid positions) are different
        # This assumes one-hot encoding where a difference means a non-zero change
        feature_names = processed_seq_A.columns
        changed_features_mask = (processed_seq_A.values[0] != processed_seq_B.values[0])
        
        comparison_df = pd.DataFrame({
            'feature': feature_names[changed_features_mask],
            'shap_difference': shap_diff[changed_features_mask]
        })
        
        # Sort by the magnitude of the contribution to the difference
        comparison_df['abs_shap_diff'] = comparison_df['shap_difference'].abs()
        comparison_df = comparison_df.sort_values(by='abs_shap_diff', ascending=False).head(20) # Top 20 drivers
        comparison_df = comparison_df.sort_values(by='shap_difference') # Sort for plotting

        # Create the comparison bar plot
        fig, ax = plt.subplots(figsize=(12, 8))
        colors = ['#ff4d4d' if x < 0 else '#4d4dff' for x in comparison_df['shap_difference']]
        ax.barh(comparison_df['feature'], comparison_df['shap_difference'], color=colors)
        
        ax.axvline(0, color='grey', linewidth=0.8)
        ax.set_xlabel(f'Contribution to λmax Difference ({sequence_name_B} vs {sequence_name_A})', fontsize=12)
        ax.set_ylabel('Amino Acid Position/Feature', fontsize=12)
        ax.set_title(f'Key Drivers of the {pred_diff:.2f} nm Shift in Predicted λmax', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()
        print(f"Comparison explanation plot saved to: {output_path}")