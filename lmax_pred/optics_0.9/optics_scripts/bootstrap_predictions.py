import os
import numpy as np
import random
import matplotlib.pyplot as plt
from sklearn.utils import resample
from deepBreaks.utils import load_obj
from matplotlib.backends.backend_pdf import PdfPages

def calculate_ensemble_CI(model_folder, query, n_iterations=100):
    
    predictions_all = np.zeros((n_iterations, query.shape[0]))  # Account for multiple models

    # Load models
    i = 0
    for filename in os.listdir(model_folder):
        if filename.endswith('.pkl'):  # Assuming you saved models as .pkl
            model_path = os.path.join(model_folder, filename)
            model = load_obj(model_path)  # You'll need a 'load_model' function
            prediction = model.predict(query)
            prediction = round(float(prediction[0]))
            predictions_all[i] = prediction
            i+=1
    # Calculate ensemble confidence intervals
    ci_lower = np.percentile(predictions_all, 2.5, axis=1)
    ci_upper = np.percentile(predictions_all, 97.5, axis=1)
    mean_predictions = np.mean(predictions_all, axis=1)

    return mean_predictions, ci_lower, ci_upper


def plot_predictions_with_CI(name, mean_preds, ci_lower, ci_upper, pdf_file):
    with PdfPages(pdf_file) as pdf:  # Open the PDF file
        plt.figure(figsize=(8, 6))
        plt.plot(mean_preds, color='blue', label='Mean Prediction')
        plt.fill_between(np.arange(len(mean_preds)), ci_lower, ci_upper, color='blue', alpha=0.2, label='95% CI')
        plt.xlabel(f'{name}')
        plt.ylabel('λmax (nm)')
        plt.legend()
        pdf.savefig()  # Save the current figure to the PDF
        plt.close()  # Close the figure 
# ... (Rest of your code where you load data, call the functions) 
