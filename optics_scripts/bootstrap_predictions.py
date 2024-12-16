import os
import numpy as np
import random
from joblib import Parallel, delayed
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from deepBreaks.utils import load_obj
#try:
#    from deepBreaks.utils import load_obj
#except:
#    from optics_scripts.deepBreaks.utils import load_obj


def calculate_ensemble_CI(model_folder, query, name, predictions_dict):
    # Load models
    predictions_all = []
    for filename in os.listdir(model_folder):
        if filename.endswith('.pkl'):  # Assuming you saved models as .pkl
            model_path = os.path.join(model_folder, filename)
            model = load_obj(model_path)  # You'll need a 'load_model' function
            prediction = model.predict(query)
            prediction = round(float(prediction[0]),1)
            predictions_all.append(prediction)
    # Calculate ensemble confidence intervals
    predictions_all = np.array(predictions_all)
    predictions_dict.update({name: predictions_all})
    confidence_level = 0.95
    alpha = 1 - confidence_level
    ci_lower = np.percentile(predictions_all, alpha / 2 * 100)
    ci_upper = np.percentile(predictions_all, 100 - alpha / 2 * 100)
    mean_predictions = np.mean(predictions_all)
    median_predictions = np.median(predictions_all)
    std_dev = np.std(predictions_all)

    return mean_predictions, ci_lower, ci_upper, predictions_dict, median_predictions, std_dev

def plot_prediction_subsets_with_CI(names, predictions, mean_preds, pdf_file, visualize_bootstrap):
    # Customize colors 
    colors = [wavelength_to_rgb(pred) for pred in mean_preds]
    color_specs = [matplotlib.colors.to_hex(color) for color in colors]

    # Confidence interval calculation
    confidence_level = 0.95
    # Number of plots to generate
    num_plots = (len(names) + 4) // 5
    
    if visualize_bootstrap == 'True' or visualize_bootstrap == True:
        # Function to generate a plot for a subset of names
        for plot_idx in range(num_plots):
            plt.rcParams["figure.autolayout"] = True
            plt.rcParams["figure.figsize"] = [11.00, 5.00]
            plt.rcParams["figure.autolayout"] = True
            plt.rcParams["figure.figsize"] = [11.00, 5.00]
            
            start_idx = plot_idx * 5
            end_idx = min(start_idx + 5, len(names))
            subset_names = names[start_idx:end_idx]
            
            
            bp_data = [predictions[seq] for seq in subset_names]
            bplot = plt.boxplot(bp_data, showfliers=False, positions=range(1, len(bp_data) + 1), vert=False, patch_artist=True, boxprops={'alpha': 0.6}, widths=0.5)

            for i, seq in enumerate(subset_names):
                y_jitter = np.random.uniform(-0.05, 0.05, size=len(predictions[seq]))
                plt.scatter(predictions[seq], np.full_like(predictions[seq], i + 1) + y_jitter, color='gray', s=20, alpha=0.5)

            for box, color in zip(bplot['boxes'], color_specs[start_idx:start_idx + len(subset_names)]):
                box.set_facecolor(color)

            for line in bplot['medians']:
                line.set_color('black')
                line.set_linewidth(1)

            # Confidence intervals
            for i, seq in enumerate(subset_names):
                alpha = 1 - confidence_level
                lower_ci = np.percentile(predictions[seq], alpha / 2 * 100)
                upper_ci = np.percentile(predictions[seq], 100 - alpha / 2 * 100)
                plt.axhline(i + 1, lower_ci, upper_ci, color='b', linestyle='--')

            # Add protein names
            concat_names = []
            for name in subset_names:
                if name.count('_') >= 3:
                    temp = name.split('_')
                    concat_names.append(f'{temp[0]}_{temp[1]}_{temp[2]}')
                else:
                    concat_names.append(name)

            for i, seq_name in enumerate(concat_names):
                median_y = bplot['medians'][i].get_ydata()[0]
                y_pos = median_y + 0.71
                x_pos = bplot['medians'][i].get_xdata()[1]
                plt.text(x_pos, y_pos, seq_name, ha='center', va='top', color='black', fontsize=10, zorder=3)

            # plotting code for labels, legend, grid, etc.
            #medians = bplot['medians']
            #median_values = [line.get_xydata()[1][0] for line in medians]  # Using list comprehension

            preds_handle = mpatches.Patch(facecolor='white', edgecolor='black', label='IQR')
            ci_handle = plt.Line2D([], [], color='black', label='95% Confidence Interval')
            
            plt.xlabel("Predicted λmax (nm)")
            plt.ylabel("Opsin Sequences")
            #plt.yticks(range(1, len(names) + 1), seq_id, rotation = 45)  # Set protein names as y-ticks
            plt.yticks([])
            plt.legend(handles=[preds_handle, ci_handle])
            #plt.legend(labels=['IQR', 'Confidence Interval']) 
            plt.grid(True, axis='x')  # Grid only on x-axis        
            # Save each plot with a unique filename
            plt.savefig(f'{pdf_file}_part{plot_idx + 1}.pdf', format='pdf', dpi=300)
            plt.savefig(f'{pdf_file}_part{plot_idx + 1}.svg', format='svg')
            plt.close()
        
    #print('\nBootstrap plots done!\n')
    return(color_specs)

def wavelength_to_rgb(wavelength, gamma=0.8):
    ''' taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an 
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).

    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    Additionally alpha value set to 0.5 outside range
    '''
    wavelength = float(wavelength)
    if wavelength >= 380 and wavelength <= 750:
        A = 1.
    else:
        A=0.5
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R,G,B,A)