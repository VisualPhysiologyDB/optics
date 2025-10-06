import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from deepBreaks.utils import load_obj

def calculate_ensemble_CI(original_prediction, loaded_bs_models, query, name, bootstrap_num=100, model_folder=None,):
    """
    Calculates ensemble predictions and confidence intervals from pre-loaded bootstrap models.

    Args:
        loaded_bs_models (list): A list of pre-loaded model objects.
        query (DataFrame): The input data for prediction.
        name (str): The name of the sequence being processed.
        predictions_dict (Manager.dict): A shared dictionary to store raw bootstrap predictions.

    Returns:
        tuple: Contains mean, CI lower, CI upper, updated predictions_dict, median,
               standard deviation, and a list of all predictions.
    """
    # Use the pre-loaded or real-time load boot-strapped models to make predictions
    
    # List for storing all the bootstrap predictions
    predictions_all = []
    # Append original prediction to list of predictions
    predictions_all.append(round(float(original_prediction[0]),1))

    if len(loaded_bs_models) == bootstrap_num:
        for model in loaded_bs_models:
            prediction = model.predict(query)
            prediction = round(float(prediction[0]),1)
            predictions_all.append(prediction)
    else:
        i = 0
        for filename in os.listdir(model_folder):
            if filename.endswith('.pkl') and i < bootstrap_num:  # Assuming you saved models as .pkl
                model_path = os.path.join(model_folder, filename)
                model = load_obj(model_path)  # You'll need a 'load_model' function
                prediction = model.predict(query)
                prediction = round(float(prediction[0]),1)
                predictions_all.append(prediction)
                
                i+= 1


    # Calculate ensemble confidence intervals
    predictions_all = np.array(predictions_all)
    confidence_level = 0.95
    alpha = 1 - confidence_level
    ci_lower = np.percentile(predictions_all, alpha / 2 * bootstrap_num)
    ci_upper = np.percentile(predictions_all, bootstrap_num - alpha / 2 * bootstrap_num)
    mean_predictions = np.mean(predictions_all)
    median_predictions = np.median(predictions_all)
    std_dev = np.std(predictions_all)

    return mean_predictions, ci_lower, ci_upper, median_predictions, std_dev, predictions_all

def plot_prediction_subsets_with_CI(names, predictions, mean_preds, pdf_file, visualize_bootstrap, save_as='svg', full_spectrum_xaxis = False):
    
    """Generates and saves boxplots of bootstrapped predictions for sequences.

    This function creates a series of plots, with each displaying up to five
    sequences. For each sequence, it visualizes the distribution of
    bootstrapped predictions as a boxplot overlaid with a scatter plot of
    the individual bootstrap values. A 95% confidence interval is also
    indicated. Plots are saved as separate image files.

    Parameters
    ----------
    names : list of str
        A list of sequence identifiers. The order determines the plotting order.
    predictions : dict
        A dictionary mapping sequence names from `names` to their corresponding
        list or array of bootstrapped prediction values (e.g., λmax).
    mean_preds : list or array-like
        A list of mean prediction values for each sequence, used to determine
        the color of the boxplots. The order must match the `names` list.
    pdf_file : str
        The base path and filename for the output plots, without the file
        extension. E.g., 'results/my_analysis'.
    visualize_bootstrap : bool or str
        If `True` or the string 'True', the function will generate and save
        the plots. Otherwise, it will do nothing.
    save_as : str, optional
        The file format for the saved plots (e.g., 'svg', 'png').
        Defaults to 'svg'.
    full_spectrum_xaxis : bool, optional
        If `True`, the x-axis of the plots is fixed to the range [300, 650] nm.
        Defaults to `False`, which allows matplotlib to auto-scale the axis.

    Notes
    -----
    - The function saves output plots as a series of files named using the
      `pdf_file` base, such as `{pdf_file}_part1.{save_as}`,
      `{pdf_file}_part2.{save_as}`, etc.
    - The color of each boxplot is determined by mapping the `mean_preds`
      value to a color via the `wavelength_to_rgb` helper function.
    """
    
    # Customize colors 
    colors = [wavelength_to_rgb(pred) for pred in mean_preds]
    color_specs = [matplotlib.colors.to_hex(color) for color in colors]

    # Number of plots to generate
    num_plots = (len(names) + 4) // 5
    
    if visualize_bootstrap == 'True' or visualize_bootstrap == True:
        plt.rcParams["figure.autolayout"] = True  
        
              
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Century Gothic', 'Avenir', 'Helvetica', 'DejaVu Sans', 'Arial']

        # Function to generate a plot for a subset of names
        for plot_idx in range(num_plots):

            # Explicitly create a Figure and Axes object for each plot
            fig, ax = plt.subplots(figsize=(11.00, 5.00))

            start_idx = plot_idx * 5
            end_idx = min(start_idx + 5, len(names))
            subset_names = names[start_idx:end_idx]

            bp_data = [predictions[seq] for seq in subset_names]
            # Use the 'ax' object for all plotting functions
            bplot = ax.boxplot(bp_data, showfliers=False, positions=range(1, len(bp_data) + 1), vert=False, patch_artist=True, boxprops={'alpha': 0.6}, widths=0.5)

            for i, seq in enumerate(subset_names):
                y_jitter = np.random.uniform(-0.05, 0.05, size=len(predictions[seq]))
                ax.scatter(predictions[seq], np.full_like(predictions[seq], i + 1) + y_jitter, color='gray', s=20, alpha=0.5)

            for box, color in zip(bplot['boxes'], color_specs[start_idx:start_idx + len(subset_names)]):
                box.set_facecolor(color)

            for line in bplot['medians']:
                line.set_color('black')
                line.set_linewidth(1)

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
                if len(concat_names) >= 3:
                    y_pos = median_y + 0.70
                elif len(concat_names) == 2:
                    y_pos = median_y + 0.60
                else:
                    y_pos = median_y + 0.57

                x_pos = bplot['medians'][i].get_xdata()[1]
                ax.text(x_pos, y_pos, seq_name, ha='center', va='top', color='black', fontsize=12, zorder=3)

            preds_handle = mpatches.Patch(facecolor='white', edgecolor='black', label='IQR')

            if full_spectrum_xaxis:
                ax.set_xlim(300, 650)
            ax.set_xlabel("Predicted λmax (nm)", fontsize=15)
            ax.set_ylabel("Opsin Sequences", fontsize=15)
            ax.tick_params(axis='x', labelsize=15)
            ax.set_yticks([])
            ax.legend(handles=[preds_handle, plt.Line2D([], [], color='black', linestyle='-', label='95% Confidence Interval')])
            ax.grid(True, axis='x')

            # Save the specific figure object
            fig.savefig(f'{pdf_file}_part{plot_idx + 1}.{save_as}', format=save_as, dpi=400)

            # Explicitly close the figure to release memory
            plt.close(fig)
        
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