import os
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#from matplotlib.backends.backend_pdf import PdfPages
try:
    from deepBreaks.utils import load_obj
except:
    from optics_scripts.deepBreaks.utils import load_obj


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

    return mean_predictions, ci_lower, ci_upper, predictions_dict, median_predictions


def plot_predictions_with_CI(names, predictions, mean_preds, pdf_file):
#    if os.path.exists(pdf_file):
        # Customization
        plt.rcParams["figure.autolayout"] = True
        plt.rcParams["figure.figsize"] = [11.00, 5.00]
        #plt.subplots_adjust(left=0.1)  # More space to the left of the subplot

        # Customize colors 
        # Color palette with hex codes 
        colors=[wavelength_to_rgb(pred) for pred in mean_preds]
        color_specs = [matplotlib.colors.to_hex(color) for color in colors]  # Convert RGBA tuples to hex codes

        #Confidence interval calculation (you can adjust)
        confidence_level = 0.95

        #Box plots (group them together for each protein)
        boxprops = {'alpha': 0.6}  # Adjust alpha as needed
        bp_data = [predictions[seq] for seq in names] 
        bplot = plt.boxplot(bp_data, showfliers=False, positions=range(1, len(bp_data) + 1), vert = False, patch_artist=True, boxprops=boxprops, widths=0.5) 
        for i, seq in enumerate(names):
            y_jitter = np.random.uniform(-0.05, 0.05, size=len(predictions[seq]))  # Add slight y-jitter
            plt.scatter(predictions[seq], np.full_like(predictions[seq], i + 1)+ y_jitter, color='gray', s=20, alpha=0.5)

        for box, color in zip(bplot['boxes'], color_specs):
            box.set_facecolor(color)  # Set median line color to black

        for line in bplot['medians']:  
            line.set_color('black')  
            line.set_linewidth(1)  

        # Confidence intervals (calculated per protein sequence)
        for i, seq in enumerate(names):
            alpha = 1 - confidence_level
            lower_ci = np.percentile(predictions[seq], alpha / 2 * 100)
            upper_ci = np.percentile(predictions[seq], 100 - alpha / 2 * 100)
            plt.axhline(i + 1, lower_ci, upper_ci, color='b', linestyle='--')

        # Add protein names above boxes
        concat_names = []
        for name in names:
            if name.count('_') >= 3:
                temp = name.split('_')
                concat_names.append(f'{temp[0]}_{temp[1]}_{temp[2]}')
            else:
                concat_names.append(name)

        for i, seq_name in enumerate(concat_names):
            #if len(concat_names) > 4:
            #    y_pos = i + 1.50  # Center the text vertically
            #else:
            #    y_pos = i + 1.40  # Center the text vertically
            
            median_y = bplot['medians'][i].get_ydata()[0]  # Get the y-coordinate of the median line
            y_pos = median_y + 0.71 # Add an offset above the median
            x_pos = bplot['medians'][i].get_xdata()[1]  # Place slightly to the right of the median
            plt.text(x_pos, y_pos, seq_name, ha='center', va='top', color='black', fontsize=10, zorder = 3)  # Adjust horizontal alignment
        
        # Adjust plot limits to ensure annotations are visible
        #plt.ylim(top=plt.ylim()[1] + 0.2) # Increase the upper y-limit
        #plt.xlim(left=plt.xlim()[0] - 0.2) # Decrease the left x-limit by 0.2
        #getting median values to return for user
        medians = bplot['medians']
        median_values = [line.get_xydata()[1][0] for line in medians]  # Using list comprehension

        preds_handle = mpatches.Patch(facecolor='white', edgecolor='black', label='IQR')
        ci_handle = plt.Line2D([], [], color='black', label='95% Confidence Interval')

        #plt.title(f"λmax Predictions ({confidence_level * 100}% CI)")
        #seq_id = []
        #for i in range(len(names)):
            #seq_id.append(f"Seq{i}")
        plt.xlabel("Predicted λmax (nm)")
        plt.ylabel("Opsin Sequences")
        #plt.yticks(range(1, len(names) + 1), seq_id, rotation = 45)  # Set protein names as y-ticks
        plt.yticks([])
        plt.legend(handles=[preds_handle, ci_handle])
        #plt.legend(labels=['IQR', 'Confidence Interval']) 
        plt.grid(True, axis='x')  # Grid only on x-axis        
        plt.savefig(pdf_file, format = 'pdf', dpi = 300)
        #plt.show()
        plt.close() 

        return(median_values)



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