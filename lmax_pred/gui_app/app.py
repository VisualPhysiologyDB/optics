#app.py

from flask import Flask, render_template, request, jsonify
from prediction_functions import process_sequence, process_sequences_from_file  # Adjust the import name if needed

app = Flask(__name__)

# Example dictionary to map model names to data paths (adjust for your setup)
model_datasets = {
    "Whole Dataset Model": "./data/whole_dataset.fasta",
    "Vertebrate Model": "./data/vertebrate_opsins.fasta",
    "Rod Model": "./data/rod_opsins.fasta",

}  
import logging 
app.logger.setLevel(logging.DEBUG)  # Set to DEBUG for detailed lo

@app.route('/')
def index():
    models = list(model_datasets.keys())  # Get model names
    return render_template('index.html', models=models)

@app.route('/predict', methods=['POST'])
def predict():
    try:
        sequence = request.form['sequence']
        selected_model = request.form['model']
        print("Received sequence:", sequence)  # Check if data is received
        print("Selected model:", selected_model) 
        
    except Exception as e:  # Catch any exceptions
        print(f"An error occurred: {e}")
        return jsonify({'error': 'An error occurred during prediction'}), 500 

    try:
        prediction = process_sequence(sequence, selected_model) 
        print("Prediction:", prediction)
        return jsonify({'prediction': prediction})
    except Exception as e:  # Catch any exceptions
        print(f"An error occurred: {e}")
        return jsonify({'error': 'An error occurred during prediction'}), 500 

@app.route('/predict_from_file', methods=['POST'])
def predict_from_file():
    app.logger.debug("File upload request received")
    if 'file' not in request.files:
        print(f"An error occurred while retireving file data - catch 1")
        return jsonify({'error': 'No file uploaded'}), 500 
    #if 'model2' not in request.form: 
        #print(f"An error occurred while retireving model form data")
        #return jsonify({'error': 'No Model Selected'}), 500

    
    try:
        file = request.files['file'] 
    except Exception as e:  # Catch any exceptions
        print(f"An error occurred while retireving file data: {e}")
        return jsonify({'error': 'An error occurred while retireving file data'}), 500 
    
    try:
        selected_model = request.form['model']
        print(selected_model)
    # Consider adding file validation (size, type) here
    except Exception as e:  # Catch any exceptions
        print(f"An error occurred while retireving model form data: {e}")
        return jsonify({'error': 'An error occurred while retireving model form data'}), 500 

    try:
        names, predictions = process_sequences_from_file(file, selected_model)  # Call your file processing function
        for i in range(len(names)-1):
            print(f'{names[i]} : {predictions[i]}')
        return jsonify({
            'names' : names, # Return a list of names corresponding to predictions 
            'predictions': predictions})  # Return a list of predictions 
    except Exception as e:  # Catch any exceptions
        print(f"An error occurred during prediction: {e}")
        return jsonify({'error': 'An error occurred during prediction'}), 500 

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)