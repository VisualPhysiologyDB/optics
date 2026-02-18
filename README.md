**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  **VPOD_v1.2 DOI**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12213246.svg)](https://doi.org/10.5281/zenodo.12213246)

Opsin Phenotype Tool for Inference of Color Sensitivity (OPTICS) [v1.3]
=======================================================================

*Example Box Plot Output for Bootstrap Predictions of Opsin λmax by OPTICS*

Description
-----------

-   **OPTICS** is an open-source tool that predicts the Opsin Phenotype (λmax) from unaligned opsin amino-acid sequences.

-   **OPTICS** leverages machine learning models trained on the Visual Physiology Opsin Database (VPOD).

-   **OPTICS** allows for **structural mapping** of prediction features, translating machine learning insights directly onto 3D protein structures (PDB).

-   **OPTICS** can be downloaded and used as a command-line or GUI tool.

-   **OPTICS** is also available as an online tool [**here**](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest "null"), hosted on our [**Galaxy Project**](https://usegalaxy.org/ "null") server.

-   **Check out our pre-print** [**Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map**](https://doi.org/10.1101/2025.08.22.671864 "null") **to read more about it!**

Key Features
------------

-   **λmax Prediction**: Predicts the peak light absorption wavelength (λmax) for opsin proteins.

-   **Model Selection**: Choose from different pre-trained models for prediction (e.g., vertebrate, invertebrate, wildtype).

-   **Encoding Methods**: Select between one-hot encoding or amino-acid property encoding.

-   **BLAST Analysis**: Optionally perform BLASTp analysis to compare query sequences against reference datasets.

-   **Bootstrap Predictions**: Enable bootstrap predictions for enhanced accuracy assessment with confidence intervals.

-   **Prediction Explanation (SHAP)**: Explains the key features driving the λmax difference between sequences using SHAP values.

-   **Structure Mapping (NEW)**: Project SHAP importance values onto 3D PDB structures to create "importance heatmaps."

-   **Structure Annotation (NEW)**: Visualize custom annotations on 3D structures using automated PyMOL or ChimeraX scripting.

Table of Contents
-----------------

1.  [Installation](https://www.google.com/search?q=%23installation "null")

2.  [Data File Structure](https://www.google.com/search?q=%23data-file-structure "null")

3.  [Usage](https://www.google.com/search?q=%23usage "null")

    -   [Prediction: `optics_predictions.py`](https://www.google.com/search?q=%231-prediction-optics_predictionspy "null")

    -   [Explanation: `optics_shap.py`](https://www.google.com/search?q=%232-explanation-optics_shappy "null")

    -   [Structure Mapping: `optics_structure_map.py`](https://www.google.com/search?q=%233-structural-mapping-optics_structure_mappy "null")

    -   [Annotation: `optics_structure_annotations.py`](https://www.google.com/search?q=%234-structure-annotation-optics_structure_annotationspy "null")

    -   [GUI: `run_optics_gui.py`](https://www.google.com/search?q=%235-using-the-optics-gui "null")

4.  [Understanding Model Choice](https://www.google.com/search?q=%23understanding-the-%CE%BBmax-prediction-models "null")

5.  [License & Citation](https://www.google.com/search?q=%23license "null")

Installation
------------

1.  **Clone the repository:**

    ```
     git clone [https://github.com/VisualPhysiologyDB/optics.git](https://github.com/VisualPhysiologyDB/optics.git)

    ```

2.  **Install dependencies:** [Make sure you are working in the repository directory from here-after]

    A. Create a Conda environment for OPTICS (make sure you have [Conda](https://www.anaconda.com/ "null") installed)

    ```
    conda create --name optics_env python=3.11

    ```

    ### THEN

    ```
    conda activate optics_env

    ```

    B. Use the 'requirements.txt' file to download base package dependencies for OPTICS

    ```
    pip install -r requirements.txt

    ```

    C. **Download MAFFT and BLAST**

    IF working on MAC or LINUX device:

    -   Install *BLAST* and *MAFFT* directly from the *bioconda* channel

        ```
        conda install bioconda::blast bioconda::mafft

        ```

    IF working on WINDOWS device:

    -   Manually install the Windows compatible [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata "null") executable on your system PATH.

    -   You DO NOT need to download MAFFT; OPTICS includes a Windows-compatible version in the `optics_scripts/mafft` folder that it will try to use automatically.

Data File Structure
-------------------

OPTICS relies on a specific directory structure to locate models, alignment files, and cache data. When you clone the repository, the structure should generally look like this:

```
optics/
├── data/
│   ├── fasta/              # Alignment files for each model version (e.g., vpod_1.3)
│   ├── blast_dbs/          # BLAST databases for sequence identity checks
│   ├── aa_property_index/  # AA property values used for feature encoding
│   ├── importance_reports/ # Site translation dictionaries (Feature Name -> True Position)
│   ├── cached_predictions/ # Stores previous predictions (JSON) to speed up runtime
│   └── cached_structures/  # Stores downloaded PDB files (e.g., 1U19.pdb)
├── models/
│   ├── reg_models/         # Regression models (XGBoost/GradientBoosting) for point predictions
│   └── bs_models/          # Bootstrap model ensembles for confidence intervals
├── optics_scripts/         # Helper modules (utils, blast, bootstrap, maft wrappers, etc.)
└── prediction_outputs/     # Default output directory for all runs

```

*Note: The `cached_predictions` folder allows OPTICS to skip re-running heavy alignment/prediction steps for sequences it has seen before. You can clear this folder to force a fresh run.*

Usage
-----

**MAKE SURE YOU HAVE ALL DEPENDENCIES DOWNLOADED AND THAT YOU ARE IN THE FOLDER DIRECTORY FOR OPTICS (or have loaded it as a module) BEFORE RUNNING ANY SCRIPTS!**

### 1\. Prediction (`optics_predictions.py`)

The main script for generating λmax predictions.

```
Required Args:
  -i, --input: Either a single sequence string or a path to a FASTA file.

General Optional Args:
  -o, --output_dir: Desired directory to save output. Default: './prediction_outputs'
  -p, --prediction_prefix: Base filename for prediction outputs.
  -v, --model_version: Version of models (e.g., vpod_1.3).
  -m, --model: Prediction model to use. Default: whole-dataset
  -e, --encoding: Encoding method (one_hot, aa_prop). Default: aa_prop
  --tolerate_non_standard_aa: Allow sequences with non-standard AAs (X, B, Z) by stripping them. Default: True
  --tolerate_incomplete_seqs: Allow sequences outside 250-650aa range. Default: False
  --n_jobs: Number of parallel processes. Default: -1 (all CPUs)

BLASTp & Bootstrap Args (Optional):
  --blastp: Enable BLASTp analysis.
  --blastp_report: Filename for BLASTp report.
  --refseq: Reference sequence for BLAST (bovine, squid, microbe, custom).
  --bootstrap: Enable bootstrap predictions (provides confidence intervals).
  --bootstrap_num: Number of replicates (max 100).
  --visualize_bootstrap: Generate box-plots for bootstrap results.
  --full_spectrum_xaxis: Force x-axis to show 300-650nm range on plots.

```

**Example Command:**

```
python optics_predictions.py -i ./examples/optics_ex_short.txt -p ex_pred -m wildtype --blastp --bootstrap --visualize_bootstrap

```

### 2\. Explanation (`optics_shap.py`)

Generates SHAP (SHapley Additive exPlanations) plots to explain *why* a model predicted a specific value, or to attribute the difference between two sequences to specific amino acid sites.

```
Required Args:
  -i, --input: Path to FASTA file (must contain at least 2 sequences for comparison mode).

Optional Args:
  -o, --output_dir: Directory to save results.
  --mode: 'single' (explain one seq), 'comparison' (explain difference between two), or 'both'.
  --n_positions: Number of top features to show on the graph. Default: 10.
  --use_reference_sites: If set, plots use Reference Numbering (e.g., Bovine Rhodopsin sites) instead of alignment feature names.
  --save_viz_as: File format (svg, png, pdf).
  -m, --model: Prediction model to use.

```

**Example Command:**

```
python optics_shap.py -i ./examples/optics_ex_short.fasta -p shap_analysis --mode comparison --use_reference_sites

```

### 3\. Structural Mapping (`optics_structure_map.py`)

**NEW in v1.3!** This script takes the output CSV from the SHAP analysis and maps the importance values onto a 3D protein structure (PDB). It modifies the B-factor column of the PDB file, allowing you to visualize "importance" as a heat map (Blue=Low, Red=High importance).

```
Required Args:
  -s, --shap_csv: Path to the SHAP analysis CSV file generated by optics_shap.py.

Optional Args:
  -p, --pdb_file: Path to a local PDB file OR a 4-letter PDB ID (e.g., 1U19). Default: 1U19 (Bovine Rhodopsin).
  -o, --output_dir: Output directory.
  --chain: Chain ID to map to. Default: A.
  --use_query_position: Check this if your CSV uses target sequence numbering rather than reference numbering.
  --map_bovine_also: If using a custom PDB, this flag forces a second output mapped to Bovine Rhodopsin (1U19) for comparison.

```

**Example Command:**

```
python optics_structure_map.py -s ./examples/shap_output/my_seq_shap_analysis.csv -p 1U19 --map_bovine_also

```

*Output: Generates a `.pdb` file with importance scores in the B-factor column and a `.pml` script to automatically visualize it in PyMOL.*

### 4\. Structure Annotation (`optics_structure_annotations.py`)

**NEW in v1.3!** A general-purpose tool to visualize arbitrary annotations (e.g., mutation sites, binding pockets) on a structure. It takes a simple CSV and creates a runnable visualization script.

```
Required Args:
  -a, --annotation_file: CSV/TSV file with columns: 'position' (required), 'color' (optional), 'style' (optional), 'label' (optional).

Optional Args:
  -p, --pdb: PDB ID or path. Default: 1U19.
  --software: Target visualization software ('pymol' or 'chimerax'). Default: chimerax.
  --chain: Chain identifier.

```

**Example Command:**

```
python optics_structure_annotations.py -a ./examples/mutations.csv -p 1U19 --software pymol

```

### 5\. Using the OPTICS GUI

OPTICS includes a graphical interface for users who prefer not to use the command line.

**To run the GUI:**

```
python run_optics_gui.py

```

<img src="https://github.com/VisualPhysiologyDB/optics/blob/main/data/logo/optics_gui_ex.png?raw=true" alt="ex optics gui" style="width:50%; height:50%;">

The GUI provides tabs/buttons for all four major pipelines:

1.  **Standard Predictions**: Run the main λmax prediction workflow.

2.  **SHAP Interpretation**: Run feature attribution analysis.

3.  **Structure Mapping**: Map SHAP values to PDB files.

4.  **Structure Annotations**: Visualize custom data on structures.

Understanding the λmax Prediction Models
----------------------------------------

The `--model` flag allows you to select a specific pre-trained model. Each is named after the data-subset it was trained on:

### **Base Model Datasets**

-   **whole-dataset**: Trained on the entire VPOD dataset. **Recommended**.

-   **wildtype**: Trained exclusively on wild-type sequences.

-   **vertebrate** / **invertebrate**: Taxonomic subsets.

-   **wildtype-vert**: Wild-type vertebrate sequences only.

### **The `-mnm` Suffix**

Models ending in **-mnm** (e.g., `wildtype-mnm`) are trained on augmented datasets.

-   **Standard models**: Trained *exclusively* on heterologous expression data (in-vitro).

-   **-mnm models**: Trained on heterologous data *plus* data inferred via our **"Mine-n-Match"** procedure (in-vivo correlations). See [Frazer et al. 2025](https://doi.org/10.1101/2025.08.22.671864 "null") for details.

License
-------

All data and code is covered under a GNU General Public License (GPL)(Version 3).

Citation
--------

-   **Code/Repository:** 10.5281/zenodo.10667840

-   **OPTICS Publication (Methodology & Tools):** Seth A. Frazer, Todd H. Oakley. Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map. bioRxiv, 2025.08.22.671864. https://doi.org/10.1101/2025.08.22.671864

-   **VPOD Publication (Database & Training Data):** Seth A. Frazer, Mahdi Baghbanzadeh, Ali Rahnavard, Keith A. Crandall, & Todd H Oakley. Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD). GigaScience, 2024.09.01. https://doi.org/10.1093/gigascience/giae073

Contact
-------

**Todd H. Oakley** - [ORCID](https://orcid.org/0000-0002-4478-915X "null") - oakley@ucsb.edu

**Seth A. Frazer** - [ORCID](https://orcid.org/0000-0002-3800-212X "null") - sethfrazer@ucsb.edu