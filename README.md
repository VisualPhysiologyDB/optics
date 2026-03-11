**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  **VPOD_v1.2 DOI**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12213246.svg)](https://doi.org/10.5281/zenodo.12213246)

# Opsin Phenotype Tool for Inference of Color Sensitivity (OPTICS) [v1.3] - 
*_Note_* - A simpler version of an intro to OPTICS is also available on our organization github.io page -> [here](https://visualphysiologydb.github.io/optics.html)

![](https://github.com/VisualPhysiologyDB/optics/blob/main/data/logo/optics_bs_fig_ex.svg?raw=true)

  _Example Box Plot Output for Bootstrap Predictions of Opsin λmax by OPTICS_

Description
-----------

-   **OPTICS** is an open-source tool that uses machine learning (ML) models to predict Opsin Phenotype (λmax) from unaligned opsin amino-acid sequences.

-   **OPTICS** leverages machine learning models trained on opsin genotype-phenotype data from the [**Visual Physiology Opsin Database (VPOD)**](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db).

-   **OPTICS** allows for **structural mapping** of sequence features important to model prediction (using [SHAP](https://shap.readthedocs.io/en/latest/)), translating machine learning insights directly onto 3D protein structures (PDB).

-   **OPTICS** can be downloaded and used as a command-line or GUI tool.

-   **OPTICS** is also available as an online tool [**here**](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest "null"), hosted on our [**Galaxy Project**](https://usegalaxy.org/ "null") server.

-   **Check out our pre-print** [**Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map**](https://doi.org/10.1101/2025.08.22.671864 "null") **to read more about it!**

Key Features
------------

-   **λmax Prediction**: Predicts the peak light absorption wavelength (λmax) for opsin proteins.

-   **Model Selection**: Choose from different pre-trained models for prediction.

-   **BLAST Analysis**: Optionally perform BLASTp analysis to compare query sequences against reference datasets.

-   **Bootstrap Predictions**: Enable bootstrap predictions for enhanced accuracy assessment with confidence intervals.

-   **Prediction Explanation (SHAP)**: Explains the key sequence features driving model predictions of λmax. This feature also allows for all-to-all pairwise comparisons of the features driving differences in predicted λmax between sequences using SHAP values.

-   **Structure Mapping**: Project SHAP importance values onto 3D PDB structures to create "importance heatmaps."

-   **Custom Structure Annotation**: Visualize custom annotations on 3D structures using automated PyMOL or ChimeraX scripting.

Table of Contents
-----------------

1.  [Installation](#installation)

2.  [Data File Structure](#data-file-structure)

3.  [Usage](#usage)

    -   [Prediction: `optics_predictions.py`](#1-λmax-prediction-optics_predictionspy)

    -   [SHAP Explanation: `optics_shap.py`](#2-explaining-model-predictions-with-shap-optics_shappy)

    -   [SHAP Structure Mapping: `optics_structure_map.py`](#3-mapping-shap-importance-to-3d-structure-optics_structure_mappy)

    -   [Custom Structure Annotation: `optics_structure_annotations.py`](#4-generate-custom-structure-annotations-optics_structure_annotationspy)

    -   [GUI: `run_optics_gui.py`](#5-using-the-optics-gui)

4.  [Understanding Model Choice](#understanding-the-λmax-prediction-models)

5. [License](#license)

6. [Citation](#citation)

7. [Contact](#contact)

8. [Additional Resources](#additional-notesresources)

Installation
------------

1.  **Clone the repository:**

    ```
     git clone https://github.com/VisualPhysiologyDB/optics.git
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

    -   DO NOT need to download MAFFT; OPTICS includes a Windows-compatible version in the `optics_scripts/mafft` folder that it will try to use automatically.
      - You can download your own version of MAFFT but it must be executable on your system path. 

---

Data File Structure
-------------------

OPTICS relies on a specific directory structure to locate models, alignment files, and cache data. When you clone the repository, the structure should generally look like this:

```
optics/
├── data/
│   ├── fasta/              # Alignment files for each model version (e.g., vpod_1.3)
│   ├── blast_dbs/          # BLAST databases for sequence identity checks
│   ├── aa_property_index/  # AA property values used for feature encoding
│   ├── importance_reports/ # Feature importance data & site translation information (Feature Name -> True Position)
│   ├── cached_structures/  # Stores downloaded PDB files (e.g., 1U19.pdb)
│   ├── cached_predictions/ # Stores previous predictions (JSON) to speed up runtime
|   └── cached_blastp_analysis/ # Stores data from previous runs of BLASTp (JSON) to speed up runtime
├── models/
│   ├── reg_models/         # Regression models (XGBoost/GradientBoosting) for point predictions
│   └── bs_models/          # Bootstrap model ensembles for confidence intervals
├── optics_scripts/         # Helper modules (utils, blast, bootstrap, maft wrappers, etc.)
├── deepBreaks/             # A key component of the OPTICS pipeline, this folder must stay here
└── prediction_outputs/     # Default output directory for all runs

```

*Note: The `cached_predictions` folder allows OPTICS to skip re-running heavy alignment/prediction steps for sequences it has seen before. You can clear this folder to force a fresh run.*

Usage
-----

**MAKE SURE YOU HAVE ALL DEPENDENCIES DOWNLOADED AND THAT YOU ARE IN THE FOLDER DIRECTORY FOR OPTICS (or have loaded it as a module) BEFORE RUNNING ANY SCRIPTS!**

### 1\. λmax Prediction (`optics_predictions.py`)

The main script for generating λmax predictions.

```
Required Args:

  -i, --input: Either a single sequence or a path to a FASTA file.

General Optional Args:

  -o, --output_dir: Desired directory to save output folder/files (optional). Default: './prediction_outputs'

  -p, --prediction_prefix: Base filename for prediction outputs. Default: 'unnamed'

  -v, --model_version: Version of models to use (optional). Based on the version of VPOD used to train models. Options/Default: vpod_1.3 (More version coming later)

  -m, --model: Prediction model to use. Options: whole-dataset, wildtype, vertebrate, invertebrate, wildtype-vert, type-one, whole-dataset-mnm, wildtype-mnm, vertebrate-mnm, invertebrate-mnm, wildtype-vert-mnm. **Default: whole-dataset** 

  -e, --encoding: Encoding method to use (optional). Options: one_hot, aa_prop. Default: aa_prop

  --tolerate_non_standard_aa: Allows OPTICS to run predictions on sequences with 'non-standard' amino-acids (e.g. - 'X','O','B', etc...)(optional). Default: True

  --tolerate_incomplete_seqs: Allows OPTICS to run predictions on sequences outside the predefined limits of 250-650 amino-acids. (optional) Default: False 
                              NOTE - if you enable this option, then you may have predictions on incomplete sequences, which should be treated as less accurate.

  --n_jobs: Number of parallel processes to run (optional). -1 is the default, utilizing all avaiable processors.


BLASTp Analysis Args (optional):

  --blastp: Enable BLASTp analysis.

  --blastp_report: Filename for BLASTp report. Default: blastp_report.txt

  --refseq: Reference sequence used for blastp analysis. Options: bovine, squid, microbe, custom. Default: bovine

  --custom_ref_file: Path to a custom reference sequence file for BLASTp.  Required if --refseq custom is selected.

Bootstrap Analysis Args (optional):

  --bootstrap: Enable bootstrap predictions.

  --visualize_bootstrap: Enable visualization of bootstrap predictions.

  --bootstrap_num: Number of bootstrap models to load for prediction replicates. Default // Maximum: 100

  --bootstrap_viz_file: Filename prefix for bootstrap visualization. Default: bootstrap_viz

  --save_viz_as: File type for bootstrap visualizations. Options: svg, png, or pdf Default: svg
  
  --full_spectrum_xaxis: Enables visualization of predictions on a full spectrum x-axis (300-650nm). Otherwise, x-axis is scaled with predictions.

```

**Example Command:**

```
  python optics_predictions.py -i ./examples/optics_ex_short.txt -o ./examples -p ex_predictions -m whole-dataset -e aa_prop --blastp --blastp_report blastp_report_ex --refseq squid --bootstrap --visualize_bootstrap --bootstrap_viz_file bootstrap_viz --save_viz_as svg
```

### Input

- **Unaligned** FASTA file containing opsin amino-acid sequences.
- Example FASTA Entry:
  ```
    >NP_001014890.1_rhodopsin_Bos_taurus
    MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRT 
    PLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVC 
    KPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVV 
    HFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQG 
    SDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA   
  ```
  
### Output

- Predictions (TSV, Excel): λmax values, BLASTp information, hex-codes (colors) corresponding to predicted λmax.
- BLAST Results (TXT, optional): Comparison of query sequences to reference datasets.
- Bootstrap Graphs (PDF, optional): Visualization of bootstrap prediction results.
- Job Log (TXT): Log file containing input command to OPTICS, including encoding method and model used.
- iTol & FigTree Annotation Files (TXT): Annotations for visualizing the λmax of opsins, with the hex-codes (colors) corresponding to predicted λmax.

  _**Note** - All outputs are written into subfolders generated based on your 'prediction-prefix' under your specified output directory, and are marked by time and date._

### 2\. Explaining Model Predictions with SHAP (`optics_shap.py`)

For users interested in the "nitty-gritty" of _why_ sequences have different predicted λmax values, we provide a specialized script that uses *SHAP* ([SHapley Additive exPlanations](https://shap.readthedocs.io/en/latest/)). 

This tool generates detailed plots and reports that attribute the difference in prediction to specific features (i.e., amino acid sites and their properties).

![](https://github.com/VisualPhysiologyDB/optics/blob/main/examples/optics_shap_on_short_ex_test_aa_prop_2026-01-09_02-28-40/Bombus_impatiens_424_individual_shap.svg?raw=true)

_Example SHAP plot for explaining individual predictions of opsin λmax by OPTICS_

![](https://github.com/VisualPhysiologyDB/optics/blob/main/examples/optics_shap_on_short_ex_test_aa_prop_2026-01-09_02-28-40/Bombus_impatiens_424_vs_Bombus_impatiens_347_viz.svg?raw=true)

_Example SHAP comparison plot for explaining pair-wise differences in predictions of opsin λmax by OPTICS_

This script requires a **FASTA file** 
- File must contain at least **two or more sequences** if you are running a SHAP comparison.
- Only a single sequence is needed for an individual SHAP explanation

```
Required Args:
  -i, --input: Path to FASTA file (must contain at least 2 sequences for comparison mode).

Optional Args:
  -o, --output_dir: Directory to save the SHAP analysis output folder.

  -p, --prediction_prefix: Base filename for the SHAP plot and data files.

  --mode: Analysis mode - 'single' (generate SHAP explanation for any number of individual sequences), 'comparison' (generate SHAP explanation for pairwise predictionn difference between any number of sequences), or 'both'. Default: 'both'

  -m, --model: Prediction model to use.

  -v, --model_version: Version of models to use (optional). Based on the version of VPOD used to train models. Options/Default: vpod_1.3 (More version coming later)

  -e, --encoding: Encoding method to use.

  --n_positions: Number of positions to show on SHAP explanation graphs. Default: 10 (to limit noisiness)

  --save_viz_as: File type for the SHAP visualization (svg, png, or pdf). Default: svg

  --use_reference_sites : Enable to use reference site numbering (i.e. - Bovine or Squid Rhodopsin), instead of feature names.

  --n_jobs: Number of parallel processes to run (optional). -1 is the default, utilizing all avaiable processors., 

```

**Example Command:**

```
python optics_shap.py -i ./examples/optics_ex_short.fasta -o ./examples -p short_ex_test_aa_prop --mode both --use_reference_sites
```

### Output

- Predictions (TSV): Single (non-bootstrapped) λmax values 
- SHAP Explanation Data (CSV): SHAP data for individual sequence explanations and/or pairwise SHAP comparisons for all sequences.
- SHAP Graphs (SVG): Visualizations for individual sequence SHAP explanations and/or pairwise SHAP comparisons for all sequences.
- Difference Matrix (CSV & Excel): An all-to-all λmax difference matrix for query sequences. The Excel version has a built in 'heat-map'.
- Job Log (TXT): Log file containing input command to OPTICS, including encoding method and model used.

_**WARNING** - Be cautious if you choose the 'comparison' or 'both' mode for SHAP with many sequences. Too many sequence can end up generating hundreds-to-thousands of comparison files (raw data and visualizations)._

### 3\. Mapping SHAP Importance to 3D Structure (`optics_structure_map.py`)

<img src="https://github.com/VisualPhysiologyDB/optics/blob/main/examples/optics_shap_on_structure_map_ex_2026-03-10_13-19-02/ex_screenshot.png?raw=true" alt="example output of OPTICS SHAP map" style="width:75%; height:75%;">

_Example output with SHAP values mapped to structure by OPTICS_

This script takes the output CSVs from the SHAP analysis script (both individual explanations AND pairwise comparison differences) and maps the importance values onto a 3D protein structure (PDB).

It modifies the B-factor column of the PDB file, allowing you to visualize "importance" as a heat map (Blue=Low, Orange=High importance). For comparison outputs, it maps the absolute difference in SHAP values to highlight the regions most responsible for the divergence in predicted λmax.

```
Required Args:
  -s, --shap_csv: Path to the SHAP analysis CSV file generated by optics_shap.py.

Optional Args:
  -p, --pdb_file: Path to PDB file(s) or ID(s). Can provide two comma-separated paths (e.g., struct1.pdb,struct2.pdb) if mapping 'both' sequences from a comparison CSV. Default: 1U19 (Bovine Rhodopsin).

  -o, --output_dir: Output directory.

  --chain: Chain ID to map to. Default: A.

  --use_query_position: Check this if your CSV uses target sequence numbering rather than reference numbering.

  --comp_target: If mapping a pairwise comparison CSV using query positions, select which sequence's numbering to map ('1', '2', or 'both'). Default: 1.

  --map_bovine_also: If using a custom PDB, this flag forces a second output mapped to Bovine Rhodopsin (1U19) for comparison.

  --top_n_labels: Number of top SHAP sites to automatically label in the generated visualization script. Default: 10.

  --software: The target software for the visualization script output ('pymol' or 'chimerax'). Default: 'chimerax'.

```

**Example Command #1: SHAP Importance Structure Mapping For Single Sequence (no comparison metrics)**

```
python optics_structure_map.py -s ./examples/optics_shap_on_structure_map_ex_2026-03-10_13-19-02/C_phantasticus_LWS1_shap_analysis.csv -p ./examples/ex_structures/C_phantasticus_LWS1_esmfold.pdb --top_n_labels 10 --software chimerax --map_bovine_also
```

**Example Command #2: SHAP Comparison Mapping on Both Sequences with Top 5 Positions Labeled**

```
python optics_structure_map.py -s ./examples/optics_shap_on_structure_map_ex_2026-03-10_13-19-02/C_phantasticus_LWS1_vs_C_phantasticus_LWS2_shap_data.csv -p ./examples/ex_structures/C_phantasticus_LWS1_esmfold.pdb,./examples/ex_structures/C_phantasticus_LWS2_esmfold.pdb --use_query_position --comp_target both --top_n_labels 10 --software chimerax
```

### Output

- SHAP Annotated Structure File (PDB): Generates a .pdb file with importance (or difference) scores in the B-factor column. If comp_target is both, generates two PDB files using the sequence numbering schemes and provided PDB templates.

- Visualization Script (.pml or .cxc): Generates a PyMOL or ChimeraX script designed to automatically color the heat-map properly and display text labels over the top SHAP sites.

### 4\. Generate Custom Structure Annotations (`optics_structure_annotations.py`)

A general-purpose tool to visualize arbitrary annotations (e.g., mutation sites, binding pockets) on a structure. It takes a simple CSV and creates a runnable visualization script.

```
Required Args:
  -a, --annotation_file: CSV/TSV file with columns: 'position' (required), 'color' (optional), 'style' (optional), 'label' (optional).

Optional Args:
  -p, --pdb: PDB ID or path. Default: 1U19.

  -o, --output_dir. Default: '.'

  --software: Target visualization software ('pymol' or 'chimerax'). Default: 'chimerax'.

  --chain: Chain identifier. Default: 'A'

```

**Example Command:**

```
python optics_structure_annotations.py -a ./examples/optics_custom_annotations_ex.csv -p 1U19 --software chimerax
```

### Output

- Custom Annotation Script (ChimeraX or PyMol): Generates a ChimeraX or PyMol specific visualization script. Typically you can just open these if your protein structure of interest is in the same folder.

### 5\. Using the OPTICS GUI

That's right! No-need for command line, OPTICS can also be used as a GUI! 
The usage is quite simple, just use the command below (with your OPTICS conda enviornment activated) and get to predicting. ;)

**To run the GUI:**

```
python run_optics_gui.py
```

<img src="https://github.com/VisualPhysiologyDB/optics/blob/main/data/logo/optics_gui_ex.png?raw=true" alt="ex optics gui" style="width:65%; height:65%;">

The GUI provides tabs/buttons for all four major pipelines:

1.  **Standard Predictions**: Run the main λmax prediction workflow.

2.  **SHAP Interpretation**: Run feature attribution analysis.

3.  **Structure Mapping**: Map SHAP values to PDB files.

4.  **Structure Annotations**: Visualize custom data on structures.

---

Understanding the λmax Prediction Models
----------------------------------------

The `--model` flag allows you to select a specific pre-trained model. Each model is named after the data-subset it was trained on. 

To keep the base installation lightweight, models are divided into **Core** and **Extra** categories.

### **Core Models (Included by Default)**

These models are included out-of-the-box when you clone this repository:

-   **```whole-dataset```**: Trained on the entire VPOD dataset. **Recommended**.

-   **```whole-dataset-mnm```**: Trained on the entire dataset including "Mine-n-Match" inferred data.

-   **```wildtype```**: Trained exclusively on wild-type sequences.

-   **```wildtype-mnm```**: Trained on wild-type sequences including "Mine-n-Match" inferred data.

-   **```type-one```**: Trained on the Type-One (Microbial) opsin dataset (previously published by [**Karyasuyama et al. 2018**](10.1038/s41598-018-33984-w))

### **Extra Models (Requires Separate Download)**

We also offer specialized taxonomic and mutational subset models. Because of file size constraints, these are hosted in a separate repository.

-   **```vertebrate```** & **```vertebrate-mnm```**

-   **```invertebrate```** & **```invertebrate-mnm```**

-   **```wildtype-vert```** & **```wildtype-vert-mnm```**

-   **```wildtype-mut```**

> **📥 How to get Extra Models:** To use any of the extra models, please visit the [**Extra OPTICS Models Repository**](https://github.com/VisualPhysiologyDB/extra_optics_models/tree/main "null"). Download the required `.pkl` files and place them in your local `models/reg_models/` and/or `models/bs_models/` directories as instructed there.

### **The `-mnm` Suffix**

Models ending in **-mnm** (e.g., `wildtype-mnm`) are trained on augmented datasets.

-   **Standard models**: Trained *exclusively* on heterologous expression data (in-vitro).

-   **-mnm models**: Trained on heterologous data *plus* data inferred via our **"Mine-n-Match"** procedure (in-vivo correlations). See [Frazer et al. 2025](https://doi.org/10.1101/2025.08.22.671864 "null") for details.
  - _**Note** - These models should be treated as secondary to heterolgous models for now. The heterolgous models are the 'gold-standard' and MNM models are the 'silver-standard'. Still useful, but not equal._       
---
## License
All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies

## Citation

- **IF citing this GitHub and its contents use the following DOI provided by Zenodo...**

      10.5281/zenodo.10667840
    
- **IF you use OPTICS in your research, please cite the following paper(s):**

  - Our more recent publication directly on the making/utility of OPTICS.

        Seth A. Frazer, Todd H. Oakley. Accessible and Robust Machine Learning Approaches to Improve the Opsin Genotype-Phenotype Map. bioRxiv, 2025.08.22.671864. https://doi.org/10.1101/2025.08.22.671864
    
  - Our original paper on the development of VPOD; the opsin genotype-phenotype database backbone for training the ML models used in OPTICS. 

        Seth A. Frazer, Mahdi Baghbanzadeh, Ali Rahnavard, Keith A. Crandall, & Todd H Oakley. Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD). GigaScience, 2024.09.01. https://doi.org/10.1093/gigascience/giae073

## Contact
Contact information for author questions or feedback.

  **Todd H. Oakley** - [ORCID ID](https://orcid.org/0000-0002-4478-915X)
    
    oakley@ucsb.edu
    
**Seth A. Frazer** - [ORCID ID](https://orcid.org/0000-0002-3800-212X)

    sethfrazer@ucsb.edu
    
---
## Additional Notes/Resources

- Want to use OPTICS without the hassle of the setup? -> [CLICK HERE](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest) to visit our Galaxy Project server and use our tool!

- *OPTICS v1.3 uses VPOD_v1.3 for training.*

- **[Here](https://tinyurl.com/u7hn9adm)** is a link to a bibliography of the publications used to build VPOD_v1.2 (VPOD_v1.3 version not yet released)
  
- If you know of publications for training opsin ML models not included in the VPOD_v1.2 database, please send them to us through **[this form](https://tinyurl.com/29afaxyr)**
  
- Check out the **[VPOD GitHub](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db)** repository to learn more about our database and ML models!
