**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  **VPOD_v1.2 DOI**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12213246.svg)](https://doi.org/10.5281/zenodo.12213246)


# Opsin Phenotype Tool for Inference of Color Sensitivity (OPTICS) [v1.1]

![](https://github.com/VisualPhysiologyDB/optics/blob/main/examples/optics_on_full_ex_test_of_optics_2024-10-15_15-18-01/ex_bs_viz_part4.svg?raw=true)

  _Example Box Plot Output for Bootstrap Predictions of Opsin λmax by OPTICS_

---
## Description

- **OPTICS** is an open-source tool that predicts the Opsin Phenotype (λmax) from unaligned opsin amino-acid sequences. 
- **OPTICS** leverages machine learning models trained on the Visual Physiology Opsin Database (VPOD).
- **OPTICS** is also avaliable as an online tool [**here**](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest), hosted on our [**Galaxy Project**](https://usegalaxy.org/) server.

## Key Features

- **λmax Prediction**: Predicts the peak light absorption wavelength (λmax) for opsin proteins.
- **Model Selection**: Choose from different pre-trained models for prediction.
- **Encoding Methods**: Select between one-hot encoding or amino-acid property encoding for model training and prediction.
- **BLAST Analysis**: Optionally perform BLASTp analysis to compare query sequences against reference datasets.
- **Bootstrap Predictions**: Optionally enable bootstrap predictions for enhanced accuracy assessment (suggested limit to 10 sequences for bootstrap visulzations).

## Installation

1. **Clone the repository:**
   ```bash
    git clone https://github.com/VisualPhysiologyDB/optics.git

2. **Install dependencies:** [Make sure you are working in the repository directory from here-after]

   A. Create a Conda environment for OPTICS (make sure you have [Conda](https://www.anaconda.com/) installed)
   ```bash
   conda create --name optics_env python=3.11 
   ```
   ### THEN
   ```bash
   conda activate optics_env
   ```
   B. Use the 'requirements.txt' file to download base package dependencies for OPTICS
   ```bash
   pip install -r requirements.txt
   ```
     

   C. **Download MAFFT and BLAST**
   
     IF working on MAC or LINUX device:

     - Install _BLAST_ and _MAFFT_ directly from the _bioconda_ channel
       ```bash
       conda install bioconda::blast bioconda::mafft
       ```
     
     IF working on WINDOWS device:
      - Manaully install the Windows compatable [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) executable on your system PATH; [the download list is here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
        - We suggest downloading '[ncbi-blast-2.16.0+-win64.exe](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-win64.exe)' 
      - You DO NOT need to download MAFFT, OPTICS should be able to run MAFFT from the files we provide when downloading this GitHub.

  
## Usage

  **MAKE SURE YOU HAVE ALL DEPENDENCIES DOWNLOADED ARE IN THE FOLDER DIRECTORY FOR OPTICS BEFORE RUNNING ANY SCRIPTS!**
  
  ### **Parameters**
     
  ```
Required Args:

-i, --input: Either a single sequence or a path to a FASTA file.

General Optional Args:

-o, --output_dir: Desired directory to save output folder/files (optional). Default: './prediction_outputs'

-p, --prediction_prefix: Base filename for prediction outputs (optional). Default: 'unnamed'

-m, --model: Prediction model to use (optional). Options: whole-dataset, wildtype, vertebrate, invertebrate, wildtype-vert, type-one, whole-dataset-mnm, wildtype-mnm, vertebrate-mnm, invertebrate-mnm, wildtype-vert-mnm. **Default: whole-dataset** 

-e, --encoding: Encoding method to use (optional). Options: one_hot, aa_prop. Default: aa_prop

BLASTp Analysis Args (optional):

--blastp: Enable BLASTp analysis.

--blastp_report: Filename for BLASTp report. Default: blastp_report.txt

--refseq: Reference sequence used for blastp analysis. Options: bovine, squid, microbe, custom. Default: bovine

--custom_ref_file: Path to a custom reference sequence file for BLASTp.  Required if --refseq custom is selected.

Bootstrap Analysis Args (optional):

--bootstrap: Enable bootstrap predictions.

--visualize_bootstrap: Enable visualization of bootstrap predictions.

--bootstrap_viz_file: Filename prefix for bootstrap visualization (PDF and SVG). Default: bootstrap_viz

  ```     
  ### Example Command Line Usage vvv
  
  ```bash
  python optics_predictions.py -i ./examples/optics_ex_short.txt -o ex_test_of_optics -p ex_predictions -m wildtype -e aa_prop --blastp -blastp_report blastp_report.txt --refseq squid --bootstrap --visualize_bootstrap --bootstrap_viz_file bootstrap_viz
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

- Predictions (TSV): λmax values, model used, and encoding method.
- BLAST Results (TXT, optional): Comparison of query sequences to reference datasets.
- Bootstrap Graphs (PDF, optional): Visualization of bootstrap prediction results.
- Job Log (TXT): Log file containing input command to OPTICS, including encoding method and model used.

  **Note** - All outputs are written into sub-folders within the 'prediction_outputs' folder, and are marked by time and date.


---
## License
All data and code is covered under a GNU General Public License (GPL)(Version 3), in accordance with Open Source Initiative (OSI)-policies

## Citation

- **IF citing this GitHub and its contents use the following DOI provided by Zenodo...**

      10.5281/zenodo.10667840
    
- **IF you use OPTICS in your research, please cite the following paper:**
  **NOTE** - We are currently working on a manuscript specific to OPTICS - so this citation will change in the near future. 

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

- *OPTICS v1.1 uses VPOD_v1.3 for training.*

- **[Here](https://tinyurl.com/u7hn9adm)** is a link to a bibliography of the publications used to build VPOD_v1.2 (Full version not yet released)
  
- If you know of publications for training opsin ML models not included in the VPOD_v1.2 database, please send them to us through **[this form](https://tinyurl.com/29afaxyr)**
  
- Check out the **[VPOD GitHub](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db)** repository to learn more about our database and ML models!
