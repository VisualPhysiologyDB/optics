**Code**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) **Data**: [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  **VPOD_1.2 DOI**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12213246.svg)](https://doi.org/10.5281/zenodo.12213246)


# Opsin Phenotype Tool for Inference of Color Sensitivity (OPTICS) 

![](https://github.com/VisualPhysiologyDB/optics/blob/main/examples/optics_on_msp_test_of_optics_2024-09-19_19-48-49/msp_bs_viz.svg?raw=true)

  _Example Box Plot Output for Bootstrap Predictions of Opsin λmax by OPTICS_

---
## Description

- **OPTICS** is an open-source tool that predicts the Opsin Phenotype (λmax) from unaligned opsin amino-acid sequences. 
- **OPTICS** leverages machine learning models trained on the Visual Physiology Opsin Database (VPOD).
- **OPTICS** is also avaliable as an online tool [**here**](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest), hosted on our [**Galaxy Project**]([https://usegalaxy.org/](http://galaxy-dev.cnsi.ucsb.edu:8080/)) server.

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
   B. Use the 'requirements.txt' file to download base package dependencies for OPTICS
   ```bash
   pip install -r requirements.txt
   ```
   - THEN
     
   ```bash
   conda activate optics_env
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
   
4. **Usage**

    - **MAKE SURE YOU HAVE ALL DEPENDENCIES DOWNLOADED ARE IN THE FOLDER DIRECTORY FOR OPTICS BEFORE RUNNING ANY SCRIPTS!**

    - Parameters
      ```
        -in - FASTA file containing unaligned opsin sequences.

        -rd - Name for job; used to create output folder. (deafult = optics_on_unamed_{date_and_time_label})

        -out - Name for output file. (deafult = 'optics_predictions.txt')

        -m - Select model to use for prediction. Options are 'whole-dataset', 'vertebrate', 'invertebrate', 'wildtype', or 'wildtype-vert'

        -e - Select preferred encoding method used to train model and make predictions. Options are 'one-hot' or 'aa_prop'

        -b - Option to enable/disable Blastp analysis on query sequences. [True/False]

        -ir - Name for the blastp report output file. (deafult = 'blastp_report.txt')

        -r - Select reference sequence used for position numbering for blastp analysis. Options are 'bovine', 'squid', or 'custom'

        -f - Custom reference sequence file used for blastp analysis - **ONLY NEED TO PROVIDE IF SELECTED 'CUSTOM' FOR REFERENCE SEQUENCE**

        -s - Option to enable/disable bootstrap predictions on query sequences. **NOTE: VISULIZATION ONLY PRODUCED FOR 10 SEQUENCES OR LESS**

        -bsv - Name for the pdf file output file for visualizing bootstrap predictions. (default = 'bootstrap_viz.pdf')
       ```     
     - Example Command Line Usage vvv
       ```bash
        python optics_predictions.py -in ./examples/msp_erg_raw.txt -rd msp_test_of_optics -out msp_predictions.tsv -m whole-dataset -e aa_prop -b True -ir msp_blastp_report.tsv -r squid -s True -bsv msp_bs_viz.pdf 

## Input

- **Unaligned** FASTA file containing opsin amino-acid sequences.
- Example FASTA Entry:
  ```
    >NP_001014890.1_rhodopsin_Bos taurus
    MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRT 
    PLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVC 
    KPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVV 
    HFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQG 
    SDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA   
  ``` 
## Output

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

      Seth A. Frazer, Mahdi Baghbanzadeh, Ali Rahnavard, Keith A. Crandall, & Todd H Oakley. (2024). Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD). GigaScience, 2024.09.01 https://doi.org/10.1093/gigascience/giae073

## Contact
Contact information for author questions or feedback.

  **Todd H. Oakley** - [ORCID ID](https://orcid.org/0000-0002-4478-915X)
    
    oakley@ucsb.edu
    
**Seth A. Frazer** - [ORCID ID](https://orcid.org/0000-0002-3800-212X)

    sethfrazer@ucsb.edu
    
---
## Additional Notes/Resources

- Want to use OPTICS without the hassle of the setup? -> [CLICK HERE](http://galaxy-dev.cnsi.ucsb.edu:8080/?tool_id=optics_1&version=latest) to visit our Galaxy Project server and use our tool!

- *OPTICS v1.0 uses VPOD v1.2 for training.*

- **[Here](https://tinyurl.com/u7hn9adm)** is a link to a bibliography of the publications used to build VPOD_1.2 (Full version not yet released)
  
- If you know of publications for training opsin ML models not included in the VPOD_1.2 database, please send them to us through **[this form](https://tinyurl.com/29afaxyr)**
  
- Check out the **[VPOD GitHub](https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db)** repository to learn more about our database and ML models!
