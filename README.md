# OPTICS - Opsin Phenotype Tool for Inference of Color Sensitivity 

## Description

OPTICS is a Galaxy tool that predicts the Opsin Phenotype (λmax) from unaligned opsin amino-acid sequences. It leverages machine learning models trained on the Visual Physiology Opsin Database (VPOD).

## Key Features

- **λmax Prediction**: Predicts the peak light absorption wavelength (λmax) for opsin proteins.
- **Model Selection**: Choose from different pre-trained models for prediction.
- **Encoding Methods**: Select between one-hot encoding or amino-acid property encoding for model training and prediction.
- **BLAST Analysis**: Optionally perform BLASTp analysis to compare query sequences against reference datasets.
- **Bootstrap Predictions**: Optionally enable bootstrap predictions for enhanced accuracy assessment (limited to 10 sequences).

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://your_github_repository_url.git

2. **Install dependencies:**
   ```bash
    pip install -r requirements.txt 

3. **Usage**
    *Upload your FASTA file containing unaligned opsin sequences.
    *Select the desired model and encoding method.
    *Optionally enable BLAST analysis and/or bootstrap predictions.
    *Run the tool and retrieve the output files.
    *Command Line Example

    ```bash
    python '$__tool_directory__/prediction_functions_galaxy.py' '$input' '$output' '$blast_report' '$boot_strap_file' -m '$model' -b '$blast_checkbox' -r '$ref_sequence' -f '$ref_seq_file' -s '$boot_strap' -e '$encode_method' 

## Input

    *Unaligned FASTA file containing opsin amino-acid sequences.
    *Example FASTA Entry:

    >NP_001014890.1_rhodopsin_Bos taurus 
    MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRT 
    PLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVC 
    KPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVV 
    HFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQG 
    SDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA   
 
## Output

    *Predictions (TSV): λmax values, model used, and encoding method.
    *BLAST Results (TXT, optional): Comparison of query sequences to reference datasets.
    *Bootstrap Graphs (PDF, optional): Visualization of bootstrap prediction results.
    
## Citation

    **If you use OPTICS in your research, please cite the following paper:**

    @article{Frazer2024.02.12.579993,
    author = {Frazer, Seth A. and Baghbanzadeh, Mahdi and Rahnavard, Ali and Crandall, Keith A. and Oakley, Todd H.},
    title = {Discovering genotype-phenotype relationships with machine learning and the Visual Physiology Opsin Database (VPOD)},
    elocation-id = {2024.02.12.579993},
    year = {2024},
    doi = {10.1101/2024.02.12.579993},
    publisher = {Cold Spring Harbor Laboratory},
    URL = {[https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579993](https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579993)},
    eprint = {[https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579993.full.pdf](https://www.biorxiv.org/content/early/2024/02/14/2024.02.12.579993.full.pdf)},
    journal = {bioRxiv}
    }

## Additional Notes

    *OPTICS v1.0 uses VPOD v1.2 for training.
    *For more information about VPOD, visit:
    *VPOD Bibliography: https://tinyurl.com/u7hn9adm
    *Contribute to VPOD: https://tinyurl.com/29afaxyr
    *VPOD GitHub: https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db
    *License -
