# GIS-weighted LASSO
Enhancing functional interpretability in gene expression data analysis by prior-knowledge incorporation


## Description 
We developed an integrative approach to feature selection that combines weighted LASSO feature selection and prior biological knowledge in a single step by means of a novel score of biological relevance that summarizes information extracted from popular biological knowledge bases.

## Application Use Cases
We compared the performance of the standard regularized LASSO model and our proposed approach on two application use cases concerning the cancer-related subtype prediction of patients based on gene expression data. The use cases concern the classification of Breast Invasive Carcinoma (BRCA) patients and Colorectal Cancer (CRC) patients in their corresponding cancer subtypes. We also performed two distinct sensitivity analyses to evaluate the impact of incorporating our proposed score of biological relevance into LASSO regularization. We used a controlled dataset with limited correlation among the features for these analyses, considering publicly available RNA-seq profiles of Kidney Renal Clear Cell Carcinoma patients from The Cancer Genome Atlas (TCGA) project. The preprocessed dataset is available [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/data/data_kidney), along with the list of features considered in the controlled dataset. 
For all datasets analysed, the data are preprocessed as described in the notebooks found [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/notebooks).



## Implementation
To perform the GIS-weighted LASSO in Python using the scikit-learn library, we modified the corresponding functions using the development version of scikit-learn. 
The modified package is available [here](https://github.com/SofSof98/scikit-learn-lasso/tree/weightedlasso).
After cloning the repository, build a dedicated environment with:
```bash
conda create -n sklearn-env -c conda-forge python=3.9 numpy scipy cython=0.29.33 
conda activate sklearn-env
```

Then, build the scikit-learn package with:
```bash
cd scikit-learn-lasso 
pip install -v --no-use-pep517 --no-build-isolation -e . 
```
Lastly, install the required packages from requirements.txt


We computed the score of biological relevance using the specific versions of the knowledge bases, which can found [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/data/prior_knowledge). These versions are the following:
- GO (format-version 1.2, release date: 2023-03-06) 
- Reactome (version V85)
- HPO (format-version 1.2, release date: 2023-09-01)

To download and use the updated versions, run [this notebook](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/blob/main/notebooks/get_updated_prior_knowldge.ipynb). We used the [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) R library (version 2.56.1) to extract all available GO annotation terms. To obtained updated GO annotation terms for each gene run the R script [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/src/R).

## Additional Information
All the code is available [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/src). To replicate the experiments run the following scripts:
- `experiments_BRCA.py` for experiments on the BRCA dataset
- `experiments_CRC.py` for experiments on the CRC dataset

The code to replicate the two sensitivity analyses is available [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/src/gis_sensitivity_analysis).
All the results from the experiments we performed can be found [here](https://github.com/DEIB-GECO/GIS-weigthed_LASSO/tree/main/results).
