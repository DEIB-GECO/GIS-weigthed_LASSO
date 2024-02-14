# GIS prior knowledge
Enhancing functional interpretability in gene expression data analysis by prior-knowledge incorporation

## Description 
We developed an integrative approach to feature selection that combines weighted LASSO feature selection and prior biological knowledge in a single step by means of a novel score of biological relevance that summarizes information extracted from popular biological knowledge bases.

## Application Use Cases
We compared the performance of the standard regularized LASSO model and our proposed approach on two application use cases concerning the cancer-related subtype prediction of patients based on gene expression data. The use cases concern the the classification of Breast Invasive Carcinoma patients and Colorectal Cancer patients in their corresponding cancer subtypes. 

## Implementation
To perform the weighted LASSO in Python using the scikit-learn library, we modified the corresponding functions using the development version of scikit-learn. 
The modified package is available at https://github.com/SofSof98/scikit-learn-lasso/tree/weightedlasso.
After cloning the repository, build a dedicated environment with:
'''
conda create -n sklearn-env -c conda-forge python=3.9 numpy scipy cython=0.29.33 
conda activate sklearn-env 
'''

Then build the scikit-learn package with:
'''
cd scikit-learn-lasso 
pip install -v --no-use-pep517 --no-build-isolation -e . 
'''

Lastly, install the required packages from requirements.txt


We computed the score of biological relevance using previous versions of the knowledge bases. To download and use the updated versions, run the 'get_updated_prior_knowldge.ipynb' notebook. 
