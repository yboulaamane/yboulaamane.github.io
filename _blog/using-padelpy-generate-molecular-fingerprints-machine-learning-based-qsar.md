---
title: 'Using PaDELPy to Generate Molecular Fingerprints for Machine Learning-Based QSAR'
date: 2023-03-15
permalink: /blog/2023/03/using-padelpy-generate-molecular-fingerprints-machine-learning-based-qsar
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - drug discovery
  - qsar
---

PaDELPy is a Python library that integrates the PaDEL-Descriptor molecular descriptor calculation software, allowing efficient generation of molecular fingerprints for machine learning-based quantitative structure-activity relationship (QSAR) models in drug discovery.

Machine learning models, created by training algorithms to recognize data patterns, can be either supervised or unsupervised, applied in classification, regression, and more. Here, we’ll explore using PaDELPy to generate fingerprints that are crucial in predictive QSAR modeling, specifically targeting molecular activity prediction within the HCV Drug dataset.
<!--more-->

![Molecular Fingerprint](https://user-images.githubusercontent.com/7014404/225259643-df0568cd-1cfe-4395-aa7e-980902108f25.png)

## What Is PaDELPy?
PaDELPy is a Python wrapper for the Java-based PaDEL-Descriptor software, streamlining molecular descriptor calculations. With PaDELPy, you can easily compute molecular fingerprints, avoiding the complexity of Java setup and reducing the time required for implementation.

## Getting Started with the Code
In this tutorial, we will create a machine learning model using Random Forest to predict molecular activity within the HCV Drug dataset. The dataset is available [here](https://github.com/chaninlab/hcvpred/blob/master/HCV_NS5B_Curated.csv). 

To install the PaDELPy library, use the following:
```
# Installing the library
!pip install padelpy
```
Download and configure XML data files required by PaDELPy:
```
# Downloading XML data files
!wget https://github.com/dataprofessor/padel/raw/main/fingerprints_xml.zip
!unzip fingerprints_xml.zip
```
```
# Listing and sorting downloaded files
import glob
xml_files = glob.glob("*.xml")
xml_files.sort()
xml_files
```
Output:
```
['AtomPairs2DFingerprintCount.xml', 'AtomPairs2DFingerprinter.xml', 'EStateFingerprinter.xml', ...]
```
Create a list of available fingerprint types:
```
# Creating a list of present files
FP_list = ['AtomPairs2DCount', 'AtomPairs2D', 'EState', 'CDKextended', 'CDK', 'CDKgraphonly', 
           'KlekotaRothCount', 'KlekotaRoth', 'MACCS', 'PubChem', 'SubstructureCount', 'Substructure']
```
Create a data dictionary of file names for easy reference:
```
# Creating Data Dictionary
fp = dict(zip(FP_list, xml_files))
fp
```
After setting up, load the dataset:
```
# Loading the dataset
import pandas as pd
df = pd.read_csv('https://raw.githubusercontent.com/dataprofessor/data/master/HCV_NS5B_Curated.csv')
df.head()
```
Prepare the data by concatenating necessary columns:
```
# Concatenating necessary columns
df2 = pd.concat([df['CANONICAL_SMILES'], df['CMPD_CHEMBLID']], axis=1)
df2.to_csv('molecule.smi', sep='\t', index=False, header=False)
```
Select a fingerprint type and calculate descriptors:
```
# Setting the fingerprint module
from padelpy import padeldescriptor
fingerprint = 'Substructure'
fingerprint_output_file = ''.join([fingerprint, '.csv'])  # Substructure.csv
fingerprint_descriptortypes = fp[fingerprint]

padeldescriptor(mol_dir='molecule.smi', 
                d_file=fingerprint_output_file,
                descriptortypes=fingerprint_descriptortypes,
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True)
```
Display calculated fingerprints:
```
descriptors = pd.read_csv(fingerprint_output_file)
descriptors.head()
```
Creating a Random Forest Model
Using the processed data, create a Random Forest model for classification.
```
X = descriptors.drop('Name', axis=1)
y = df['Activity']  # Target variable

# Removing low variance features
from sklearn.feature_selection import VarianceThreshold

def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]

X = remove_low_variance(X, threshold=0.1)
```
Split data into training and testing sets:
```
# Splitting into Train and Test sets
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
```
Train the Random Forest model:
```
from sklearn.ensemble import RandomForestClassifier
model = RandomForestClassifier(n_estimators=500, random_state=42)
model.fit(X_train, y_train)
```
Evaluate model performance using Matthews Correlation Coefficient (MCC):
```
# Making predictions
y_train_pred = model.predict(X_train)
y_test_pred = model.predict(X_test)
```
```
# Calculating MCC for train and test
from sklearn.metrics import matthews_corrcoef
mcc_train = matthews_corrcoef(y_train, y_train_pred)
mcc_test = matthews_corrcoef(y_test, y_test_pred)
```
Perform 5-fold cross-validation:
```
import pandas as pd
performance_metrics = pd.DataFrame({
    'Model': ['Random Forest'],
    'MCC_Train': [mcc_train],
    'MCC_CV': [mcc_cv],
    'MCC_Test': [mcc_test]
})
performance_metrics
```
Consolidate performance metrics into a single DataFrame:
```
import pandas as pd
performance_metrics = pd.DataFrame({
    'Model': ['Random Forest'],
    'MCC_Train': [mcc_train],
    'MCC_CV': [mcc_cv],
    'MCC_Test': [mcc_test]
})
performance_metrics
```
## Conclusion
In this tutorial, we explored using PaDELPy to calculate molecular fingerprints, then developed a Random Forest model to predict molecular drug activity. The high Matthews Correlation Coefficient values suggest that this model is effective on the current dataset, though other algorithms could also be evaluated for further optimization.




