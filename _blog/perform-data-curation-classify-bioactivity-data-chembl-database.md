---
title: 'How to perform data curation and classify bioactivity data on ChEMBL database'
date: 2022-09-07
permalink: /blog/2022/09/perform-data-curation-classify-bioactivity-data-chembl-database
excerpt_separator: <!--more-->
toc: true
tags:
  - machine learning
  - cheminformatics
  - career
---

The ChEMBL database is a manually curated resource of bioactive molecules with drug-like properties, integrating chemical, bioactivity, and genomic data to facilitate the translation of genomic information into new drugs.
<!--more-->
As of February 2022, ChEMBL version 30 contains over 2.2 million compound records, 1.5 million assays, and spans 14,000 targets, 2,000 cells, and 43,000 indications. The database has continued to grow significantly in both scope and scale, with the most recent release, ChEMBL 33, in May 2023, containing 2,786,911 compound records. ChEMBL is an essential resource for the scientific community, enabling the investigation of various health-related and scientific questions [1-3].

## 1. Data search
The first step to use this database is to search the ChEMBL database using keywords of a target protein of interest, it is possible to run a search using other keywords related to diseases, compounds or assays. In this tutorial, we are going to search for Acetylcholinesterase as illustrated in Figure 1.

| ![Figure1](https://user-images.githubusercontent.com/7014404/225256089-5c72e4cc-f77c-4ee8-90ea-f548343b5ac4.png) |
|:--:|
| <b>Figure 1: ChEMBL search result example</b> |


Notice that our search resulted on 24 targets, it is important to choose the right protein for the right organism of the study of interest. For this example, we are interested in human Acetylcholinesterase corresponding to the ID: CHEMBL220. 
After clicking on the target, we will be sent to another page containing all the data concerning the selected target such as: name and classification, drugs and clinical candidates, activity charts, etc. 
Scroll down to activity charts and notice the pie chart on the left concerning all the associated bioactivity data compiled from the literature and their distribution according to the activity type.  
 
 
| ![Figure2](https://user-images.githubusercontent.com/7014404/225256100-351cffcb-dc10-4d4a-949c-07397f6b4bb6.png) |
|:--:|
| <b>Figure 2: Activity charts and distribution of activity types of the selected target, CHEMBL220</b> |

Upon observation of the activity chart, we can quickly determine which activity type is the most reported in the literature, in this case it refers to half-maximal inhibitory concentrations (IC50) which have been reported 8205 times. 
Once we click on the desired activity type, we can download the entire dataset in CSV or TSV Format containing various informations such as ChEMBL ID for each compound, SMILES, Standard Type and Standard Value referring to the activity type and value respectively.## 

## 2. Data curation
Note that it is necessary to remove any unwanted data before proceeding with data curation. In this case, we are only interested in the compound’s IDs, Smiles, Standard Type and Standard Value. It is possible to perform this task with any CSV reader such as Google Sheets or Microsoft Excel. 
Once we have performed the primary cleaning on our data, we can import it on Google Colab or Jupyter Notebook using the code below:


**Import necessary libraries**
```
import pandas as pd
```
**Read the dataset**
```
x=pd.read_csv('ache.csv')
```
**Display the dataset**
```
x
```

**Output**

| ![Figure3](https://user-images.githubusercontent.com/7014404/225256101-c4f9bfc6-652b-46da-af63-dbcab7c91255.png) |
|:--:|
| <b>Figure 3: AChE curated dataset output.</b> |
 


### 2.1. Remove duplicate compounds

When dealing with a large dataset of compounds, it is very likely to find a great deal of duplicates with different reported activities due to different conditions of each laboratory. However, it is possible to deal with this issue by averaging all reported activity by calculating their mean values using the code below:

```
x["mean_value"]=x[['Molecule ChEMBL ID', 'Smiles','Standard Type','Standard Value']].groupby(['Molecule ChEMBL ID'])['Standard Value'].transform("mean")
```
The next step is to merge all the duplicate compounds into one, for this reason we can use the code below to remove all duplicates while keeping only the first one.
```
x=x.drop_duplicates("Molecule ChEMBL ID", keep="first")
```

**Output:**

| ![Figure4](https://user-images.githubusercontent.com/7014404/225256109-4ce6979f-0e42-469a-8191-459dc6530c3b.png) |
|:--:|
| <b>Figure 4: Pandas output of AChE dataset after removing duplicate compounds.</b> |
 
 
It is possible to find some compounds on the dataset with no available activity, notice on the sheet above that some activities are marked with “NaN” which stands for “Not a Number” in computer science, it is therefore necessary to remove them before proceeding. We can simply run the code below:

```
x=x.dropna()
x
```

**Output:**

| ![Figure5](https://user-images.githubusercontent.com/7014404/225256112-91eb87c4-cbd7-424c-80fe-bb45ed014f31.png) |
|:--:|
| <b>Figure 5: Final curated dataset.</b> |

## 3. Data classification

Once we have curated our data, now it is possible to classify compounds in order to apply it for machine learning classification models. For this reason, we need to define an activity cutoff to define our active and inactive compounds. In the case of enzyme inhibition, the literature indicates that most potent enzyme inhibitors have activities in the nanomolar range. For this reason, we can proceed by setting a threshold of 1000 nM corresponding to 1 μM or lower for defining our active compounds. 

Run the code below to create a variable with all active compounds:

```
active=x.loc[x['mean_value']<=1000]
```

We do the same for inactive compounds by setting a cutoff of 10 000 nM (10 μM) or higher using the code below:

```
inactive=x.loc[x['mean_value']>10000]
```

## 4. Data labelling

Now that we have defined our active and inactive compounds, it is necessary to label the data in order to combine the entire dataset. We will simply refer to active compounds as “1” and inactive compounds as “0”.
Run the code below:

```
active["Class"]=1
inactive["Class"]=0
```


Now we can proceed to combining the entire dataset.

Run the code below:

```
combined=pd.concat([active,inactive],axis=0)
combined
```


**Output:**
 
| ![Figure6](https://user-images.githubusercontent.com/7014404/225256117-a1a5f883-f733-4679-9ac5-668d4a82e708.png) |
|:--:|
| <b>Figure 6: Curated dataset with label column indicating whether the compound is active or inactive.</b> |

Finally, we can save our dataset for further use.
Run the code below:

```
combined.to_csv("ache_labelled.csv", index=None)
```


## 5. Bottom line

This article’s aim was to demonstrate an alternative way to retrieve bioactivity data from ChEMBL without using code. Furthermore, data curation and data classification was covered in detail as it is a necessary step and can highly impact the performance of machine learning models. If you found this article useful, follow the blog for more tutorials in the future.

## 6. References

[1] https://chembl.blogspot.com/2022/03/chembl-30-released.html

[2] https://www.ebi.ac.uk/chembl/id_lookup/CHEMBL14116/

[3] https://academic.oup.com/nar/article/52/D1/D1180/7337608
