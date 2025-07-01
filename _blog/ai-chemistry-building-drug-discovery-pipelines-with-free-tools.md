---
title: 'AI + Chemistry: Building Drug Discovery Pipelines with Free Tools'
date: 2025-07-01
permalink: /blog/2025/07/ai-chemistry-building-drug-discovery-pipelines-with-free-tools
excerpt_separator: <!--more-->
toc: true
tags:
  - drug discovery
  - cheminformatics
---

The future of drug discovery is open, smart, and community-driven. Gone are the days of relying solely on expensive, proprietary platforms, today’s open-source AI tools are transforming how we explore chemical space, and they’re accessible to all.  
<!--more-->

Whether you're screening billions of compounds or predicting molecular properties, a rich ecosystem of open-source libraries and data is leveling the playing field for startups, academics, and even “indie” scientists. To illustrate this revolution, let’s walk through a drug discovery pipeline for a real target – the ABL kinase (the c-Abl tyrosine kinase, infamous as the target of the leukemia drug imatinib). We’ll see how open tools empower every step, from data to models to visualization, with code snippets showing these tools in action.

## From Data to Discovery: The Open-Source Pipeline (Case Study: ABL Kinase)

To ground things, imagine we’re hunting for new inhibitors of the ABL kinase, a critical enzyme in cancer (Bcr-Abl causes chronic myeloid leukemia when mutated). We want to:  
1. Gather known bioactivity data for ABL  
2. Featurize and analyze molecules  
3. Train an AI model to predict new inhibitors  
4. Convert and prepare compounds for simulation  
5. Visualize how they bind  

Open resources make this feasible:  

### Open Data (ChEMBL)  
We can retrieve ABL bioactivity data from ChEMBL, a large open database of drug-like molecules and their biological activities. (ChEMBL contains millions of measured compound-target activities – over 5.4 million bioactivity data points for 1M+ compounds as of 2012, and it’s grown even larger since!). This gives us a training dataset of known ABL inhibitors and non-inhibitors without any paywalls. We can use ChEMBL’s web services or downloads to get, say, all compounds tested against ABL (IC50, Ki values, etc.), then use Pandas to filter and tabulate the data.  

### Data Handling (Pandas & NumPy)  
With our ABL dataset in hand (e.g. as a CSV of SMILES and activity labels), we use Pandas to clean and manipulate it and NumPy for any numerical computing. These “classics” form the backbone of any custom pipeline – e.g., grouping data, normalizing values, splitting into train/test sets. They might not be drug discovery-specific, but their flexibility is indispensable. We might do:

```python
import pandas as pd
df = pd.read_csv("ABL_bioactivity.csv")
```

---

## RDKit – The Unsung Hero of Cheminformatics

If you're doing anything with molecules in Python, chances are RDKit is working behind the scenes. RDKit is an open-source cheminformatics toolkit widely used for tasks like generating molecular fingerprints, performing substructure searches, computing descriptors, and manipulating chemical structures.  

In our ABL example, we use RDKit to:  
- Generate Morgan fingerprints  
- Perform substructure searches  
- Compute molecular descriptors  

### Code Example:

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

imatinib_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
nilotinib_smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)C(=O)Nc4cc(cc(c4)n5cc(nc5)C)C(F)(F)F"

mol1 = Chem.MolFromSmiles(imatinib_smiles)
mol2 = Chem.MolFromSmiles(nilotinib_smiles)

fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)

sim = DataStructs.TanimotoSimilarity(fp1, fp2)
print(f"Tanimoto similarity: {sim:.2f}")

query = Chem.MolFromSmarts("c1ccc(cc1)Nc2nccc(n2)")
match = mol1.HasSubstructMatch(query)
print("Substructure match:", match)
```

---

## DeepChem – AI Made Beautifully Simple

DeepChem is an open-source library that brings advanced models like GCNs and multitask networks to your fingertips. It’s built on TensorFlow/PyTorch but hides the complexity.

### Code Example:

```python
import deepchem as dc
import numpy as np

smiles_list = [
    "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C",  # active
    "CCOC(=O)c1ccc(cc1)N",  # inactive
]
labels = np.array([1, 0])

featurizer = dc.feat.MolGraphConvFeaturizer()
X = featurizer.featurize(smiles_list)
y = labels

dataset = dc.data.NumpyDataset(X, y)

model = dc.models.GraphConvModel(n_tasks=1, mode='classification', metrics=[dc.metrics.Metric(dc.metrics.roc_auc_score)])
model.fit(dataset, nb_epoch=20)

pred_probs = model.predict(dataset)
print(pred_probs)
```

---

## Open Babel – Convert Like a Pro

Open Babel helps switch between formats like SMILES, SDF, PDB, etc.

### Code Example:

```python
from openbabel import pybel

smiles = "Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C"
mol = pybel.readstring("smi", smiles)
mol.addh()
mol.make3D()
mol.write("sdf", "imatinib_3D.sdf", overwrite=True)
```

---

## PyMOL – Visualize What AI Discovers

PyMOL is great for inspecting protein–ligand complexes and generating publication-quality figures.

### Code Example:

```python
import pymol2

with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd
    cmd.load("ABL_kinase.pdb", "protein")
    cmd.load("imatinib_3D.sdf", "ligand")
    cmd.hide("everything")
    cmd.show("cartoon", "protein")
    cmd.show("sticks", "ligand")
    cmd.zoom("ligand", 5)
    cmd.png("abl_imatinib.png", width=800, height=600, ray=1)
```

---

## Chemprop – Graph Neural Networks Made Easy

Chemprop offers fast training of MPNNs for tasks like QSAR and virtual screening.

### CLI Example:

```bash
chemprop_train --data_path abl_activity.csv --smiles_column smiles --target_columns active \
               --dataset_type classification --save_dir abl_model
```

### Python Example:

```python
from chemprop.train import run_training

params = {
    "data_path": "abl_activity.csv",
    "smiles_column": "smiles",
    "target_columns": ["active"],
    "dataset_type": "classification",
    "save_dir": "abl_model",
    "epochs": 30
}
run_training(params)
```

---

## The Classics: Pandas, NumPy, Scikit-Learn – Data Science Backbone

These libraries handle everything from preprocessing to baseline models.

### Example Random Forest:

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

X = [list(map(int, fp.ToBitString())) for fp in [fp1, fp2]]
y = labels

X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=0)
model = RandomForestClassifier(n_estimators=100).fit(X_train, y_train)
print("Validation accuracy:", model.score(X_val, y_val))
```

---

## The Open-Source Revolution in Action

As we’ve seen, an end-to-end drug discovery pipeline can now be constructed with open-source tools at every step, with each tool excelling in its domain:  

- Data acquisition from open databases (ChEMBL, PubChem, etc.) provides the fuel.  
- RDKit ensures we can manipulate and understand chemical structures easily.  
- DeepChem and Chemprop bring powerful AI models to make predictions or generate new hypotheses.  
- Open Babel makes sure our data can move anywhere it needs to (no format silos).  
- PyMOL lets us visualize and validate the AI’s suggestions in the context of 3D biology.  
- Pandas/NumPy/Sci-kit tie everything together with data handling and auxiliary analyses.  

Importantly, all of these are free and open. There are no license fees or onerous contracts – you can install them on your laptop right now and get to work. This open-source ecosystem is accelerating innovation by putting incredible power in the hands of anyone with an idea and a bit of coding knowledge. It lowers the barrier to entry for drug discovery projects: a small startup or an academic lab can now deploy workflows that rival those in big pharma, without spending a fortune on software.  

The open-source revolution is not just about cost savings; it’s about community and collaboration. These tools are constantly improving through contributions from users worldwide. For example, new algorithms and best practices in chemoinformatics are often rapidly integrated into RDKit or DeepChem by community members. If a feature is missing, you can add it or request it. This fosters an environment where researchers share not just data, but also the methods to analyze that data – leading to more reproducible and transparent science.  

What’s in your AI drug discovery toolkit? Chances are, if you start exploring, you’ll end up with many of the open-source tools above in your repertoire. And you’ll be joining a movement that is driving the future of pharmaceutical innovation. No paywalls, no red tape – just raw potential and a community eager to push the boundaries of what AI can do for medicine. The open-source revolution in drug discovery is here, and it’s incredibly exciting. Get involved, experiment with these tools, and who knows – you might discover the next Halicin or imatinib, and you’ll have the open-source community cheering you on every step of the way.  
