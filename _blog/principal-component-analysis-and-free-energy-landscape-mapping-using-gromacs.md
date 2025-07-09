---
title: 'Principal Component Analysis and Free Energy Landscape Mapping Using GROMACS'
date: 2025-07-09
permalink: /blog/2024/01/principal-component-analysis-and-free-energy-landscape-mapping-using-gromacs
excerpt_separator: <!--more-->
toc: true
tags:
  - computational chemistry
  - drug discovery
---

This guide outlines how to perform Principal Component Analysis (PCA) and compute Free Energy Landscapes (FEL) from molecular dynamics (MD) simulations using GROMACS. These analyses are useful to capture dominant motions and identify energetically favorable states in biomolecular systems.
<!--more-->
---

## Part 1: Principal Component Analysis (PCA)

### Step 1: Covariance Matrix and Eigenvector Calculation

Run the following command to compute the covariance matrix and extract eigenvectors:

```bash
gmx covar -s md_0_100.tpr -f md_0_100.xtc -o eigenvalues.xvg -v eigenvectors.trr -xpma covar.xpm
```

> Note: If the command fails using `md_0_100.gro`, use `md_0_100.tpr`.

Select **protein** and **ligand** when prompted.

---

### Step 2: Determine the Number of Principal Components

Analyze `eigenvalues.xvg` to compute:

- **Total Variance**  
  `Total Variance = Σ λᵢ`

- **Explained Variance (%)**  
  `Explained Variance (%) = (λᵢ / Σ λᵢ) × 100`

- **Cumulative Variance**  
  `Cumulative Variance = Σ (λᵢ / Σ λᵢ) × 100`

Select the minimum number of components that cumulatively explain more than **50%** of the total variance.

---

### Step 3: Plot Principal Component Projections

For projecting motion along PC1 to PC5, run:

```bash
gmx anaeig -f md_0_100.xtc -s md_0_100.tpr -v eigenvectors.trr -first 1 -last 5 \
-proj pc15_lovastatina.xvg -2d project2d_s100_lovastatina.xvg -tu ns
```

- `-proj`: 1D projection along each selected eigenvector  
- `-2d`: Combined 2D projection using the first 5 components

---

## Part 2: Free Energy Landscape (FEL) Calculation

### Step 1: Generate Principal Component Projections

```bash
gmx anaeig -f md_0_100.xtc -s md_0_100.tpr -v eigenvectors.trr -last 1 -proj pc1.xvg
gmx anaeig -f md_0_100.xtc -s md_0_100.tpr -v eigenvectors.trr -first 2 -last 2 -proj pc2.xvg
```

Merge PC1 and PC2 into a single file:

```bash
paste pc1.xvg pc2.xvg | awk '{print $1, $2, $4}' > PC1PC2.xvg
```

---

### Step 2: Compute Free Energy Surface

Use GROMACS SHAM module:

```bash
gmx sham -f PC1PC2.xvg -ls FES.xpm
```

Convert `.xpm` to `.dat` using a Python script:

```bash
python2.7 xpm2txt.py -f FES.xpm -o fel.dat
```

*Python script `xpm2txt.py` converts the XPM matrix to a 3-column text file (X, Y, Energy).*

---

### Step 3: Extract Minimum Energy Conformation

After plotting `fel.dat`, identify the time of the minimum energy, then extract the conformation:

```bash
gmx trjconv -s md_0_100.tpr -f md_0_100.xtc -o min_energy.pdb -dump 520
```

Replace `520` with the appropriate time (in ps) corresponding to the energy minimum.

---

## Summary

This workflow provides a systematic approach to:
- Identify dominant motions in your trajectory via PCA
- Visualize structural variation through projection plots
- Map the free energy landscape based on PC space
- Extract the most stable conformations for further analysis

These methods are widely used for post-simulation analysis in structural biology and drug discovery.
