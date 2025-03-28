---
title: 'How to Validate AlphaFold Structures'
date: 2024-03-28
permalink: /blog/2024/03/how-to-validate-alphafold-structures
excerpt_separator: <!--more-->
toc: true
tags:
  - drug discovery
  - computational chemistry
---

AlphaFold has revolutionized protein structure prediction, but **validating** its output is essential before using the model in downstream applications like docking or molecular dynamics. Here's a step-by-step guide to assess the quality and reliability of AlphaFold-predicted structures.

<!--more-->

# 🧬 How to Validate AlphaFold Structures

## 🔍 1. Check Prediction Confidence

AlphaFold provides two key metrics:

- **pLDDT (per-residue confidence)**:
  - >90: Very high confidence
  - 70–90: Confident
  - <70: Low confidence (often in loops or disordered regions)

- **PAE (Predicted Aligned Error)**:
  - Lower is better.
  - Useful to assess domain-level flexibility or uncertainty in multi-domain proteins.

> 📌 Tip: You can visualize both scores in the AlphaFold DB or with tools like PyMOL or ChimeraX.

---

## 🧫 2. Visual Inspection

Use molecular visualization tools to manually inspect the structure:
- **PyMOL**
- **UCSF Chimera / ChimeraX**
- **Mol\*** (web-based viewer)

Look for:
- Disordered loops
- Steric clashes
- Abnormal bond angles
- Missing secondary structure elements

---

## 🧪 3. Structure Quality Assessment Tools

You can evaluate the stereochemistry and structural quality using these tools:

| Tool | Purpose |
|------|---------|
| **MolProbity** | Clash score, Ramachandran outliers |
| **PROCHECK** | Geometry, dihedrals |
| **WHAT_CHECK** | Stereochemical checks |
| **SwissModel Structure Assessment** | QMEAN Z-score |
| **Verify3D** | Environment profile check |
| **ERRAT** | Non-bonded interaction analysis |

> 📦 Most of these tools accept PDB files and provide web or local versions.

---

## 🔬 4. Compare with Experimental Structures

If homologous structures exist:
- **Align** the AlphaFold model with known crystal/NMR structures using **TM-align**, **PyMOL**, or **Chimera**.
- Compute **RMSD** to assess structural similarity.

---

## ⚙️ 5. Functional and Dynamic Validation

- **Dock known ligands** to check active/binding site geometry.
- Perform **Molecular Dynamics (MD) simulations** to test structure stability over time.
  - Use **GROMACS**, **AMBER**, or **NAMD**.

---

## ✅ Summary

| Validation Step | Why It Matters |
|-----------------|----------------|
| pLDDT/PAE Scores | Confidence in prediction |
| Visual Inspection | Catch obvious errors |
| Geometry Tools | Check stereochemical quality |
| Structural Comparison | Cross-validate with known structures |
| Functional Tests | Assess biological relevance |

---

## 📚 References
- Jumper et al., *Nature* (2021): AlphaFold paper
- AlphaFold Protein Structure Database: https://alphafold.ebi.ac.uk/
- MolProbity: http://molprobity.biochem.duke.edu/

---
