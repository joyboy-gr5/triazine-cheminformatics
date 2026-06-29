# Triazine Cheminformatics

A cheminformatics workflow for the design, analysis, and virtual screening of triazine derivatives using Python and RDKit.

## Project Overview

This project aims to combine synthetic organic chemistry with cheminformatics to prioritize novel triazine derivatives before synthesis.

The workflow includes:

- Molecular descriptor calculation
- Structural similarity analysis
- Chemical space visualization
- Compound prioritization
- Virtual screening of aromatic hydrazides from PubChem
- Selection of promising candidates for future triazine synthesis

---

## Project Structure

```
triazine-cheminformatics/
│
├── data/
│   ├── triazine_dataset.csv
│   └── hydrazide_library.csv
│
├── notebooks/
│   ├── 01_descriptor_calculation.ipynb
│   ├── 02_priority_scoring.ipynb
│   ├── 03_similarity_analysis.ipynb
│   ├── 04_Virtual_Library_Generation_and_Screening.ipynb
│
├── results/
│   ├── descriptors.csv
│   ├── filtered_hydrazide.csv
│   ├── similarity_matrix.csv
│   └── triazine_dataset.csv
│
├── src/
└── README.md
```

## Current Progress

- ✅ Descriptor calculation using RDKit
- ✅ Molecular fingerprints (Morgan Fingerprints)
- ✅ Tanimoto similarity analysis
- ✅ Similarity matrix generation
- ✅ PCA visualization
- ✅ Hierarchical clustering
- ✅ Compound ranking
- ✅ PubChem hydrazide library preparation
- ✅ Drug-like filtering of aromatic hydrazides

Current filtered library:

**6855 aromatic hydrazides**

---

## Technologies

- Python
- RDKit
- Pandas
- NumPy
- Matplotlib
- Seaborn
- Scikit-learn

---

## Future Work

- SMARTS-based functional group annotation
- Virtual triazine generation
- Chemical diversity selection
- Novelty analysis
- Machine learning models
- QSPR modelling

---

## Author
Ajay K
