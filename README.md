#  Triazine Cheminformatics Platform

> **A computational workflow for designing, curating and prioritizing novel 1,2,4-triazine derivatives using an experimentally validated one-pot synthetic strategy.**

---

## Overview

This repository presents a complete **cheminformatics pipeline** for the virtual design of **3-substituted 5,6-dimethyl-1,2,4-triazine derivatives**.

Starting from publicly available hydrazides, the workflow performs rigorous structure curation, virtual reaction enumeration based on an experimentally validated one-pot synthesis, molecular descriptor calculation, similarity analysis, and candidate prioritization.

The ultimate goal is to accelerate the discovery of synthetically accessible triazine derivatives for medicinal chemistry and metal sensing applications.

---

# Workflow

<p align="center">
<img src="figures/workflow.png" width="900">
</p>

---

# Experimental Inspiration

The virtual reaction implemented in this repository is based on our experimentally validated one-pot synthesis.

```
Hydrazide
      +
2,3-Butanedione
      +
NH4OAc
          │
      DMF, 110 °C
          ▼
3-Substituted
5,6-Dimethyl-1,2,4-Triazine
```

---

# Project Statistics

| Stage | Molecules |
|--------|----------:|
| Raw PubChem hydrazides | ~10,000 |
| Curated aromatic hydrazides | **5,971** |
| Mono-triazine derivatives | **4,981** |
| Bis-triazine derivatives | **58** |
| Total virtual triazine library | **5,039** |

---

# Repository Structure

```
TRIAZINE_QSPR_PROJECT
│
├── data
│   ├── clean_hydrazides.csv
│   ├── hydrazide_library.csv
│   ├── triazine_data.csv
│   └── virtual_triazines_filtered.csv
│
├── notebooks
│   ├── 01_descriptor_calculation.ipynb
│   ├── 02_priority_scoring.ipynb
│   ├── 03_similarity_analysis.ipynb
│   ├── 04_hydrazide_cleaning.ipynb
│   ├── 05_virtual_reaction.ipynb
│   └── 06_metal_sensing_analysis.ipynb
│
├── results
│   ├── descriptors.csv
│   ├── similarity_matrix.csv
│   ├── triazine_dataset.csv
│   └── virtual_triazine_library.csv
│
├── figures
│   ├── dendrogram.png
│   ├── pca_plot.png
│   └── similarity_heat_map.png
│
└── src
    └── read_data.py
```

---

# Computational Pipeline

### 1. Hydrazide Library Construction

- PubChem data collection
- Invalid SMILES removal
- Hydrazide verification
- Salt and mixture removal
- Aromatic compound filtering

Final curated library:

**5,971 aromatic hydrazides**

---

### 2. Virtual Reaction Enumeration

Each curated hydrazide was computationally transformed into the corresponding **3-substituted 5,6-dimethyl-1,2,4-triazine** based on the experimentally validated synthetic reaction.

Generated library:

- **4,981 mono-triazines**
- **58 bis-triazines**

Total:

**5,039 virtual triazine derivatives**

---

### 3. Descriptor Calculation

Molecular descriptors were calculated using RDKit, providing structural and physicochemical properties for downstream analysis.

---

### 4. Candidate Prioritization

Virtual compounds were filtered using:

- Lipinski Rule of Five
- Molecular descriptors
- Structural similarity

---

### 5. Similarity Analysis

The virtual library was analyzed using molecular fingerprints and visualized through:

- Similarity matrix
- Hierarchical clustering
- PCA

---

### 6. Scaffold Investigation

A focused search identified

**119 virtual triazines**

sharing the same **2-hydroxy substitution pattern** as the experimentally synthesized lead scaffold.

---

### 7. Metal Sensing Analysis

A dedicated notebook investigates the applicability of selected virtual triazines toward metal sensing studies.

---

# Technologies

- Python
- RDKit
- Pandas
- NumPy
- Matplotlib
- Jupyter Notebook

---

# Current Status

✅ Hydrazide curation

✅ Virtual reaction generation

✅ Descriptor calculation

✅ Lipinski filtering

✅ Similarity analysis

✅ Priority scoring

✅ Metal sensing analysis

---

# Future Work

- QSPR model development
- ADMET prediction
- Molecular docking
- Synthetic accessibility scoring
- Biological activity prediction
- Experimental validation of prioritized candidates

---

# Author

Ajay K
