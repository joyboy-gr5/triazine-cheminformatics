#  Triazine Cheminformatics Platform

> **A computational workflow for designing, curating and prioritizing novel 1,2,4-triazine derivatives using an experimentally validated one-pot synthetic strategy.**

---

## Overview

This repository presents a complete **cheminformatics pipeline** for the virtual design of **3-substituted 5,6-dimethyl-1,2,4-triazine derivatives**.

Starting from publicly available hydrazides, the workflow performs rigorous structure curation, virtual reaction enumeration based on an experimentally validated one-pot synthesis, molecular descriptor calculation, similarity analysis, and candidate prioritization.

The ultimate goal is to accelerate the discovery of synthetically accessible triazine derivatives for medicinal chemistry and metal sensing applications.

---

##  Project Highlights

- Curated **10,000+ PubChem hydrazides** into a high-quality aromatic hydrazide library.
- Developed a **reaction-aware virtual synthesis pipeline** to enumerate 1,2,4-triazine derivatives using an experimentally validated one-pot cyclization.
- Generated **5,039 virtual triazines** from curated hydrazides.
- Applied **Lipinski's Rule of Five**, yielding a drug-like virtual library.
- Removed reaction-incompatible molecules containing additional hydrazine/acyl hydrazide functionalities.
- Annotated the library using **PAINS, Brenk, NIH, and ZINC structural alerts**.
- Evaluated **Synthetic Accessibility (SA Score)** for every compound.
- Final curated library contains **4,982 reaction-compatible, drug-like triazines**.
- Identified **119 virtual analogues** containing the same **2-hydroxy substitution pattern** as the experimentally synthesized lead compounds.

---

# Workflow

```text
10,000+ PubChem Hydrazides
            │
            ▼
Data Curation
• Remove invalid structures
• Remove salts and charged species
• Keep aromatic hydrazides
            │
            ▼
Reaction Enumeration
Virtual one-pot synthesis
Hydrazide → 1,2,4-Triazine
            │
            ▼
Drug-likeness Filtering
• Lipinski Rule of Five
            │
            ▼
Reaction-specific Filtering
• Remove residual hydrazine/acyl hydrazide motifs
            │
            ▼
Structural Quality Assessment
• PAINS
• Brenk
• NIH
• ZINC
            │
            ▼
Synthetic Feasibility
• Synthetic Accessibility (SA Score)
            │
            ▼
Final Curated Library
4,951 Virtual Triazines
```

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
3-Substituted-5,6-Dimethyl-1,2,4-Triazine
```

---

# Repository Structure

```
TRIAZINE_QSPR_PROJECT
│
├── data
│   ├── hydrazide_library.csv
│   └── triazine_data.csv
│
├── notebooks
│   ├── 01_experimental_compounds.ipynb
│   ├── 02_hydrazide_cleaning.ipynb
│   ├── 03_virtual_reaction.ipynb
│   ├── 04_library_triage.ipynb
│   ├── 05_SA_Scoring.ipynb
│   └── a_metal_sensing_analysis.ipynb
│
├── results
│   ├── clean_hydrazides.csv
│   ├── exp_dataset.csv
│   ├── pro_cleaned_triazine.csv
│   ├── virtual_triazine_filtered.csv
│   └── figures
│       ├── PCA.png
│       ├── QED.png
│       ├── SA_score.png
│       └── similarity_matrix.png
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

**5,941 aromatic hydrazides**

---

### 2. Virtual Reaction Enumeration

Each curated hydrazide was computationally transformed into the corresponding **3-substituted 5,6-dimethyl-1,2,4-triazine** based on the experimentally validated synthetic reaction.

Generated library:

- **4,939 mono-triazines**
- **60 bis-triazines**

Total:

**4,999 virtual triazine derivatives**

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

### 6. Scaffold Investigation for metal sensing

A focused search identified

**120 virtual triazines**

sharing the same **2-hydroxy substitution pattern** as the experimentally synthesized lead scaffold.

---

### 7. PAIN, Brenk, NIH, Zinc Annotation

A dedicated notebook investigates the quantity of compounds that give these alerts. 

---

# Technologies

- Python
- RDKit
- Pandas
- NumPy
- Matplotlib
- Jupyter Notebook

---

## Current Status

- Virtual reaction engine completed

- Library curation completed

- Drug-likeness filtering completed

- Structural alert annotation completed

- Synthetic accessibility analysis completed

# Next Steps

- Quantitative Estimate of Drug-likeness (QED)
- Chemical diversity analysis
- Bemis–Murcko scaffold analysis
- Similarity search against synthesized compounds
- Multi-parameter compound prioritization
---

# Future Work

- QSPR model development
- ADMET prediction
- Synthetic accessibility scoring
- Experimental validation of prioritized candidates
- DFT analysis of top molecules

---

# Author

Ajay K
