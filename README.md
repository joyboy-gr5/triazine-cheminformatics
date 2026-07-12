#  Triazine Cheminformatics Platform

> **A computational workflow for designing, curating and prioritizing novel 1,2,4-triazine derivatives using an experimentally validated one-pot synthetic strategy.**

---

## Overview

This repository presents a complete **cheminformatics pipeline** for the virtual design of **3-substituted 5,6-dimethyl-1,2,4-triazine derivatives**.

Starting from publicly available hydrazides, the workflow performs rigorous structure curation, virtual reaction enumeration based on an experimentally validated one-pot synthesis, molecular descriptor calculation, similarity analysis, and candidate prioritization.

The ultimate goal is to accelerate the discovery of synthetically accessible triazine derivatives for medicinal chemistry and metal sensing applications.

---

##  Project Highlights

- **High-Throughput Curation:** Standardized and filtered **10,000+ raw PubChem hydrazides** into an electronically pristine aromatic precursor collection.
- **Reaction-Aware Synthesis:** Developed an in silico chemical transformation framework simulating an experimentally proven one-pot multicomponent cyclization.
- **Physicochemical Triaging:** Generated **5,039 virtual triazines** subjected to Lipinski’s Rule of Five boundaries to secure highly drug-like chemical spaces.
- **Coordination-Pocket Screening:** Built structural pincer algorithms tracking donor-atom proximity relative to the triazine core to extract **3,263 metal-coordinating scaffolds.**
- **Topology & Denticity Pruning:** Filtered out strained, non-planar aliphatic systems, mapping **1,901 aromatic candidates.** Used custom SMARTS strings to capture high-affinity tri- and tetradentate pincer patterns, pruning the set to **183 distinct sensors.**
- **Tanimoto Similarity Lead Generation:** Screened the high-affinity subspace against targeted lead chemistries using RDKit Morgan Fingerprints to narrow down **30 highly promising candidates.**
- **Quantum Mechanical Sorting:** Deployed Grimme's semi-empirical tight-binding framework (xTB) to optimize 3D conformers and extract explicit HOMO-LUMO gap topologies, mapping ideal Cu2+ energy decay properties.
---

## Workflow

```text
       [ 10,000+ Raw Precursors ]
                    │
                    ▼ Data Curation & Standardisation
       [ 5,941 Aromatic Hydrazides ]
                    │
                    ▼ Virtual One-Pot Synthesis
       [ 4,999 Total Enumerated Triazines ] (4,939 Mono- / 60 Bis-Triazines)
                    │
                    ▼ Lipinski's Rule of Five Filter
       [ 4,951 Drug-Like Virtual Molecules ]
                    │
                    ▼ Coordination-Pocket Vector Proximity
       [ 3,263 Core Coordination Scaffolds ]
                    │
                    ▼ Geometric Pruning (Aromaticity & Planarity Focus)
       [ 1,901 Rigid Planar Chelators ]
                    │
                    ▼ SMARTS Pattern Identification (Tridentate / Tetradentate Focus)
       [  183 Pincer Ligand Architectures ]
                    │
                    ▼ Fingerprint Similarity Profiling (Morgan / Tanimoto)
       [   30 High-Priority Sensor Candidates ]
                    │
                    ▼ Semi-Empirical Quantum Screening (GFN2-xTB Geometry & Orbitals)
       [ Target Prioritized Lead Portfolio ]

```
---

## Experimental Inspiration

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

## Repository Structure

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
│   ├── Metal_sensing_App.ipynb
│   └── Quantum_Screening.ipynb
│
├── results
│   ├── clean_hydrazides.csv
│   ├── exp_dataset.csv
│   ├── pro_cleaned_triazine.csv
│   ├── virtual_triazine_filtered.csv
│   ├── figures
│   │   ├── PCA.png
│   │   ├── similarity_matrix.png
│   │   ├── top_mols.png
│   │   └── top5_radar.png
│   ├── metal_sensing_app
│   │   ├── Metal_sensing_mol.csv
│   │   └── top_MS_cand.csv
│   └── quantum_screening_data
│       ├── xtb_inputs
│       ├── xtb_outputs
│       ├── xtb_scratch
│       └── xtb_baseline_screeing.csv
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
---

### 5. Metal Sensing Application

To discover an optimal Cu2+ fluorescent "Turn-Off" sensor, the architecture evaluates ligand structures against strict stereochemical and electronic design criteria:
- **The Spatial Constraint:** Pruned out aliphatic and sterically hindered pathways to allow the triazine pincer to transition into the stable, coplanar cis-cis configuration required for strong chelation.
- **Hard-Soft Acid-Base (HSAB) Balancing:** Favored borderline nitrogen-rich donor sets matching the borderline acid profile of Cu2+ to drive high thermodynamic binding affinity.
- **Orbital Quenching Mechanics:** Deployed GFN2-xTB to screen for a targeted optical HOMO-LUMO gap range 1.5 eV - 2.2 eV in the resulting complexes. This narrow gap facilitates non-radiative energy dissipation through direct Ligand-to-Metal Charge Transfer (LMCT) and inner d-d orbital pathways, confirming a highly sensitive optical response. 

---

## Technologies

- Python
- RDKit
- Pandas
- NumPy
- Matplotlib
- Jupyter Notebook
- xtb

---

## Current Status

- Precursor database sanitization and curation engines
- Multi-component reaction simulation compiler
- Multi-tier structural filters (Lipinski, Heteroatom Proximity, Aromaticity)
- Sub-structure SMARTS denticity categorization
- Semi-empirical (GFN2-xTB) optimization and frontier orbital extraction

## Next Steps

- Ab initio Density Functional Theory (DFT) confirmation for Top 5 systems
- Integrated ADMET (SwissADME) ADMET profiling
- Target synthetic accessibility verification score indexing

---

## Future Work

- **Quantum Precision Upgrades:** Migrating the top 5 optimized structures into high-accuracy DFT functional frameworks (e.g., B3LYP or M06-2X using ORCA) to calculate exact absorption spectrum changes.
- **Bench Synthesis Validation:** Subjecting the highest-ranked candidates to our physical, one-pot condensation setup to verify optical Cu2+ sensing properties. 
- **Multi-Metal Vector Tuning:** Shifting structural parameters on the 4'-position to evolve sensors optimized for different transition metal ions.

---

## Author

Ajay K
