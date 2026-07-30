#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd

hydrazide_df = pd.read_csv("data/hydrazide_library.csv")
len(hydrazide_df)


# In[2]:


'''removing invalid smiles from the source'''

from rdkit import Chem

hydrazide_df["mol"] = hydrazide_df["SMILES"].apply(Chem.MolFromSmiles)

invalid_count = hydrazide_df["mol"].isna().sum()
print(f"Invalid SMILES: {invalid_count}")

invalid_smiles = hydrazide_df[hydrazide_df["mol"].isna()]

hydrazide_df = hydrazide_df[
    hydrazide_df["mol"].notna()
].reset_index(drop=True)

len(hydrazide_df)


# In[3]:


'''making sure all molecules have valid hydrazide moeity'''

hydrazide_pattern = Chem.MolFromSmarts("[CX3](=O)[NH][NH2]")

true_hydrazide_df = hydrazide_df[
    hydrazide_df["mol"].apply(
        lambda mol: mol.HasSubstructMatch(hydrazide_pattern)
    )
].reset_index(drop=True)

len(true_hydrazide_df)


# In[4]:


'''removing salts from the library'''

neutral_hydrazide_df = true_hydrazide_df[
    true_hydrazide_df["mol"].apply(
        lambda mol: Chem.GetFormalCharge(mol) == 0
    )
].reset_index(drop=True)

len(neutral_hydrazide_df)


# In[5]:


'''calculating no. of aromatic rings present in the molecule'''

from rdkit.Chem import rdMolDescriptors

neutral_hydrazide_df["Aromatic_Ring_Count"] = hydrazide_df["mol"].apply(
    rdMolDescriptors.CalcNumAromaticRings
)


# In[6]:


'''Applying Lipinski's rule'''

filtered_hydrazide_df = neutral_hydrazide_df[
    (neutral_hydrazide_df["Molecular_Weight"] < 500) &
    (neutral_hydrazide_df["XLogP"] < 5) &
    (neutral_hydrazide_df["Polar_Area"] < 140) &
    (neutral_hydrazide_df["Rotatable_Bond_Count"] < 10) &
    (neutral_hydrazide_df["Aromatic_Ring_Count"] >= 1)
].reset_index(drop=True)

len(filtered_hydrazide_df)


# In[7]:


'''catagorizing molecules based on the func. groups, heterocycles, halogens'''

functional_groups = {
    # Nitrogen heterocycles
    "Pyridine": Chem.MolFromSmarts("n1ccccc1"),
    "Pyrrole": Chem.MolFromSmarts("[nH]1cccc1"),
    "Imidazole": Chem.MolFromSmarts("c1ncc[nH]1"),
    "Pyrazole": Chem.MolFromSmarts("c1n[nH]cc1"),
    "1,2,3-Triazole": Chem.MolFromSmarts("n1nncc1"),
    "1,2,4-Triazole": Chem.MolFromSmarts("n1ncnc1"),
    "1,2,3,4-Tetrazole": Chem.MolFromSmarts("c1nnnn1"),
    "1-H-Tetrazole": Chem.MolFromSmarts("[nH]1nnnc1"),
    "Quinoline": Chem.MolFromSmarts("c1ccc2ncccc2c1"),
    "Isoquinoline": Chem.MolFromSmarts("c1ccc2ccnc2c1"),

    # Oxygen/Sulfur heterocycles
    "Thiophene": Chem.MolFromSmarts("c1ccsc1"),
    "Furan": Chem.MolFromSmarts("c1ccoc1"),
    "Benzofuran": Chem.MolFromSmarts("c1ccc2occc2c1"),

    # Functional groups
    "Phenol": Chem.MolFromSmarts("c[OH]"),
    "Methoxy": Chem.MolFromSmarts("[OX2][CH3]"),
    "Nitro": Chem.MolFromSmarts("[N+](=O)[O-]"),
    "Cyano": Chem.MolFromSmarts("C#N"),

    # Halogens
    "Fluoro": Chem.MolFromSmarts("[F]"),
    "Chloro": Chem.MolFromSmarts("[Cl]"),
    "Bromo": Chem.MolFromSmarts("[Br]"),
    "Iodo": Chem.MolFromSmarts("[I]"),

    # Carbonyl derivatives
    "Carboxylic_Acid": Chem.MolFromSmarts("C(=O)[OH]"),
    "Ester": Chem.MolFromSmarts("C(=O)O[#6]"),
    "Aldehyde": Chem.MolFromSmarts("[CX3H](=O)"),
    "Ketone": Chem.MolFromSmarts("[#6][CX3](=O)[#6]"),

    # Sulfur groups
    "Sulfone": Chem.MolFromSmarts("S(=O)(=O)"),
    "Sulfonamide": Chem.MolFromSmarts("S(=O)(=O)N"),
    "Trifluoromethyl": Chem.MolFromSmarts("C(F)(F)F")
}

for name, pattern in functional_groups.items():
    filtered_hydrazide_df[name] = filtered_hydrazide_df["mol"].apply(
        lambda mol: mol.HasSubstructMatch(pattern)
    )

fg_counts = filtered_hydrazide_df[
    functional_groups.keys()
].sum().sort_values(ascending=False)



# In[8]:


'''visualising the func. group distribution in our library'''

import matplotlib.pyplot as plt

fg_counts.plot.bar(figsize=(12,5))

plt.ylabel("Number of molecules")
plt.title("Functional Group Distribution in Curated Hydrazide Library")
plt.xticks(rotation=45)

plt.tight_layout()
plt.show()


# In[9]:


'''identifying hydrazide counts'''

hydrazide = Chem.MolFromSmarts("[CX3](=O)[NH][NH2]")
filtered_hydrazide_df["Hydrazide_Count"] = filtered_hydrazide_df["mol"].apply(
    lambda m: len(m.GetSubstructMatches(hydrazide))
)

filtered_hydrazide_df["Hydrazide_Count"].value_counts().sort_index()


# In[10]:


from rdkit.Chem import Draw

mols_2 = filtered_hydrazide_df[
    filtered_hydrazide_df["Hydrazide_Count"] == 1
]["mol"].tail(10).tolist()

legends_2 = filtered_hydrazide_df[
    filtered_hydrazide_df["Hydrazide_Count"] == 1
]["Name"].astype(str).tail(21).tolist()

Draw.MolsToGridImage(
    mols_2,
    legends=legends_2,
    molsPerRow=5,
    subImgSize=(300,300)
)


# In[11]:


'''Identifying Aromaticity of the molecules'''

def has_aromatic_ring(mol):
    return any(atom.GetIsAromatic() for atom in mol.GetAtoms())

filtered_hydrazide_df["Aromatic"] = filtered_hydrazide_df["mol"].apply(has_aromatic_ring)

filtered_hydrazide_df["Aromatic"].value_counts()


# In[12]:


'''Removing aliphatic molecules from the library'''

def has_aromatic_ring(mol):
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            return True
    return False

filtered_hydrazide_df = filtered_hydrazide_df[
    filtered_hydrazide_df["mol"].apply(has_aromatic_ring)
].copy()

len(filtered_hydrazide_df)


# In[13]:


filtered_hydrazide_df.to_csv(
    "results/hydrazides.csv",
    index=False
)

#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd

hydrazide_df = pd.read_csv("results/hydrazides.csv")


# In[5]:


'''reaction enumeration, converting all hydrazides into 5,6-dimethyl-1,2,4-triazines using Smart pattern'''

from rdkit import Chem
from rdkit.Chem import AllChem

hydrazide = Chem.MolFromSmarts('[CX3:1](=[OX1])[NX3;H1:2][NX3;H2:3]')

rxn = AllChem.ReactionFromSmarts(
    '[CX3:1](=[OX1])[NX3;H1:2][NX3;H2:3]>>[c:1]1[n:2][n:3]c(C)c(C)n1'
)

def to_triazine(mol):
    """Cyclize every hydrazide group present (handles mono- and di-hydrazides)."""
    working = mol
    for _ in range(5):  # safety cap
        if not working.HasSubstructMatch(hydrazide):
            break
        prods = rxn.RunReactants((working,))
        if not prods:
            break
        new_mol = prods[0][0]
        try:
            Chem.SanitizeMol(new_mol)
        except Exception:
            break
        working = new_mol
    return working

out = []
for smi in hydrazide_df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol is None or not mol.HasSubstructMatch(hydrazide):
        out.append(None)
        continue
    prod = to_triazine(mol)
    out.append(Chem.MolToSmiles(prod))

hydrazide_df['Triazine_SMILES'] = out

print(hydrazide_df['Triazine_SMILES'] )


# In[6]:


'''visual check of random mono-hydrazides and their corresponding mono-triazines'''

from rdkit.Chem import Draw

sample_df = hydrazide_df.dropna(subset=['Triazine_SMILES']).sample(6, random_state=42)

mols = []
legends = []

for _, row in sample_df.iterrows():
    hydrazide_mol = Chem.MolFromSmiles(row['SMILES'])
    triazine_mol  = Chem.MolFromSmiles(row['Triazine_SMILES'])
    mols.extend([hydrazide_mol, triazine_mol])
    legends.extend(['hydrazide', 'triazine'])

img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(300, 300),
    legends=legends,
    useSVG=False
)
img


# In[7]:


'''random dihydrazides and bis-triazines'''

dihydrazide = Chem.MolFromSmarts('[CX3](=[OX1])[NX3;H1][NX3;H2]')

def count_hydrazide_groups(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return 0
    return len(mol.GetSubstructMatches(dihydrazide))

hydrazide_df['n_hydrazide_groups'] = hydrazide_df['SMILES'].apply(count_hydrazide_groups)
di_sample = hydrazide_df[hydrazide_df['n_hydrazide_groups'] == 2].dropna(subset=['Triazine_SMILES']).sample(6, random_state=1)

mols, legends = [], []
for _, row in di_sample.iterrows():
    mols.append(Chem.MolFromSmiles(row['SMILES']))
    mols.append(Chem.MolFromSmiles(row['Triazine_SMILES']))
    legends += ['dihydrazide', 'bis-triazine']

Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(320,320), legends=legends)


# In[8]:


'''separate dataframe for generated triazines'''

triazine_df = hydrazide_df.loc[hydrazide_df['Triazine_SMILES'].notna(), ['Triazine_SMILES']].copy()
triazine_df = triazine_df.rename(columns={'Triazine_SMILES': 'SMILES'})
triazine_df = triazine_df.reset_index(drop=True)
triazine_df.insert(0, 'Triazine_ID', ['TRZ' + str(i+1).zfill(5) for i in range(len(triazine_df))])

if 'Hydrazide_ID' in hydrazide_df.columns:
    triazine_df['Source_Hydrazide_ID'] = hydrazide_df.loc[hydrazide_df['Triazine_SMILES'].notna(), 'Hydrazide_ID'].reset_index(drop=True)
else:
    triazine_df['Source_Hydrazide_SMILES'] = hydrazide_df.loc[hydrazide_df['Triazine_SMILES'].notna(), 'SMILES'].reset_index(drop=True)

triazine_df.head()


# In[9]:


'''generating descriptors for filtration'''

from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

descriptor_list = []

for index, row in triazine_df.iterrows():
    Triazine_ID = row["Triazine_ID"]
    smiles = row["SMILES"]

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("Invalid SMILES:", Triazine_ID, smiles)
        continue

    descriptor_data = {
        "Triazine_ID": Triazine_ID,
        "SMILES": smiles,
        "MolWt": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol),
        "AromaticRings": Lipinski.NumAromaticRings(mol),
    }

    descriptor_list.append(descriptor_data)

triazine_descriptor_df = pd.DataFrame(descriptor_list)
triazine_descriptor_df.head(10)


# In[10]:


'''counting number of 1,2,4-triazine ring present in a molecule'''

triazine = Chem.MolFromSmarts('c1(C)nnc(nc1C)')

def count_triazine_rings(smiles, query):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    matches = mol.GetSubstructMatches(query, uniquify=True)
    return len(matches)

# apply to your dataframe
triazine_descriptor_df['n_triazine_rings'] = triazine_descriptor_df['SMILES'].apply(
    lambda s: count_triazine_rings(s, triazine)
)

triazine_descriptor_df[['SMILES', 'n_triazine_rings']].head(10)


# In[11]:


'''filtration based on lipinski rule'''

filtered_triazine_df = triazine_descriptor_df[
    (triazine_descriptor_df["MolWt"] < 500) &
    (triazine_descriptor_df["LogP"] < 5) &
    (triazine_descriptor_df["TPSA"] < 140) &
    (triazine_descriptor_df["RotatableBonds"] < 10)
].reset_index(drop=True)

print(len(filtered_triazine_df))


# In[12]:


ring_counts = filtered_triazine_df['n_triazine_rings'].value_counts().sort_index()
print(ring_counts)


# In[13]:


'''Molecules containing additional hydrazine or acyl hydrazide moieties outside the reacting hydrazide center were removed because these
functional groups are expected to undergo competing condensation with 2,3-butanedione, generating stable hydrazone intermediates that
cannot undergo the desired triazine cyclization under the optimized one-pot reaction conditions.'''


hydrazine = Chem.MolFromSmarts("[NX3][NX3]")
acyl_hydrazide = Chem.MolFromSmarts("[CX3](=O)[NX3][NX3]")

filtered_triazine_df["mol"] = filtered_triazine_df["SMILES"].apply(
    Chem.MolFromSmiles
)

filtered_triazine_df["Has_Hydrazine"] = (
    filtered_triazine_df["mol"]
    .apply(lambda m: m.HasSubstructMatch(hydrazine))
)

filtered_triazine_df["Has_Acyl_Hydrazide"] = (
    filtered_triazine_df["mol"]
    .apply(lambda m: m.HasSubstructMatch(acyl_hydrazide))
)

filtered_triazine_df = filtered_triazine_df[
    ~(filtered_triazine_df["Has_Hydrazine"] |
      filtered_triazine_df["Has_Acyl_Hydrazide"])
].copy()

filtered_triazine_df.to_csv('results/triazines.csv', index=False)


# In[14]:


len(filtered_triazine_df)

#!/usr/bin/env python
# coding: utf-8

# In[32]:


import pandas as pd

triazine_df = pd.read_csv("results/triazines.csv")


# In[33]:


'''Let's find molecules having donar atoms at pos1, pos2 and pos3. '5,6-dimethyl-1,2,4-triazine-3-pos1-pos2-pos3'''

from rdkit import Chem

# Donor at position 1 (directly on C3)
pos1 = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#7,#8,#16]')

# Donor at position 2 (1 carbon away)
pos2 = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#7,#8,#16]')

# Donor at position 3 (2 carbons away)
pos3 = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#6]~[#7,#8,#16]')



triazine_df['has_donor_pos1'] = triazine_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos1)
    if Chem.MolFromSmiles(smi) else False
)

triazine_df['has_donor_pos2'] = triazine_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos2)
    if Chem.MolFromSmiles(smi) else False
)

triazine_df['has_donor_pos3'] = triazine_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos3)
    if Chem.MolFromSmiles(smi) else False
)

donor_df = triazine_df[
    (triazine_df['has_donor_pos1']) |
    (triazine_df['has_donor_pos2']) |
    (triazine_df['has_donor_pos3'])
].reset_index(drop=True)

print(f"Before: {len(triazine_df)}")
print(f"After:  {len(donor_df)}")
print(f"Removed: {len(triazine_df) - len(donor_df)}")
print(f"Pos1 donors: {triazine_df['has_donor_pos1'].sum()}")
print(f"Pos2 donors: {triazine_df['has_donor_pos2'].sum()}")
print(f"Pos3 donors: {triazine_df['has_donor_pos3'].sum()}")



# In[35]:


both = triazine_df[
    (triazine_df['has_donor_pos2']) &
    (triazine_df['has_donor_pos3'])
]
print(f"Molecules with donors at BOTH pos2 and pos3: {len(both)}")


# In[36]:


'''keeping molecules, that are having donor atoms at any one position'''

valid_df = triazine_df[
    (triazine_df['has_donor_pos2']) |
    (triazine_df['has_donor_pos3'])
].reset_index(drop=True)

print(f"Valid chelating candidates: {len(valid_df)}")


# In[37]:


'''Obviously, Bis-triazine will be the top candidate for metal sensing in this library. So, let's find in both Bis- and Mono-'''

bis_df = valid_df[valid_df['n_triazine_rings'] == 2].reset_index(drop=True)

mono_df = valid_df[valid_df['n_triazine_rings'] == 1].reset_index(drop=True)

print(f"Bis-triazines:  {len(bis_df)}")
print(f"Mono-triazines: {len(mono_df)}")


# In[38]:


'''segregating molecules donor atoms position, either it is at pos2 or pos3 or both'''

only_pos2 = mono_df[
    (mono_df['has_donor_pos2']) & (~mono_df['has_donor_pos3'])
]

only_pos3 = mono_df[
    (mono_df['has_donor_pos3']) & (~mono_df['has_donor_pos2'])
] 

both_pos2_pos3 = mono_df[
    (mono_df['has_donor_pos2']) & (mono_df['has_donor_pos3'])
] 

print(f"Only pos2:      {len(only_pos2)}")
print(f"Only pos3:      {len(only_pos3)}")
print(f"Both pos2+pos3: {len(both_pos2_pos3)}")
print(f"Total:          {len(only_pos2) + len(only_pos3) + len(both_pos2_pos3)}")


# In[39]:


'''Bridged and non-bridged compounds can have different geometry, even if their donor position matches. So, let's catagorize it.'''

pos2_bridged    = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X4;!R]~[#7,#8,#16]')
a = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X3]~[#7,#8,#16]')
b = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X4;R]~[#7,#8,#16]')
pos3_bridged    = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X4;!R]~[#6]~[#7,#8,#16]')
c = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X3]~[#6]~[#7,#8,#16]')
d = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6;X4;R]~[#6]~[#7,#8,#16]')


mono_df['has_pos2_bridged'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos2_bridged)
    if Chem.MolFromSmiles(smi) else False
)

mono_df['has_pos3_bridged'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos3_bridged)
    if Chem.MolFromSmiles(smi) else False
)

def check_bridged(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    return mol.HasSubstructMatch(a) or mol.HasSubstructMatch(b)

mono_df['has_pos2_notbridged'] = mono_df['SMILES'].apply(check_bridged)

def check_bridged(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    return mol.HasSubstructMatch(c) or mol.HasSubstructMatch(d)

mono_df['has_pos3_notbridged'] = mono_df['SMILES'].apply(check_bridged)

bridge_df = mono_df[
    (mono_df['has_pos2_bridged']) |
    (mono_df['has_pos2_notbridged']) |
    (mono_df['has_pos3_bridged']) |
    (mono_df['has_pos3_notbridged'])
].reset_index(drop=True)

print(f"Pos2 bridged: {mono_df['has_pos2_bridged'].sum()}")
print(f"Pos2 notbridged: {mono_df['has_pos2_notbridged'].sum()}")
print(f"Pos3 bridged: {mono_df['has_pos3_bridged'].sum()}")
print(f"Pos3 notbridged: {mono_df['has_pos3_notbridged'].sum()}")
len(bridge_df)


# In[40]:


'''narrowing further to know.. which position has how many bridged/non-bridged compounds'''

only_pos2_bridged = mono_df[
    (mono_df['has_pos2_bridged']) & 
    (~mono_df['has_pos3_bridged'])&
    (~mono_df['has_pos3_notbridged'])
]

only_pos2_notbridged = mono_df[
    (mono_df['has_pos2_notbridged']) & 
    (~mono_df['has_pos3_notbridged'])&
    (~mono_df['has_pos3_bridged'])
]  
only_pos3_bridged    = mono_df[
    (mono_df['has_pos3_bridged']) & 
    (~mono_df['has_pos2_bridged'])&
    (~mono_df['has_pos2_notbridged'])
]  
only_pos3_notbridged = mono_df[
    (mono_df['has_pos3_notbridged']) & 
    (~mono_df['has_pos2_bridged']) & 
    (~mono_df['has_pos2_notbridged'])
]  

both_bridged = mono_df[
    (mono_df['has_pos2_bridged'])&
    (mono_df['has_pos3_bridged'])
]

both_notbridged = mono_df[
    (mono_df['has_pos2_notbridged'])&
    (mono_df['has_pos3_notbridged'])
]

print(f"Only pos2 bridged:    {len(only_pos2_bridged)}")
print(f"Only pos2 notbridged: {len(only_pos2_notbridged)}")
print(f"Only pos3 bridged:    {len(only_pos3_bridged)}")
print(f"Only pos3 notbridged: {len(only_pos3_notbridged)}")
print(f"Both Bridged:         {len(both_bridged)}")
print(f"Both Not Bridged:     {len(both_notbridged)}")
print(f"Total:                {len(only_pos2_bridged)+len(only_pos2_notbridged)+len(only_pos3_bridged)+len(only_pos3_notbridged)+len(both_bridged)+len(both_notbridged)}")


# In[41]:


'''Let's store this data in a separate dataframe'''

base_cols = ['Triazine_ID', 'SMILES', 'MolWt', 'LogP', 'TPSA', 'HBD', 'HBA',
             'RotatableBonds', 'AromaticRings', 'n_triazine_rings', 'mol']

groups = {
    'only_pos2_bridged': only_pos2_bridged,
    'only_pos2_notbridged': only_pos2_notbridged,
    'only_pos3_bridged': only_pos3_bridged,
    'only_pos3_notbridged': only_pos3_notbridged,
    'both_bridged': both_bridged,
    'both_notbridged': both_notbridged,
}

pieces = []
for name, df in groups.items():
    tmp = df[base_cols].copy()
    tmp['Donor_pos'] = name
    pieces.append(tmp)

geometry_df = pd.concat(pieces, ignore_index=True)

print(geometry_df['Donor_pos'].value_counts())
print("Total:", len(geometry_df)) 


# In[42]:


'''Knowing, whether the donor is in the ring or not can help us know the rigidity of the molecule upon chelation'''

pos2_inring    = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#7,#8,#16;R]')
pos2_notinring = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#7,#8,#16;!R]')
pos3_inring    = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#6]~[#7,#8,#16;R]')
pos3_notinring = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#6]~[#7,#8,#16;!R]')


mono_df['pos2_inring'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos2_inring)
    if Chem.MolFromSmiles(smi) else False
)

mono_df['pos2_notinring'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos2_notinring)
    if Chem.MolFromSmiles(smi) else False
)

mono_df['pos3_inring'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos3_inring)
    if Chem.MolFromSmiles(smi) else False
)

mono_df['pos3_notinring'] = mono_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos3_notinring)
    if Chem.MolFromSmiles(smi) else False
)


def assign_ring_status(row):
    p2in  = row['pos2_inring']
    p2out = row['pos2_notinring']
    p3in  = row['pos3_inring']
    p3out = row['pos3_notinring']

    # Build descriptive label showing ALL donor information
    parts = []

    if p2in and p2out: parts.append('pos2_both')
    elif p2in:         parts.append('pos2_inring')
    elif p2out:        parts.append('pos2_notinring')

    if p3in and p3out: parts.append('pos3_both')
    elif p3in:         parts.append('pos3_inring')
    elif p3out:        parts.append('pos3_notinring')

    if parts: return '+'.join(parts)
    return 'no_donor'

mono_df['ring_status'] = mono_df.apply(assign_ring_status, axis=1)
print(mono_df['ring_status'].value_counts())


# In[43]:


'''clean dataset for metal sensing'''

fina_df = pd.merge(geometry_df, mono_df)

final_df = fina_df.drop(columns=[ 'Has_Hydrazine', 'Has_Acyl_Hydrazide'])


final_df.head()


# In[44]:


'''removing invalid combos, eg; donor atom at pos3 only and the compound is not bridged molecules can't act as bidentate ligand as their geometry won't allow it.'''

# Invalid combinations using Donor_pos + ring_status together
invalid1 = (final_df['Donor_pos'] == 'only_pos3_notbridged') & (final_df['pos3_inring'])
invalid2 = (final_df['Donor_pos'] == 'only_pos2_bridged') & (final_df['pos2_inring'])

valid_final_df = final_df[
    ~invalid1 & ~invalid2
].reset_index(drop=True)

print(f"Before: {len(final_df)}")
print(f"After:  {len(valid_final_df)}")
print(f"Removed: {len(final_df) - len(valid_final_df)}")


# In[45]:


'''molecules having more no. of rings with donor atoms may have higher spectroscopic properties and also can help to find multidentate ligands'''

from rdkit.Chem import rdMolDescriptors

valid_final_df['total_rings'] = valid_final_df['SMILES'].apply(
    lambda smi: rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(smi))
    if Chem.MolFromSmiles(smi) else 0
)

print(valid_final_df['total_rings'].value_counts().sort_index())


# In[46]:


'''Let's remove non-aromatic rings, As aromatic rings can give distinct fluorescene change upon metal chelation'''

fully_aromatic_df = valid_final_df[
    valid_final_df['total_rings'] == valid_final_df['AromaticRings']
].reset_index(drop=True)

print(f"Fully aromatic: {len(fully_aromatic_df)}")
print(fully_aromatic_df['total_rings'].value_counts().sort_index())


# In[47]:


'''removed 5,6 ring molecules after inspected manually...their geometry is not efficient for metal sensing'''

aromati_filtered_df = fully_aromatic_df[
    fully_aromatic_df['total_rings'] <= 4
].reset_index(drop=True)

print(f"Before: {len(fully_aromatic_df)}")
print(f"After:  {len(aromati_filtered_df)}")


# In[48]:


aromati_filtered_df['extra_donors'] = aromati_filtered_df['HBA'] - 3

print(aromati_filtered_df['extra_donors'].value_counts().sort_index())


# In[49]:


'''remove molecules that are having only one donor atoms apart from triazne N's.'''

aromatic_filtered_df = aromati_filtered_df[
    aromati_filtered_df['extra_donors'] >= 2
].reset_index(drop=True)

print(f"After:  {len(aromatic_filtered_df)}")
print(aromatic_filtered_df['extra_donors'].value_counts().sort_index())


# In[50]:


'''removing molecules that are having plain phenyl rings, rings with donor atoms can be more efficient for metal binding'''

plain_phenyl = Chem.MolFromSmarts('c1ccccc1')

aromatic_filtered_df['has_plain_phenyl'] = aromatic_filtered_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pos3_notinring)
    if Chem.MolFromSmiles(smi) else False
)

no_phenyl_df = aromatic_filtered_df[
    ~aromatic_filtered_df['has_plain_phenyl']
].reset_index(drop=True)

print(f"Before: {len(aromatic_filtered_df)}")
print(f"After:  {len(no_phenyl_df)}")


# In[51]:


direct = no_phenyl_df[
    (no_phenyl_df['total_rings'] >= 3) &
    (no_phenyl_df['pos2_inring']) &
    (no_phenyl_df['extra_donors'] >= 2)
]

bridged = no_phenyl_df[
    (no_phenyl_df['total_rings'] >= 3) &
    (no_phenyl_df['pos3_inring']) &
    (~no_phenyl_df['pos2_inring']) &
    (no_phenyl_df['extra_donors'] >= 2)
]

print(f"Direct connection: {len(direct)}")
print(f"CH2 bridged:       {len(bridged)}")


# In[52]:


'''Let's find molecules, that are having potential to be tridentate ligand'''

pattern = Chem.MolFromSmarts('c1(C)nnc(nc1C)~[#6]~[#7,#8,#16]~[#6]~[#6]~[#7,#8,#16]')

no_phenyl_df['tridentate'] = no_phenyl_df['SMILES'].apply(
    lambda smi: Chem.MolFromSmiles(smi).HasSubstructMatch(pattern) if Chem.MolFromSmiles(smi) else False
)

new_df = no_phenyl_df[no_phenyl_df['tridentate'] == True].reset_index(drop=True)

print(f"Total: {no_phenyl_df['tridentate'].sum()}")


# In[53]:


'''Through this process, few molecules have been found that can be potential tridentate ligand'''

reference_ids = [
    'TRZ02396', 'TRZ02390', 'TRZ00399',
    'TRZ00375', 'TRZ05498', 'TRZ01544',
    'TRZ02392', 'TRZ00628', 'TRZ00117', 
    'TRZ01080', 'TRZ00112', 'TRZ05734'
]

reference_df = valid_final_df[
    valid_final_df['Triazine_ID'].isin(reference_ids)
].reset_index(drop=True)

print(f"Reference compounds: {len(reference_df)}")
reference_df.to_csv('data/reference_compounds.csv', index=False)


# In[54]:


'''let's generate fingerprints of the reference and our libary molecules'''

from rdkit.Chem import rdFingerprintGenerator

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

ref_fps = [morgan_gen.GetFingerprint(Chem.MolFromSmiles(smi)) 
           for smi in reference_df['SMILES']]

lib_fps = [morgan_gen.GetFingerprint(Chem.MolFromSmiles(smi)) 
           for smi in new_df['SMILES']]

print(f"Reference fingerprints: {len(ref_fps)}")
print(f"Library fingerprints:   {len(lib_fps)}")


# In[55]:


'''Using TanimotoSimilarity, we can find molecules that are similar to reference compounds'''

from rdkit import DataStructs

def max_similarity_to_references(lib_fp):
    similarities = DataStructs.BulkTanimotoSimilarity(lib_fp, ref_fps)
    return max(similarities)

new_df['max_similarity'] = [
    max_similarity_to_references(fp) for fp in lib_fps
]

print(new_df['max_similarity'].describe())
print(f"\nDistribution:")
print(f"Similarity > 0.8: {(new_df['max_similarity'] > 0.8).sum()}")
print(f"Similarity > 0.6: {(new_df['max_similarity'] > 0.6).sum()}")
print(f"Similarity > 0.4: {(new_df['max_similarity'] > 0.4).sum()}")


# In[56]:


'''visualize molecules that are having maximum similarity with our references'''
from rdkit.Chem import Draw

similar_df = new_df[
    (new_df['max_similarity'] >= 0.6) &
    (new_df['max_similarity'] <= 1.0)
].reset_index(drop=True)

print(f"Similar molecules: {len(similar_df)}")

sample = similar_df.sample(min(30, len(similar_df)))
mols = [Chem.MolFromSmiles(s) for s in sample['SMILES']]
legends = [f"{tid}\n{sim:.2f}" for tid, sim in 
           zip(sample['Triazine_ID'], sample['max_similarity'])]

img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(350,300), legends=legends)


# In[57]:


bis_ids = [
    'TRZ00686', 'TRZ00727', 'TRZ01194', 'TRZ01492', 'TRZ02619', 'TRZ02802'
]

top_bis_df = bis_df[
    bis_df['Triazine_ID'].isin(bis_ids)
].reset_index(drop=True)

len(top_bis_df)


# In[58]:


mono_ids = [
    'TRZ05503', 'TRZ05722', 'TRZ05495', 'TRZ05736'
]

top_mono_df = similar_df[
    ~similar_df['Triazine_ID'].isin(mono_ids)
].reset_index(drop=True)

len(top_mono_df)


# In[59]:


top_df = pd.concat([top_bis_df, top_mono_df])
top_df.to_csv('results/Metal_sensing_mol.csv', index=False)     
top_df.head()


#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd

top_df = pd.read_csv("results/Metal_sensing_mol.csv")


# In[3]:


import os
import subprocess

RESULTS_DIR = os.path.join("results","quantum_screening_data")
INPUT_DIR = os.path.join(RESULTS_DIR, "xtb_inputs")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "xtb_outputs")
SCRATCH_DIR = os.path.join(RESULTS_DIR, "xtb_scratch")

for folder in [INPUT_DIR, OUTPUT_DIR, SCRATCH_DIR]:
    os.makedirs(folder, exist_ok=True)

XTB_PATH = r"D:\xtb-6.6.1\bin\xtb.exe"


# In[4]:


'''Generating XYZ coordinates for top molecules using molecular mechanics'''

import shutil
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_3d_xyz(smiles, filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    mol = Chem.AddHs(mol)  # Add explicit hydrogens required for 3D coordinates

    # Embed using the ETKDGv3 algorithm
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if status != 0:
        return False

    # Quick structural minimization via MMFF94 force field
    AllChem.MMFFOptimizeMolecule(mol)

    # Save coordinates as an XYZ file
    Chem.MolToXYZFile(mol, filename)
    return True

success_count = 0

for idx, row in top_df.iterrows():
    mol_id = row['Triazine_ID']
    smi = row['SMILES']

    filepath = os.path.join(INPUT_DIR, f"{mol_id}.xyz")

    if generate_3d_xyz(smi, filepath):
        success_count += 1
    else:
        print(f"Failed to generate 3D structure for {mol_id}")


# In[5]:


'''Running Semi-Empirical quantum calculation using xtb to find band gap, HOMO, LUMO and Tot. energy'''

xtb_data = []

# Loop through the generated XYZ inputs
for filename in sorted(os.listdir(INPUT_DIR)):
    if not filename.endswith('.xyz'):
        continue

    mol_id = filename.split('.')[0]
    input_path = os.path.abspath(os.path.join(INPUT_DIR, filename))
    log_filepath = os.path.abspath(os.path.join(OUTPUT_DIR, f"{mol_id}.log"))

    # Run command: optimize geometry using GFN2 settings
    cmd = f'"{XTB_PATH}" "{input_path}" --opt --gfn 2'

    # Execute the command inside SCRATCH_DIR to keep it clean
    with open(log_filepath, "w", encoding="utf-8") as log_file:
        subprocess.run(
            cmd, 
            shell=True, 
            cwd=SCRATCH_DIR, 
            stdout=log_file, 
            stderr=subprocess.PIPE
        )

    homo = None
    lumo = None
    total_energy = None

    if os.path.exists(log_filepath):
        with open(log_filepath, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue

                if "(HOMO)" in line:
                    homo = float(parts[-2]) 
                elif "(LUMO)" in line:
                    lumo = float(parts[-2])
                elif "TOTAL ENERGY" in line:
                    try:
                        energy_idx = parts.index("ENERGY") + 1
                        total_energy = float(parts[energy_idx])
                    except (ValueError, IndexError):
                        # Fallback if the table format changes slightly
                        total_energy = float(parts[-3])

        if homo is not None and lumo is not None:
            gap = lumo - homo
            xtb_data.append({
                "Triazine_ID": mol_id,
                "HOMO_eV": homo,
                "LUMO_eV": lumo,
                "Gap_eV": gap,
                "E_ligand_au": total_energy
            })

quantum_screening_df = pd.DataFrame(xtb_data)
quantum_df_path = os.path.join(RESULTS_DIR, "xtb_baseline_screening.csv")
quantum_screening_df.to_csv(quantum_df_path, index=False)

print(quantum_screening_df.head())


# In[6]:


df = pd.merge(
    top_df, 
    quantum_screening_df, 
    on="Triazine_ID"
)
print(df.head())


# In[7]:


# Add donor set column manually based on your earlier analysis
donor_sets = {
    'TRZ00686': 'N,N,N,N',  # tetradentate
    'TRZ00727': 'N,O,N',
    'TRZ01194': 'N,N,N',
    'TRZ01492': 'N,S,N',
    'TRZ02619': 'N,S,N',
    'TRZ02802': 'N,N,N',
    'TRZ00111': 'N,O,O',
    'TRZ00112': 'N,O,O',
    'TRZ00116': 'N,O,O,O', # tetradentate
    'TRZ00117': 'N,O,O',
    'TRZ00127': 'N,O,O',
    'TRZ00259': 'N,O,O',
    'TRZ00296': 'N,O,O',
    'TRZ00371': 'N,O,O',
    'TRZ00374': 'N,O,O',
    'TRZ00375': 'N,O,O,O',  # tetradentate
    'TRZ00376': 'N,O,N',
    'TRZ00432': 'N,O,O',
    'TRZ01251': 'N,O,O',
    'TRZ01544': 'N,O,S',
    'TRZ02396': 'N,N,N',
    'TRZ00399': 'N,S,S',
    'TRZ00628': 'N,N,N',
    'TRZ01080': 'N,N,O',
    'TRZ01527': 'N,N,N',
    'TRZ02000': 'N,N,O',
    'TRZ02390': 'N,S,O',
    'TRZ02392': 'N,N,S',
    'TRZ05498': 'N,N,O',
    'TRZ05734': 'N,N,O'
}

df['Donor_Set'] = df['Triazine_ID'].map(donor_sets)
print(df['Donor_Set'].value_counts())


# In[8]:


def calculate_cu2_score_v2(row):

    # Factor 1 — Donor Set (50 points max)
    donor_scores = {
        'N,N,N,N': 50,
        'N,O,O,O': 42,
        'N,N,N':   45,
        'N,N,O':   40,
        'N,N,S':   35,
        'N,O,N':   40,
        'N,O,O':   35,
        'N,O,S':   30,
        'N,S,O':   30,
        'N,S,N':   35,
        'N,S,S':   25
    }
    donor_score = donor_scores.get(row['Donor_Set'], 0)

    # Factor 2 — Band Gap (30 points max)
    gap = row['Gap_eV']
    if 1.0 <= gap <= 1.5:
        gap_score = 30        # ideal range
    elif gap < 1.0:
        gap_score = 15        # too small = unstable
    else:
        gap_score = 9         # too large = weak sensing

    # Factor 3 — Bis bonus (10 points max)
    bis_score = 10 if row['n_triazine_rings'] == 2 else 5

    total = donor_score + gap_score + bis_score

    return round(total, 2)

df['Cu2_Score'] = df.apply(calculate_cu2_score_v2, axis=1)

# Rank
ranked = df[['Triazine_ID','SMILES', 'Donor_Set', 'HOMO_eV','MolWt', 'LogP', 'TPSA', 'HBA',
                        'LUMO_eV', 'Gap_eV', 'n_triazine_rings','Cu2_Score']] \
            .sort_values('Cu2_Score', ascending=False) \
            .reset_index(drop=True)

ranked['Rank'] = ranked.index + 1

print(ranked[['Rank', 'Triazine_ID', 'Donor_Set',
                  'Gap_eV', 'n_triazine_rings',
                  'Cu2_Score']].to_string())


# In[9]:


from rdkit.Chem.Draw import MolsToGridImage
from IPython.display import display

mols = [Chem.MolFromSmiles(s) for s in ranked['SMILES']]

# Legend showing rank, score and donor set
legends = [
    f"#{row['Rank']} {row['Triazine_ID']}\n{row['Donor_Set']} | Score:{row['Cu2_Score']}\n{row['Gap_eV']:.4f}"
    for _, row in ranked.iterrows()
]

img = MolsToGridImage(mols, molsPerRow=5,
                      subImgSize=(350,300),
                      legends=legends)

img.save('results/figures/top30_molecules.png')

# In[10]:


import matplotlib.pyplot as plt
import numpy as np

top5 = ranked.head(5)

categories = ['Cu2_Score', 'HOMO_eV', 'LUMO_eV', 'Gap_eV', 
              'MolWt', 'LogP', 'TPSA', 'HBA']

fig, axes = plt.subplots(1, 5, figsize=(25,5),
                         subplot_kw=dict(polar=True))

for idx, (_, row) in enumerate(top5.iterrows()):
    values = [row[c] for c in categories]

    # Normalize each value 0-1
    values_norm = [(v - ranked[c].min()) / 
                   (ranked[c].max() - ranked[c].min() + 1e-9)
                   for v, c in zip(values, categories)]

    angles = np.linspace(0, 2*np.pi, len(categories), 
                         endpoint=False).tolist()
    values_norm += values_norm[:1]
    angles += angles[:1]

    axes[idx].plot(angles, values_norm, 'b-', linewidth=2)
    axes[idx].fill(angles, values_norm, alpha=0.25)
    axes[idx].set_xticks(angles[:-1])
    axes[idx].set_xticklabels(categories, size=8)
    axes[idx].set_title(f"#{row['Rank']} {row['Triazine_ID']}\n{row['Donor_Set']}", 
                        size=9)

plt.tight_layout()
plt.savefig('results/figures/radar_for_top5.png', dpi=150)
ranked.to_csv('results/electronic_prop.csv', index=False) 
plt.show()

