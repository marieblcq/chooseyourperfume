import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFingerprintGenerator
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from .dataset import load_smiles_odors, scent_categories

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None

# --- Assign Primary Scent Category ---
scent_to_category = {scent: category for category, scent_list in scent_categories.items() for scent in scent_list}

def assign_primary_category(descriptors):
    if pd.isna(descriptors):
        return None

    descriptor_list = descriptors.split(';')
    category_counts = {}

    for desc in descriptor_list:
        desc = desc.strip()
        if desc in scent_to_category:
            category = scent_to_category[desc]
            category_counts[category] = category_counts.get(category, 0) + 1

    if not category_counts:
        return None

    max_count = max(category_counts.values())
    tied_categories = [cat for cat, count in category_counts.items() if count == max_count]

    return None if len(tied_categories) > 1 else tied_categories[0]

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def generate_fingerprint(mol):
    return morgan_gen.GetFingerprint(mol) if mol else None

def fingerprint_to_array(fp):
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def run_pipeline():
    global data, tsne_result  # So tests and script can access them

    data = load_smiles_odors()
    data['canonical_smiles'] = data['nonStereoSMILES'].apply(canonicalize_smiles)
    data = data.dropna(subset=['canonical_smiles'])

    data['mol'] = data['canonical_smiles'].apply(Chem.MolFromSmiles)
    data['primary_category'] = data['descriptors'].apply(assign_primary_category)
    data = data.dropna(subset=['primary_category'])

    data['fingerprint'] = data['mol'].apply(generate_fingerprint)
    data = data.dropna(subset=['fingerprint'])

    fingerprints_array = np.array([fingerprint_to_array(fp) for fp in data['fingerprint']])

    print("Running PCA...")
    pca = PCA(n_components=50)
    pca_result = pca.fit_transform(fingerprints_array)

    print("Running t-SNE...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    tsne_result = tsne.fit_transform(pca_result)

    return data, tsne_result

if __name__ == "__main__":
    data, tsne_result = run_pipeline()

    # --- Plot Results ---
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        x=tsne_result[:, 0],
        y=tsne_result[:, 1],
        hue=data['primary_category'],
        palette='tab20',
        alpha=0.7,
        legend='full'
    )
    plt.title("Molecule Clustering by Primary Scent Category (t-SNE Projection)")
    plt.xlabel("t-SNE Dimension 1")
    plt.ylabel("t-SNE Dimension 2")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()