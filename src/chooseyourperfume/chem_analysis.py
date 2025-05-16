# --- Import Libraries ---
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFingerprintGenerator
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

# --- Import Dataset Loader ---
from .dataset import load_smiles_odors

# --- Load Data ---
data = load_smiles_odors()

# --- Canonicalize SMILES ---
def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)

data['canonical_smiles'] = data['nonStereoSMILES'].apply(canonicalize_smiles)
data = data.dropna(subset=['canonical_smiles'])

# --- Create Mol Objects Column ---
data['mol'] = data['canonical_smiles'].apply(Chem.MolFromSmiles)

# --- Define Scent Categories ---
scent_categories = {
    "Floral": ['floral', 'jasmin', 'rose', 'violet', 'lily', 'hyacinth', 'lavender', 'muguet', 'chamomile', 'orangeflower', 'geranium'],
    "Fruity": ['fruity', 'apple', 'apricot', 'banana', 'berry', 'black currant', 'citrus', 'grape', 'grapefruit', 'melon', 'orange', 'peach', 'pear', 'pineapple', 'plum', 'raspberry', 'strawberry', 'tropical'],
    "Vegetal / Herbal": ['green', 'grassy', 'herbal', 'leafy', 'celery', 'cucumber', 'hay', 'hawthorn', 'weedy', 'vegetable'],
    "Sweet / Gourmand": ['sweet', 'creamy', 'coconut', 'vanilla', 'chocolate', 'caramellic', 'buttery', 'honey', 'almond', 'milky'],
    "Nutty / Seed": ['nutty', 'hazelnut', 'almond'],
    "Lactonic / Milky": ['lactonic', 'dairy', 'milky', 'creamy'],
    "Spicy / Aromatic": ['spicy', 'clove', 'cinnamon', 'anisic', 'mint', 'coumarinic', 'camphoreous'],
    "Fresh / Volatile": ['fresh', 'cooling', 'ethereal', 'ozone', 'clean', 'citrus'],
    "Woody / Resinous": ['woody', 'sandalwood', 'cedar', 'pine', 'cortex', 'amber', 'vetiver'],
    "Animalic / Meaty": ['animal', 'meaty', 'beefy', 'musk', 'leathery', 'sweaty', 'savory'],
    "Smoky / Roasted / Burnt": ['smoky', 'burnt', 'roasted', 'cooked', 'popcorn', 'coffee'],
    "Fermented / Cheesy": ['cheesy', 'fermented', 'dairy'],
    "Sulfurous / Allium-like": ['sulfurous', 'garlic', 'onion', 'alliaceous', 'gassy', 'musty', 'cabbage', 'radish'],
    "Alcohol / Solvent": ['alcoholic', 'brandy', 'winey', 'cognac', 'solvent', 'rummy'],
    "Earthy / Mineral / Metallic": ['earthy', 'mushroom', 'metallic', 'oily', 'powdery'],
    "Chemical / Medicinal": ['medicinal', 'ketonic', 'phenolic', 'soapy', 'sharp', 'bitter'],
    "Dry / Bitter / Neutral": ['dry', 'bitter', 'natural', 'odorless']
}

# --- Reverse Mapping Scent to Category ---
scent_to_category = {}
for category, scent_list in scent_categories.items():
    for scent in scent_list:
        scent_to_category[scent] = category

# --- Assign Primary Scent Category ---
def assign_primary_category(descriptors):
    if pd.isna(descriptors):
        return None  # Mark for removal if no descriptors.

    descriptor_list = descriptors.split(';')
    category_counts = {}

    for desc in descriptor_list:
        desc = desc.strip()
        if desc in scent_to_category:
            category = scent_to_category[desc]
            category_counts[category] = category_counts.get(category, 0) + 1

    if not category_counts:
        return None  # Mark for removal if no known categories found.

    max_count = max(category_counts.values())
    tied_categories = [cat for cat, count in category_counts.items() if count == max_count]

    if len(tied_categories) > 1:
        return None  # Tie detected — remove this molecule.
    else:
        return tied_categories[0]  # Single dominant category, assign it.

data['primary_category'] = data['descriptors'].apply(assign_primary_category)
data = data.dropna(subset=['primary_category'])  # Completely remove molecules with no assigned category.

data['primary_category'] = data['descriptors'].apply(assign_primary_category)

# --- Generate Morgan Fingerprints ---
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def generate_fingerprint(mol):
    if mol is None:
        return None
    return morgan_gen.GetFingerprint(mol)

data['fingerprint'] = data['mol'].apply(generate_fingerprint)
data = data.dropna(subset=['fingerprint'])

# --- Convert Fingerprints to Numpy Arrays ---
def fingerprint_to_array(fp):
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fingerprints_array = np.array([fingerprint_to_array(fp) for fp in data['fingerprint']])

# --- Dimensionality Reduction ---
print("Running PCA...")
pca = PCA(n_components=50)
pca_result = pca.fit_transform(fingerprints_array)

print("Running t-SNE...")
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
tsne_result = tsne.fit_transform(pca_result)

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



