import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import os

from .dataset import scent_categories
from .dataset import(
    load_csv,
    load_perfume_descriptions,
    load_fragrantica_data,
    load_extended_perfume_set,
    load_smiles_odors
)
def load_data():
    perfume_to_scent_df = load_perfume_descriptions()       # Dataset 1
    perfume_clean_df = load_fragrantica_data()              # Dataset 2
    perfume_df = load_extended_perfume_set()                # Dataset 3
    scent_to_smiles_df = load_smiles_odors()                # Dataset 4

    return perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df

def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))



def ask_preferences():
    return scent_categories 


def score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights=None):
    perfume_scores = perfume_to_scent_df.copy()
    perfume_scores["score"] = 0.0  # Ensure floating point for calculations

    # Ensure weights are properly initialized
    if weights is None:
        weights = {scent: 1.0 for scent in selected_scents}

    # Normalize all column names and selected scents for consistent matching
    perfume_scores.columns = perfume_scores.columns.str.strip().str.lower()
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()
    selected_scents = [scent.strip().lower() for scent in selected_scents]
    weights = {k.strip().lower(): v for k, v in weights.items()}

    # Debugging: Check what columns are available for scoring
    print("Available scent columns in perfume_to_scent_df:", list(perfume_scores.columns))
    print("Selected scents (normalized):", selected_scents)
    print("Weights being applied:", weights)

    # Use only scents that exist in the DataFrame columns
    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    print("Valid scents found in DataFrame:", valid_scents)

    if not valid_scents:
        print("No valid scent columns found for selected scents. Returning empty DataFrame.")
        return pd.DataFrame(columns=["score"])

    # Apply weighted scoring
    for scent in valid_scents:
        weight = weights.get(scent, 1.0)
        perfume_scores["score"] += perfume_scores[scent] * weight

    # Debugging: Show scores before merging
    print("Top scores after scoring:\n", perfume_scores[['score']].sort_values(by='score', ascending=False).head())

    # Explicitly define merge key (adjust if needed)
    merge_key = 'perfumename'  # Ensure this matches exactly in both DataFrames

    if merge_key not in perfume_scores.columns or merge_key not in perfume_df.columns:
        raise KeyError(f"'{merge_key}' column must exist in both DataFrames for merging. "
                       f"perfume_scores columns: {perfume_scores.columns}, perfume_df columns: {perfume_df.columns}")

    # Merge to add metadata
    result = perfume_scores.merge(perfume_df, on=merge_key, how="left")

    # Sort and handle empty result
    result = result.sort_values(by="score", ascending=False)

    if result.empty:
        print("No matching perfumes found after scoring and merging. Final DataFrame is empty.")
    else:
        print("Top perfumes after merging:\n", result[[merge_key, 'score']].head())

    return result




# Get molecules (SMILES) related to selected scent categories
import pandas as pd

def get_molecules_for_scents(selected_scents, scent_to_smiles_df):
    # Ensure selected_scents are valid column names
    valid_scents = [scent for scent in selected_scents if scent in scent_to_smiles_df.columns]

    if not valid_scents:
        return pd.DataFrame(columns=scent_to_smiles_df.columns)  # empty result if nothing matches

    # Keep rows where **any** selected scent is present (value > 0)
    mask = scent_to_smiles_df[valid_scents].sum(axis=1) > 0
    filtered = scent_to_smiles_df[mask]

    return filtered
