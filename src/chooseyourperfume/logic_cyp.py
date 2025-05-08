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


def score_perfumes(selected_scents, perfume_to_scent_df, perfume_df):
    perfume_scores = perfume_to_scent_df.copy()
    perfume_scores["score"] = 0

    # Use only scents that exist in the scoring DataFrame
    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    if not valid_scents:
        return pd.DataFrame(columns=["score"])

    # Add up the score for each valid scent
    for scent in valid_scents:
        perfume_scores["score"] += perfume_scores[scent]

    # Detect the shared key column for merge
    merge_key = None
    for col in perfume_scores.columns:
        if col in perfume_df.columns:
            merge_key = col
            break

    if not merge_key:
        raise KeyError("No common column found between perfume_scores and perfume_df to merge on.")

    # Merge to add perfume metadata
    result = perfume_scores.merge(perfume_df, on=merge_key, how="left")

    # Sort by descending score
    result = result.sort_values(by="score", ascending=False)

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
