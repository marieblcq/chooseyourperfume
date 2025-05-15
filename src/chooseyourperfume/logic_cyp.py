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

    # Get the complete list of scent notes used in the app
    from .dataset import scent_categories
    all_scent_notes = [note.strip().lower() for notes in scent_categories.values() for note in notes]

    # Enrich with scent presence columns if they don't exist already
    perfume_to_scent_df = enrich_with_scent_columns(perfume_to_scent_df, all_scent_notes, text_column='description')

    return perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df


def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))



def ask_preferences():
    return scent_categories 


def score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights=None):
    import streamlit as st  # Import here to allow Streamlit notifications

    perfume_scores = perfume_to_scent_df.copy()
    perfume_scores["score"] = 0.0  # Ensure floating point for calculations

    # Normalize column names for consistent matching
    perfume_scores.columns = perfume_scores.columns.str.strip().str.lower()
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()
    selected_scents = [scent.strip().lower() for scent in selected_scents]
    if weights is None:
        weights = {scent: 1.0 for scent in selected_scents}
    else:
        weights = {k.strip().lower(): v for k, v in weights.items()}

    # Debugging: Check available columns and selected scents
    print("DEBUG ➤ Columns in perfume_scores:", list(perfume_scores.columns))
    print("DEBUG ➤ Columns in perfume_df:", list(perfume_df.columns))
    print("DEBUG ➤ Selected scents (normalized):", selected_scents)
    print("DEBUG ➤ Weights being applied:", weights)

    # Check valid scents that exist as columns
    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    print("DEBUG ➤ Valid scents found in DataFrame:", valid_scents)

    if not valid_scents:
        st.warning("🚫 None of the selected scents matched our dataset. Try different notes.")
        return pd.DataFrame(columns=["score"])

    # Apply weighted scoring
    for scent in valid_scents:
        weight = weights.get(scent, 1.0)
        perfume_scores["score"] += perfume_scores[scent] * weight

    # Debugging: Check scores before merging
    print("DEBUG ➤ Top scores before merging:\n", perfume_scores[['score']].sort_values(by='score', ascending=False).head())

    # Explicitly define merge key (confirm it exists!)
    possible_keys = ['perfumename', 'perfume_name', 'name']
    merge_key = next((key for key in possible_keys if key in perfume_scores.columns and key in perfume_df.columns), None)

    if not merge_key:
        st.error("❌ Merge key missing. Check that perfume names are consistently labeled in datasets.")
        raise KeyError(f"No valid merge key found. Tried {possible_keys}. "
                       f"perfume_scores columns: {perfume_scores.columns}, perfume_df columns: {perfume_df.columns}")

    print(f"DEBUG ➤ Using merge key: {merge_key}")

    # Merge and sort results
    result = perfume_scores.merge(perfume_df, on=merge_key, how="left")
    result = result.sort_values(by="score", ascending=False)

    if result.empty:
        st.warning("🚫 No matching perfumes found. Try adjusting your scent preferences or weights.")
        print("DEBUG ➤ Final DataFrame is empty after merging.")
    else:
        print("DEBUG ➤ Top perfumes after merging:\n", result[[merge_key, 'score']].head())

    return result



# Get molecules (SMILES) related to selected scent categories

def get_molecules_for_scents(selected_scents, scent_to_smiles_df):
    # Ensure selected_scents are valid column names
    valid_scents = [scent for scent in selected_scents if scent in scent_to_smiles_df.columns]

    if not valid_scents:
        return pd.DataFrame(columns=scent_to_smiles_df.columns)  # empty result if nothing matches

    # Keep rows where **any** selected scent is present (value > 0)
    mask = scent_to_smiles_df[valid_scents].sum(axis=1) > 0
    filtered = scent_to_smiles_df[mask]

    return filtered

def enrich_with_scent_columns(perfume_df, scent_list, text_column='description'):
    # Normalize column names
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()

    # Ensure description column exists
    if text_column not in perfume_df.columns:
        print(f"Column '{text_column}' not found in perfume DataFrame.")
        return perfume_df

    for scent in scent_list:
        scent_clean = scent.strip().lower()
        perfume_df[scent_clean] = perfume_df[text_column].str.contains(scent, case=False, na=False).astype(int)

    return perfume_df  # Final correct return

