import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, DataStructs
import streamlit as st
from .dataset import scent_categories
from .dataset import (
    load_perfume_descriptions,
    load_fragrantica_data,
    load_extended_perfume_set,
    load_smiles_odors
)

def load_data():
    perfume_clean_df = load_perfume_descriptions()  # switched to dataset with 'description'
    scent_to_smiles_df = load_smiles_odors()

    # Standardize column names
    perfume_clean_df.columns = perfume_clean_df.columns.str.strip().str.lower()
    scent_to_smiles_df.columns = scent_to_smiles_df.columns.str.strip().str.lower()

    all_scent_notes = [note.strip().lower() for notes in scent_categories.values() for note in notes]
    perfume_to_scent_df = enrich_with_scent_columns(perfume_clean_df.copy(), all_scent_notes, text_column='description')

    return perfume_to_scent_df, perfume_clean_df, perfume_clean_df, scent_to_smiles_df

def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))

def ask_preferences():
    return scent_categories

def enrich_with_scent_columns(perfume_clean_df, scent_list, text_column='description'):
    if text_column not in perfume_clean_df.columns:
        print(f"\u26a0\ufe0f Text column '{text_column}' not found in perfume_clean_df columns: {perfume_clean_df.columns}")
        return perfume_clean_df
    for scent in scent_list:
        scent_clean = scent.strip().lower()
        perfume_clean_df[scent_clean] = perfume_clean_df[text_column].str.contains(scent, case=False, na=False).astype(int)
    return perfume_clean_df

def score_perfumes(selected_scents, perfume_to_scent_df, perfume_clean_df, weights=None):
    import re

    perfume_scores = perfume_to_scent_df.copy()

    # Normalize column names
    perfume_scores.columns = perfume_scores.columns.str.strip().str.lower()
    perfume_clean_df.columns = perfume_clean_df.columns.str.strip().str.lower()

    selected_scents = [scent.strip().lower() for scent in selected_scents]

    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    if not valid_scents:
        st.warning("\ud83d\udeab None of the selected scents matched our dataset. Try different notes.")
        return pd.DataFrame(columns=["score"])

    # Percentage scoring: count how many valid scents match per perfume
    perfume_scores['match_count'] = perfume_scores[valid_scents].sum(axis=1)
    perfume_scores['score'] = (perfume_scores['match_count'] / len(valid_scents)) * 100

    def clean_name(name):
        name = str(name).lower()
        name = re.sub(r'\b(for men|for women|for women and men|pour femme|pour homme|by|pour)\b', '', name)
        name = re.sub(r'\s+', ' ', name)
        return name.strip()

    perfume_scores["name"] = perfume_scores["name"].astype(str).str.lower().str.strip()
    perfume_clean_df["name"] = perfume_clean_df["name"].astype(str).str.lower().str.strip()
    perfume_scores["name_clean"] = perfume_scores["name"].apply(clean_name)
    perfume_clean_df["name_clean"] = perfume_clean_df["name"].apply(clean_name)

    result = perfume_scores.merge(perfume_clean_df, on="name_clean", how="left")

    result = result.sort_values(by="score", ascending=False)


    return result

def get_molecules_for_scents(selected_scents, scent_to_smiles_df):
    valid_scents = [scent for scent in selected_scents if scent in scent_to_smiles_df.columns]
    if not valid_scents:
        return pd.DataFrame(columns=scent_to_smiles_df.columns)
    mask = scent_to_smiles_df[valid_scents].sum(axis=1) > 0
    return scent_to_smiles_df[mask]

def split_name_brand(name_string):
    if not isinstance(name_string, str):
        return "Unknown", "Unknown"
    parts = name_string.strip().split()
    if len(parts) >= 2:
        brand = parts[-1]
        perfume_name = " ".join(parts[:-1])
    else:
        brand = "Unknown"
        perfume_name = name_string
    return perfume_name, brand

def display_results(result):
    st.markdown("## Perfume Matches")
    if result.empty:
        st.warning(" No matching perfumes found.")
    else:
        weights = {k.strip().lower(): v for k, v in weights.items()}

    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]

    if not valid_scents:
        st.warning("🚫 None of the selected scents matched our dataset. Try different notes.")
        return pd.DataFrame(columns=["score"])

    for scent in valid_scents:
        weight = weights.get(scent, 1.0)
        perfume_scores["score"] += perfume_scores[scent] * weight

    possible_keys = ['name', 'perfumename', 'perfume_name']
    merge_key = next((key for key in possible_keys if key in perfume_scores.columns and key in perfume_df.columns), None)
    print("Selected merge key:", merge_key)
    if not merge_key:
        st.error("❌ Merge key missing. Check that perfume names are consistently labeled in datasets.")
        raise KeyError("No valid merge key found.")

    perfume_scores[merge_key] = perfume_scores[merge_key].astype(str).str.strip().str.lower()
    perfume_df[merge_key] = perfume_df[merge_key].astype(str).str.strip().str.lower()

    # Debug: Verify sample keys after normalization
    print("Sample Keys in perfume_scores (after cleanup):", perfume_scores[merge_key].dropna().unique()[:5])
    print("Sample Keys in perfume_df (after cleanup):", perfume_df[merge_key].dropna().unique()[:5])

    result = perfume_scores.merge(perfume_df, on=merge_key, how="left")
    print("Merge result sample:\n", result.head())
    result = result.sort_values(by="score", ascending=False)

    return result

def get_molecules_for_scents(selected_scents, scent_to_smiles_df):
    valid_scents = [scent for scent in selected_scents if scent in scent_to_smiles_df.columns]
    if not valid_scents:
        return pd.DataFrame(columns=scent_to_smiles_df.columns)
    mask = scent_to_smiles_df[valid_scents].sum(axis=1) > 0
    return scent_to_smiles_df[mask]