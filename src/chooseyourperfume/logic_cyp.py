import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, DataStructs
import streamlit as st
import re

# Dataset imports
try:
    from dataset import scent_categories
    from dataset import (
        load_perfume_descriptions,
        load_fragrantica_data,
        load_extended_perfume_set,
        load_smiles_odors
    )
except ImportError:
    from .dataset import scent_categories
    from .dataset import (
        load_perfume_descriptions,
        load_fragrantica_data,
        load_extended_perfume_set,
        load_smiles_odors
    )


def load_data():
    perfume_desc_df = load_perfume_descriptions()
    perfume_df = load_extended_perfume_set()
    scent_to_smiles_df = load_smiles_odors()

    # Standardize columns
    perfume_desc_df.columns = perfume_desc_df.columns.str.strip().str.lower()
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()
    scent_to_smiles_df.columns = scent_to_smiles_df.columns.str.strip().str.lower()

    # Extract all scent keywords
    all_scent_notes = [note.strip().lower() for notes in scent_categories.values() for note in notes]

    # Enrich with scent indicator columns
    perfume_to_scent_df = enrich_with_scent_columns(perfume_desc_df.copy(), all_scent_notes, text_column='description')

    return perfume_to_scent_df, perfume_desc_df, perfume_df, scent_to_smiles_df


def enrich_with_scent_columns(df, scent_list, text_column='description'):
    if text_column not in df.columns:
        print(f"⚠️ Text column '{text_column}' not found in columns: {df.columns}")
        return df
    for scent in scent_list:
        scent_clean = scent.strip().lower()
        df[scent_clean] = df[text_column].str.contains(scent, case=False, na=False).astype(int)
    return df


def score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights=None):
    perfume_scores = perfume_to_scent_df.copy()
    perfume_scores.columns = perfume_scores.columns.str.strip().str.lower()
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()

    selected_scents = [scent.strip().lower() for scent in selected_scents]

    # Validate and apply weights
    if weights is None:
        weights = {scent: 1.0 for scent in selected_scents}
    else:
        weights = {k.strip().lower(): v for k, v in weights.items()}

    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    if not valid_scents:
        st.warning("🚫 None of the selected scents matched our dataset. Try different notes.")
        return pd.DataFrame(columns=["score"])

    # Percentage-based scoring
    perfume_scores['match_count'] = perfume_scores[valid_scents].sum(axis=1)
    perfume_scores['score'] = (perfume_scores['match_count'] / len(valid_scents)) * 100

    # Clean perfume names for robust merging
    def clean_name(name):
        name = str(name).lower()
        name = re.sub(r'\s+', ' ', name)
        return name.strip()

    perfume_scores["name_clean"] = perfume_scores["name"].astype(str).apply(clean_name)
    perfume_df["name_clean"] = perfume_df["name"].astype(str).apply(clean_name)

    result = perfume_scores.merge(perfume_df, on="name_clean", how="left")
    result = result.sort_values(by="score", ascending=False)
    return result


def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))


def ask_preferences():
    return scent_categories


def get_molecules_for_scents(selected_scents, scent_to_smiles_df):
    valid_scents = [scent for scent in selected_scents if scent in scent_to_smiles_df.columns]
    if not valid_scents:
        return pd.DataFrame(columns=scent_to_smiles_df.columns)
    mask = scent_to_smiles_df[valid_scents].sum(axis=1) > 0
    return scent_to_smiles_df[mask]


def avg_similarity(smiles_list, radius: int = 2, n_bits: int = 1024) -> float | None:
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), radius, n_bits)
           for smi in smiles_list if Chem.MolFromSmiles(smi)]
    if len(fps) < 2:
        return None
    sims = [DataStructs.TanimotoSimilarity(fps[i], fps[j])
            for i in range(len(fps)) for j in range(i + 1, len(fps))]
    return sum(sims) / len(sims)

