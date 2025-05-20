import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import streamlit as st
from .dataset import get_scent_categories

scent_categories = get_scent_categories()

from .dataset import (
    load_perfume_descriptions,
    load_fragrantica_data,
    load_extended_perfume_set,
    load_smiles_odors
)

def load_data():
    perfume_clean_df = load_fragrantica_data()
    perfume_df = load_extended_perfume_set()
    scent_to_smiles_df = load_smiles_odors()

    # Standardize column names
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()
    perfume_clean_df.columns = perfume_clean_df.columns.str.strip().str.lower()
    scent_to_smiles_df.columns = scent_to_smiles_df.columns.str.strip().str.lower()

    all_scent_notes = [note.strip().lower() for notes in scent_categories.values() for note in notes]
    perfume_to_scent_df = enrich_with_scent_columns(perfume_df.copy(), all_scent_notes, text_column='description')

    return perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df

def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=(200, 200))

def ask_preferences():
    return scent_categories


def enrich_with_scent_columns(perfume_df, scent_list, text_column='description'):
    if text_column not in perfume_df.columns:
        print(f"⚠️ Text column '{text_column}' not found in perfume_df columns: {perfume_df.columns}")
        return perfume_df
    for scent in scent_list:
        scent_clean = scent.strip().lower()
        perfume_df[scent_clean] = perfume_df[text_column].str.contains(scent, case=False, na=False).astype(int)
    return perfume_df

def score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights=None, gender_preference="Any", perfume_clean_df=None):
    import re
    import pandas as pd
    import streamlit as st

    perfume_scores = perfume_to_scent_df.copy()
    perfume_scores["score"] = 0.0

    # Normalize column names
    perfume_scores.columns = perfume_scores.columns.str.strip().str.lower()
    perfume_df.columns = perfume_df.columns.str.strip().str.lower()

    selected_scents = [scent.strip().lower() for scent in selected_scents]
    st.write("📊 Gender preference selected:", gender_preference)

    # Apply weights
    if weights is None:
        weights = {scent: 1.0 for scent in selected_scents}
    else:
        weights = {k.strip().lower(): v for k, v in weights.items()}

    valid_scents = [scent for scent in selected_scents if scent in perfume_scores.columns]
    if not valid_scents:
        st.warning("🚫 None of the selected scents matched our dataset. Try different notes.")
        return pd.DataFrame(columns=["score"])

    for scent in valid_scents:
        weight = weights.get(scent, 1.0)
        perfume_scores["score"] += perfume_scores[scent] * weight

    # Clean names
    def clean_name(name):
        name = str(name).lower()
        name = re.sub(r'\b(for men|for women|for women and men|pour femme|pour homme|by|pour)\b', '', name)
        name = re.sub(r'\s+', ' ', name)
        return name.strip()

    perfume_scores["name"] = perfume_scores["name"].astype(str).str.lower().str.strip()
    perfume_df["name"] = perfume_df["name"].astype(str).str.lower().str.strip()
    perfume_scores["name_clean"] = perfume_scores["name"].apply(clean_name)
    perfume_df["name_clean"] = perfume_df["name"].apply(clean_name)

    st.write("Sample perfume_scores['name_clean']:", perfume_scores['name_clean'].unique()[:5])
    st.write("Sample perfume_df['name_clean']:", perfume_df['name_clean'].unique()[:5])

    result = perfume_scores.merge(perfume_df, on="name_clean", how="left")

    # --- 🔎 Fix: Identify correct gender column after merge
    gender_col = None
    for col in result.columns:
        if col.lower() == "gender":
            gender_col = col
            break
        elif "gender_" in col:
            gender_col = col
            break

    if gender_col:
        result[gender_col] = result[gender_col].astype(str).str.strip().str.lower()
        gender_filter = gender_preference.strip().lower()

        if gender_filter == "women":
            result = result[result[gender_col].isin(["for women", "for women and men"])]
        elif gender_filter == "men":
            result = result[result[gender_col].isin(["for men", "for women and men"])]
        elif gender_filter == "unisex":
            result = result[result[gender_col] == "for women and men"]

        st.write(f"✅ Number of results after filtering for '{gender_preference}':", len(result))
    else:
        st.warning("⚠️ 'gender' column not found in the merged data. Gender filtering skipped.")

    result = result.sort_values(by="score", ascending=False)
    st.write("🔎 Final result preview:", result.head())

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
    st.markdown("## ♾️ Perfume Matches")
    if result.empty:
        st.warning("😔 No matching perfumes found.")
    else:
        top_n = min(5, len(result))
        for i, (_, row) in enumerate(result.head(top_n).iterrows()):
            perfume_name, brand = split_name_brand(row.get("name_x", "Unknown"))
            score = row.get("score", 0.0)
            description = row.get("description_x", "No description available")

            if not perfume_name:
                perfume_name = "Unknown"
            if not brand:
                brand = "Unknown"

            st.subheader(f"{perfume_name} by {brand}")
            st.markdown(f"**Score:** {score:.2f}")
            st.markdown(f"**Description:** {description}")
            st.markdown("---")
