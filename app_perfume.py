import streamlit as st
import sys
import os
import io
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
import base64
import pandas as pd

# Add the src directory to the path for logic_cyp
sys.path.append('src/')

from logic_cyp import (
    load_data,
    ask_preferences,
    score_perfumes,
    get_molecules_for_scents
)

# --- Streamlit App ---
st.set_page_config(page_title="Choose Your Perfume", layout="centered")
st.title("🌸 Choose Your Perfume")

@st.cache_data
def cached_load_data():
    return load_data()

perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df = cached_load_data()

# --- Step 1: Ask Preferences ---
st.header("Tell us about the scents you love")
scent_options = ask_preferences()

selected_scents = st.multiselect(
    "What scent types do you enjoy?",
    options=scent_options
)

# --- Step 2: Show Molecules for Selected Scents ---
if selected_scents:
    st.subheader("🔬 Molecules related to your scent preferences")
    molecule_df = get_molecules_for_scents(selected_scents, scent_to_smiles_df)
    st.dataframe(molecule_df)



# Create a table with SMILES and structure image
rows = []

for smiles in molecule_df["nonStereoSMILES"]:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(100, 100))  # Smaller image
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode()
        img_html = f'<img src="data:image/png;base64,{b64}" width="100">'
        rows.append({"SMILES": smiles, "Structure": img_html})
    else:
        rows.append({"SMILES": smiles, "Structure": "Invalid SMILES"})

# Create and render the styled table
df = pd.DataFrame(rows)
st.write("🧬 Molecule Structures Table")
st.write(df.to_html(escape=False), unsafe_allow_html=True)


    
# --- Step 3: Show Recommended Perfumes ---
if selected_scents:
    st.subheader("✨ Perfumes that match your preferences")
    top_perfumes = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df)

    preferred_columns = ["PerfumeID", "PerfumeName", "id", "name", "score"]
    available_columns = [col for col in preferred_columns if col in top_perfumes.columns]

    if not {"score"}.intersection(available_columns):
        st.warning("No score or perfume name columns found in the data. Please check your dataset.")
    else:
        st.write("Here are your top matches:")
        st.dataframe(top_perfumes[available_columns].head(5))
