import streamlit as st
import sys
import os
import io
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

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

    if molecule_df.empty:
        st.warning("No molecules found for the selected scents.")
    else:
        st.dataframe(molecule_df)

        st.subheader("🧪 Molecule Structures")
        for smiles in molecule_df["nonStereoSMILES"]:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(200, 200))
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                st.image(buf.getvalue(), caption=smiles, use_column_width=True)
            else:
                st.warning(f"Invalid SMILES: {smiles}")

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
