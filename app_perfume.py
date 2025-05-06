import streamlit as st
import sys
import os

# Add the src directory to the path
sys.path.append('src/')

# Correct import from logic_cyp directly in src
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

# --- Step 3: Show Recommended Perfumes ---
if selected_scents:
    st.subheader("✨ Perfumes that match your preferences")
    top_perfumes = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df)
    
    # Show top 5 perfumes
    st.write("Here are your top matches:")
    # Figure out which columns to display safely
preferred_columns = ["PerfumeID", "PerfumeName", "id", "name", "score"]
available_columns = [col for col in preferred_columns if col in top_perfumes.columns]

# Show a warning if expected columns are missing
if not {"score"}.intersection(available_columns):
    st.warning("No score or perfume name columns found in the result. Please check your dataset.")
else:
    st.dataframe(top_perfumes[available_columns].head(5))


