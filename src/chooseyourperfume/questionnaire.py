import streamlit as st
from dataset import load_smiles_odors, load_perfume_descriptions
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

def render_smiles(smiles_str):
    """Convert a SMILES string into an RDKit molecule image."""
    mol = Chem.MolFromSmiles(smiles_str)
    if mol:
        img = Draw.MolToImage(mol, size=(200, 200))
        return img
    return None

def scent_questionnaire():
    st.header("🔍 Explore Perfumes by Scent")

    df_smells = load_smiles_odors()
    df_perfumes = load_perfume_descriptions()

    all_scents = sorted({scent.strip() for scents in df_smells['Odor'].dropna() for scent in scents.split(',')})

    selected_scent = st.selectbox("Choose a scent note you're interested in:", all_scents)

    matched_molecules = df_smells[df_smells['Odor'].str.contains(selected_scent, case=False, na=False)]

    st.write(f"### Molecules with the scent '{selected_scent}':")

    for _, row in matched_molecules.iterrows():
        col1, col2 = st.columns([1, 3])
        with col1:
            img = render_smiles(row['SMILES'])
            if img:
                st.image(img)
            else:
                st.write("Invalid SMILES")
        with col2:
            st.markdown(f"**SMILES:** `{row['SMILES']}`")
            st.markdown(f"**Odor Tags:** {row['Odor']}")

    matching_perfumes = df_perfumes[df_perfumes['accords'].str.contains(selected_scent, case=False, na=False)]

    st.write(f"### Perfumes that likely feature '{selected_scent}':")
    st.dataframe(matching_perfumes[['name', 'brand', 'accords', 'description']])