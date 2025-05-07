import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
import io

st.title("✅ RDKit Molecule Rendering Test")

# Force script to confirm execution
st.write("🔍 App loaded successfully!")

# Hardcoded SMILES strings
smiles_list = ["CCO", "O=C(O)Cc1ccccc1", "c1ccccc1"]

# Display each molecule
for smiles in smiles_list:
    st.write(f"Rendering: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(200, 200))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        st.image(buf.getvalue(), caption=smiles, use_column_width=True)
    else:
        st.warning(f"Invalid SMILES: {smiles}")
