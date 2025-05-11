import streamlit as st
import io, random, base64
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

from src.chooseyourperfume.logic_cyp import (
    load_data, ask_preferences, score_perfumes, get_molecules_for_scents
)

# --- Config & CSS ---
st.set_page_config(page_title="Choose Your Perfume", layout="wide")
st.markdown(
    """
    <style>
      .stApp { background-color: #F4EDDE; }
      .css-1d391kg { background-color: #F4EDDE; }
    </style>
    """, unsafe_allow_html=True
)

# --- Data Loading ---
@st.cache_data
def cached_load_data():
    return load_data()
perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df = cached_load_data()

# --- Header ---
col1, col2 = st.columns([1, 4])
with col1:
    st.image("assets/logo.png", width=180)
with col2:
    st.markdown(
        "<h1>✨ Choose Your Perfume</h1><p style='color:gray;'>Find your signature scent</p>",
        unsafe_allow_html=True
    )

# --- Step 1: Preferences Selection ---
st.header("1. Tell us about the scents you love")

scent_dict = ask_preferences()
categories = list(scent_dict.keys())

# Surprise Me button
if st.button("🎲 Surprise Me!"):
    all_notes = [note for subs in scent_dict.values() for note in subs]
    picks = random.sample(all_notes, min(3, len(all_notes)))
    for cat in categories:
        st.session_state[f"sel_{cat}"] = [n for n in picks if n in scent_dict[cat]]
    st.success("✨ Surprise picks loaded! Scroll down to see and adjust.")

# Initialize session state
for cat in categories:
    st.session_state.setdefault(f"sel_{cat}", [])

# Two-column expanders for all categories
col_a, col_b = st.columns(2, gap="medium")
for idx, cat in enumerate(categories):
    target = col_a if idx % 2 == 0 else col_b
    with target.expander(cat, expanded=False):
        st.multiselect(
            cat,
            scent_dict[cat],
            key=f"sel_{cat}"
        )

# Show current picks
selected_scents = [note for cat in categories for note in st.session_state.get(f"sel_{cat}", [])]
st.write("**Your picks:**", ", ".join(selected_scents) if selected_scents else "None yet.")

# --- Step 2: Weight Sliders ---
weights = {}
if selected_scents:
    st.subheader("2. Rate how much you love each note")
    for note in selected_scents:
        weights[note] = st.slider(
            label=note,
            min_value=0.0,
            max_value=1.0,
            value=0.5,
            key=f"w_{note}"
        )
else:
    weights = {}

# --- Step 3: Generate Recommendations ---
if st.button("🔍 Generate Recommendations"):
    if not selected_scents:
        st.warning("Please pick at least one note before generating.")
    else:
        # Results panes
        col_mol, col_perf = st.columns(2, gap="large")
        with col_mol:
            st.subheader("🔬 Molecules")
            # Build a flat table: each selected scent paired with its molecules
            mol_entries = []
            for scent in selected_scents:
                if scent not in scent_to_smiles_df.columns:
                    continue
                # select rows where the scent flag is 1
                subset = scent_to_smiles_df[scent_to_smiles_df[scent] == 1]
                for smi in subset['nonStereoSMILES']:
                    img_html = 'Invalid SMILES'
                    m = Chem.MolFromSmiles(smi)
                    if m:
                        img = Draw.MolToImage(m, (80,80))
                        buf = io.BytesIO(); img.save(buf, 'PNG')
                        img_html = f'<img src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}" />'
                    mol_entries.append({'Odor Note': scent, 'SMILES': smi, 'Structure': img_html})
            mol_df = pd.DataFrame(mol_entries)
            if mol_df.empty:
                st.write("No molecules found for selected scents.")
            else:
                st.write(mol_df.to_html(escape=False), unsafe_allow_html=True)
        with col_perf:
            st.subheader("✨ Perfume Matches")
            top = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights)
            cols = [c for c in ['PerfumeName','brand','score'] if c in top.columns]
            st.dataframe(top[cols].head(5), use_container_width=True)


