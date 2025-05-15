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

# CSS + Google Fonts Injection
st.markdown(
    """
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond&display=swap" rel="stylesheet">
    <style>
        /* Background */
        .stApp { 
            background-color: #F4EDDE !important; 
        }

        /* Police */
        html, body, [class*="css"]  {
            font-family: 'Cormorant Garamond', serif !important;
            color: #4C3A32 !important;
        }
        
        /* Titres */
        h1, h2, h3, h4, h5, h6 {
            font-family: 'Cormorant Garamond', serif !important;
            color: #4C3A32 !important;
        }
    </style>
    """,
    unsafe_allow_html=True
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
        "<h1> ⋆ CHOOSE YOUR PERFUME ⋆ </h1><p style='color:gray;'>Find your signature scent</p>",
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
selected_scents = list(set(
    note for cat in categories for note in st.session_state.get(f"sel_{cat}", [])
))
st.write("**Your picks:**", ", ".join(selected_scents) if selected_scents else "None yet.")

# --- Step 2: Weight Sliders ---
weights = {}
if selected_scents:
    st.subheader("2. How much do you love each scent?")
    
    slider_labels = {
        0.1: "It's okay",
        0.3: "I like",
        0.7: "I love",
        1.0: "I adore",
        1.5: "I only want this scent!"
    }
    slider_steps = list(slider_labels.keys())

    for note in set(selected_scents):  # Ensure each scent appears only once
        val = st.select_slider(
            f"{note}",
            options=slider_steps,
            value=0.7,
            format_func=lambda x: slider_labels[x],
            key=f"w_{note}"  # Now this key is unique per scent
        )
        weights[note] = val
else:
    weights = {}

# --- Step 3: Generate Recommendations ---
if st.button("🔍 Generate Recommendations"):
    if not selected_scents:
        st.warning("Please pick at least one note before generating.")
    else:
        # Results Panes: Perfume Recommendations First
        col_perf, col_mol = st.columns(2, gap="large")

        with col_perf:
            st.subheader("⚭ Perfume Matches")
            top = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights)
            cols = [c for c in ['PerfumeName', 'brand', 'score'] if c in top.columns]
            st.dataframe(top[cols].head(5), use_container_width=True)

        # Now show Molecules based on Weights after Recommendations
    with col_mol:
        st.subheader("⌬ Explore Molecules Based on Your Preferences")

        for scent in selected_scents:
            if scent not in scent_to_smiles_df.columns:
                continue

            with st.expander(f"Molecules carrying the scent note '{scent}'", expanded=False):
                subset = scent_to_smiles_df[scent_to_smiles_df[scent] == 1]
                mol_entries = []
                for smi in subset['nonStereoSMILES']:
                    m = Chem.MolFromSmiles(smi)
                    if m:
                        img = Draw.MolToImage(m, (120, 120))  # Adjust size as needed
                        buf = io.BytesIO(); img.save(buf, 'PNG')
                        img_base64 = base64.b64encode(buf.getvalue()).decode()
                        mol_entries.append({'img_base64': img_base64})

                        
                if not mol_entries:
                    st.write("No molecules found for this scent.")
                else:
                    html_content = """
                    <style>
                    /* Custom scrollbar styling */
                    .scroll-container {
                        display: flex;
                        flex-wrap: nowrap;  /* Force horizontal layout */
                        overflow-x: auto;
                        gap: 20px;
                        padding: 10px;
                        white-space: nowrap;
                        scrollbar-width: thin;
                        scrollbar-color: #c2b280 #f4eddd;
                    }

                    .scroll-container::-webkit-scrollbar {
                        height: 8px;
                    }

                    .scroll-container::-webkit-scrollbar-track {
                        background: #f4eddd;
                        border-radius: 4px;
                    }

                    .scroll-container::-webkit-scrollbar-thumb {
                        background-color: #c2b280;
                        border-radius: 4px;
                    }
                    </style>

                    <div class="scroll-container">
                    """
                    for mol in mol_entries:
                        html_content += f"""
                    <div style="text-align: center; min-width: 140px;">
                        <img src="data:image/png;base64,{mol['img_base64']}" 
                            style="border-radius:8px; box-shadow:0px 4px 6px rgba(0,0,0,0.1);"/>
                    </div>
                    """
                    html_content += "</div>"

                    st.markdown(html_content, unsafe_allow_html=True)