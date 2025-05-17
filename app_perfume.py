import streamlit as st
import io, random, base64
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from src.chooseyourperfume.logic_cyp import (
    load_data, ask_preferences, score_perfumes, get_molecules_for_scents
)

# --- Config & CSS ---
st.set_page_config(page_title="Choose Your Perfume", layout="wide")

st.markdown(
    """
    <link href="https://fonts.googleapis.com/css2?family=Cormorant+Garamond&display=swap" rel="stylesheet">
    <style>
        .stApp { background-color: #F4EDDE !important; }
        html, body, [class*="css"] {
            font-family: 'Cormorant Garamond', serif !important;
            color: #4C3A32 !important;
        }
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

if st.button("🎲 Surprise Me!"):
    all_notes = [note for subs in scent_dict.values() for note in subs]
    picks = random.sample(all_notes, min(3, len(all_notes)))
    for cat in categories:
        st.session_state[f"sel_{cat}"] = [n for n in picks if n in scent_dict[cat]]
    st.success("✨ Surprise picks loaded! Scroll down to see and adjust.")

for cat in categories:
    st.session_state.setdefault(f"sel_{cat}", [])

col_a, col_b = st.columns(2, gap="medium")
for idx, cat in enumerate(categories):
    target = col_a if idx % 2 == 0 else col_b
    with target.expander(cat, expanded=False):
        st.multiselect(
            " ",
            scent_dict[cat],
            key=f"sel_{cat}"
        )

selected_scents = list(set(
    note for cat in categories for note in st.session_state.get(f"sel_{cat}", [])
))
st.write("**Your picks:**", ", ".join(selected_scents) if selected_scents else "None yet.")

# --- Step 1.5: Gender Preference ---
st.subheader("Optional: Do you have a gender preference for your perfume?")
gender_preference = st.radio(
    "Choose Gender:",
    options=["Any", "Women", "Men", "Unisex"],
    index=0,
    horizontal=True
)

# --- Step 2: Weight Sliders ---
weights = {}
if selected_scents:
    st.subheader("2. How much do you love each scent?")
    slider_labels = {0.1: "It's okay", 0.3: "I like", 0.7: "I love", 1.0: "I adore", 1.5: "Obsessed!"}
    slider_steps = list(slider_labels.keys())

    for note in selected_scents:
        val = st.select_slider(
            f"{note}",
            options=slider_steps,
            value=0.7,
            format_func=lambda x: slider_labels[x],
            key=f"w_{note}"
        )
        weights[note] = val

# --- Step 3: Generate Recommendations ---
if st.button("🔍 Generate Recommendations"):
    if not selected_scents:
        st.warning("🚨 Please pick at least one note before generating recommendations.")
    else:
        with st.spinner("🔬 Finding your perfect perfumes..."):
            top = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights, gender_preference, perfume_clean_df)

        col_perf, col_mol = st.columns(2, gap="large")

        with col_perf:
            st.subheader("⚭ Perfume Matches")
            if top.empty:
                st.warning("🚫 No matching perfumes found. Try selecting different notes or adjusting weights.")
            else:
                for idx, row in top.head(5).iterrows():
                    perfume_name = row.get('name', 'Unknown')
                    brand = row.get('brand', 'Unknown')
                    score = row.get('score', 0)
                    description = row.get('description_x') or row.get('description_y') or 'No description available.'
                    ingredients = row.get('notes')

                    if pd.isna(description) or str(description).strip() == '':
                        description = 'No description available.'
                    if pd.isna(ingredients) or str(ingredients).strip() == '':
                        ingredients = 'No ingredients listed.'
                    else:
                        if isinstance(ingredients, str):
                            ingredients = ingredients.replace('[', '').replace(']', '').replace("'", '').replace(';', ',').strip()

                    image_url = row.get('image url')
                    with st.container():
                        left, right = st.columns([1, 3])
                        if pd.notna(image_url) and isinstance(image_url, str) and image_url.strip():
                            left.image(image_url.strip(), width=70)
                        with right:
                            st.markdown(f"""
                                <h4>{perfume_name} <small style="color:gray;">by {brand}</small></h4>
                                <p><strong>Score:</strong> {score}</p>
                            """, unsafe_allow_html=True)

                            if len(description) > 300:
                                with st.expander(" Full Description"):
                                    st.markdown(f"<p>{description}</p>", unsafe_allow_html=True)
                            else:
                                st.markdown(f"<p><strong>Description:</strong> {description}</p>", unsafe_allow_html=True)

                            st.markdown(f"<p><strong>Ingredients:</strong> {ingredients}</p>", unsafe_allow_html=True)

        with col_mol:
            st.subheader("⌬ Explore Molecules Based on Your Preferences")
            for scent in selected_scents:
                if scent not in scent_to_smiles_df.columns:
                    continue

                with st.expander(f"Molecules carrying the scent note '{scent}'", expanded=False):
                    subset = scent_to_smiles_df[scent_to_smiles_df[scent] == 1]
                    mol_entries = []
                    for smi in subset['nonstereosmiles']:
                        m = Chem.MolFromSmiles(smi)
                        if m:
                            img = Draw.MolToImage(m, (120, 120))
                            buf = io.BytesIO(); img.save(buf, 'PNG')
                            img_base64 = base64.b64encode(buf.getvalue()).decode()
                            mol_entries.append({'img_base64': img_base64})

                    if not mol_entries:
                        st.write("No molecules found for this scent.")
                    else:
                        html_content = """
                        <style>
                        .scroll-container { display: flex; flex-wrap: nowrap; overflow-x: auto; gap: 20px; padding: 10px; white-space: nowrap; scrollbar-width: thin; scrollbar-color: #c2b280 #f4eddd; }
                        .scroll-container::-webkit-scrollbar { height: 8px; }
                        .scroll-container::-webkit-scrollbar-track { background: #f4eddd; border-radius: 4px; }
                        .scroll-container::-webkit-scrollbar-thumb { background-color: #c2b280; border-radius: 4px; }
                        </style>
                        <div class="scroll-container">
                        """
                        for mol in mol_entries:
                            html_content += f"""
                        <div style="text-align: center; min-width: 140px;">
                            <img src="data:image/png;base64,{mol['img_base64']}" style="border-radius:8px; box-shadow:0px 4px 6px rgba(0,0,0,0.1);"/>
                        </div>
                        """
                        html_content += "</div>"
                        st.markdown(html_content, unsafe_allow_html=True)