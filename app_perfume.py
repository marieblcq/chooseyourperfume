import streamlit as st
import io, random, base64
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from src.chooseyourperfume.logic_cyp import (
    load_data, ask_preferences, score_perfumes, get_molecules_for_scents, avg_similarity
)

# --- Config & CSS ---
st.set_page_config(page_title="Choose Your Perfume", layout="wide")

st.markdown("""
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
""", unsafe_allow_html=True)

# --- Data ---
@st.cache_data
def cached_load_data():
    return load_data()

perfume_to_scent_df, perfume_clean_df, perfume_df, scent_to_smiles_df = cached_load_data()

# --- Header ---
col1, col2 = st.columns([1, 4])
with col1:
    st.image("assets/logo.png", width=180)
with col2:
    st.markdown("<h1> ⋆ CHOOSE YOUR PERFUME ⋆ </h1><p style='color:gray;'>Find your signature scent</p>", unsafe_allow_html=True)

# --- Step 1: Preferences ---
st.header("1. Tell us about the scents you love")
scent_dict = ask_preferences()
categories = list(scent_dict.keys())

if st.button("🎲 Surprise Me!"):
    all_notes = [note for sublist in scent_dict.values() for note in sublist]
    picks = random.sample(all_notes, min(3, len(all_notes)))
    for cat in categories:
        st.session_state[f"sel_{cat}"] = [n for n in picks if n in scent_dict[cat]]
    st.success("✨ Surprise picks loaded!")

for cat in categories:
    st.session_state.setdefault(f"sel_{cat}", [])

col_a, col_b = st.columns(2)
for idx, cat in enumerate(categories):
    target = col_a if idx % 2 == 0 else col_b
    with target.expander(cat, expanded=False):
        st.multiselect(
            " ",
            scent_dict[cat],
            key=f"sel_{cat}"
        )

selected_scents = list(set(note for cat in categories for note in st.session_state.get(f"sel_{cat}", [])))
st.write("**Your picks:**", ", ".join(selected_scents) if selected_scents else "None yet.")

# --- Step 2: Weights ---
weights = {}
if selected_scents:
    st.subheader("2. How much do you love each scent?")
    slider_labels = {0.1: "It's okay", 0.3: "I like", 0.7: "I love", 1.0: "I adore", 1.5: "Obsessed!"}
    for note in selected_scents:
        weights[note] = st.select_slider(
            f"{note}",
            options=list(slider_labels.keys()),
            value=0.7,
            format_func=lambda x: slider_labels[x],
            key=f"w_{note}"
        )

# --- Step 3: Generate Recommendations ---
if st.button("🔍 Generate Recommendations"):
    if not selected_scents:
        st.warning("🚨 Please pick at least one note.")
    else:
        with st.spinner("🔬 Finding your perfect perfumes..."):
            top = score_perfumes(selected_scents, perfume_to_scent_df, perfume_df, weights)

        col_perf, col_mol = st.columns(2, gap="large")

        with col_perf:
            st.subheader("⚭ Perfume Matches")
            if top.empty:
                st.warning("🚫 No matches found.")
            else:
                for _, row in top.head(5).iterrows():
                    perfume_name = str(row.get('name_x') or row.get('name_y') or row.get('name', 'Unknown')).title()
                    brand = row.get('brand_x') or row.get('brand', 'Unknown')
                    score = row.get('score', 0)
                    description = row.get('description_x') or row.get('description_y') or 'No description available.'
                    notes = row.get('main accords_x') or row.get('notes') or row.get('main accords_y') or 'No ingredients listed.'
                    image_url = row.get('image url_x') or row.get('image url')

                    if isinstance(notes, str):
                        notes = notes.replace('[', '').replace(']', '').replace(';', ',').replace("'", '').strip()

                    with st.container():
                        left, right = st.columns([1, 3])
                        if image_url and isinstance(image_url, str):
                            left.image(image_url.strip(), width=70)
                        with right:
                            st.markdown(f"<h4>{perfume_name} <small style='color:gray;'>by {brand}</small></h4><p><strong>Score:</strong> {score:.2f}%</p>", unsafe_allow_html=True)
                            st.markdown(f"<p><strong>Ingredients:</strong> {notes}</p>", unsafe_allow_html=True)
                        max_chars = 150
                        if len(description) > max_chars:
                            short_desc = description[:max_chars].rsplit(' ', 1)[0] + "..."
                            st.markdown(f"<p><strong>Description:</strong> {short_desc}</p>", unsafe_allow_html=True)
                            with st.expander("Read more"):
                                st.markdown(f"<p>{description}</p>", unsafe_allow_html=True)
                        else:
                            st.markdown(f"<p><strong>Description:</strong> {description}</p>", unsafe_allow_html=True)

        st.markdown(f"<p><strong>Ingredients:</strong> {notes}</p>", unsafe_allow_html=True)
        with col_mol:
            st.subheader("⌬ Explore Molecules Based on Your Preferences")
            for scent in selected_scents:
                if scent not in scent_to_smiles_df.columns:
                    continue
                with st.expander(f"Molecules for scent '{scent}'"):
                    subset = scent_to_smiles_df[scent_to_smiles_df[scent] == 1]
                    smiles_list = subset["nonstereosmiles"].tolist()
                    mol_entries = []
                    for smi in smiles_list:
                        m = Chem.MolFromSmiles(smi)
                        if m:
                            img = Draw.MolToImage(m, (120, 120))
                            buf = io.BytesIO(); img.save(buf, 'PNG')
                            img_base64 = base64.b64encode(buf.getvalue()).decode()
                            mol_entries.append(f'<img src="data:image/png;base64,{img_base64}" style="border-radius:8px; box-shadow:0px 4px 6px rgba(0,0,0,0.1);"/>')

                    if mol_entries:
                        st.markdown(f"""
                        <div class="scroll-container" style="display:flex;overflow-x:auto;gap:20px;padding:10px;">
                        {''.join(f'<div style="min-width:140px;text-align:center;">{mol}</div>' for mol in mol_entries)}
                        </div>
                        """, unsafe_allow_html=True)
                        avg_sim = avg_similarity(smiles_list)
                        st.markdown(f"**Average Tanimoto similarity:** {avg_sim:.2f}" if avg_sim else "**Average Tanimoto similarity:** N/A (need ≥2 molecules)")
                    else:
                        st.write("No molecules found for this scent.")

            st.markdown("""
            <details style="background-color:#fff;padding:12px;border-radius:8px;margin-top:8px;border:1px solid #ddd;">
              <summary style="cursor:pointer;font-weight:600;">🔍 What is Tanimoto similarity?</summary>
              <div style="margin-top:8px;line-height:1.5;">
                Tanimoto similarity compares two molecular fingerprints—it's like comparing barcodes of chemical features.
                <ul>
                  <li><strong>1.0</strong> → perfect match</li>
                  <li><strong>0.0</strong> → no shared features</li>
                </ul>
                Use it to gauge how alike molecules really are! For example if the tanimoto simlarity is of 0.7, the moelcules are 70% alike ! 
              </div>
            </details>
            """, unsafe_allow_html=True)

            st.markdown("""
        <details style="background-color:#fff;padding:12px;border-radius:8px;margin-top:8px;border:1px solid #ddd;">
        <summary style="cursor:pointer;font-weight:600;">💡 How are your perfume matches scored?</summary>
        <div style="margin-top:8px;line-height:1.5;">
            Your perfume recommendations are scored based on how well each scent in a perfume matches your personal preferences.
            You selected the notes you like and used sliders to show how much you love each one — from "just okay" to "obsessed."
            <br><br>
            Each perfume earns points for containing your selected notes, and the more important a note is to you (based on the slider), the more it contributes.
            We total these points and calculate a <strong>percentage match</strong> — so a perfume that includes all your top notes could score close to <strong>100%</strong>,
            while one that only matches a few will score lower. The higher the percentage, the better the match for your personal scent profile.
        </div>
        </details>
        """, unsafe_allow_html=True)
