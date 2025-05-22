import os
import sys
import pandas as pd
from rdkit import Chem

from chooseyourperfume.logic_cyp import (
    enrich_with_scent_columns,
    score_perfumes,
    render_molecule,
    ask_preferences,
    get_molecules_for_scents,
    avg_similarity
)

# --- Sample Data ---

sample_perfume_desc_df = pd.DataFrame({
    "name": ["Perfume A", "Perfume B"],
    "description": ["A lovely rose scent", "Woody tones of cedar and sandalwood"]
})

sample_perfume_df = pd.DataFrame({
    "name": ["Perfume A", "Perfume B"],
    "brand": ["Brand A", "Brand B"]
})

sample_scent_to_smiles_df = pd.DataFrame({
    "nonstereosmiles": ["CCO", "CCN"],
    "rose": [1, 0],
    "cedar": [0, 1],
    "sandalwood": [0, 1]
})

# --- Tests ---

def test_enrich_with_scent_columns():
    enriched = enrich_with_scent_columns(sample_perfume_desc_df.copy(), ["rose", "cedar", "sandalwood"])
    assert "rose" in enriched.columns
    assert enriched["rose"].sum() == 1
    assert "cedar" in enriched.columns
    assert enriched["cedar"].sum() == 1

def test_score_perfumes():
    enriched = enrich_with_scent_columns(sample_perfume_desc_df.copy(), ["rose", "cedar"])
    result = score_perfumes(["rose"], enriched, sample_perfume_df)
    assert "score" in result.columns
    assert not result.empty

def test_render_molecule():
    img = render_molecule("CCO")
    assert img is not None

def test_ask_preferences():
    prefs = ask_preferences()
    assert isinstance(prefs, dict)
    assert len(prefs) > 0

def test_get_molecules_for_scents():
    result = get_molecules_for_scents(["rose", "cedar"], sample_scent_to_smiles_df)
    assert not result.empty
    assert "nonstereosmiles" in result.columns

def test_avg_similarity_valid():
    sim = avg_similarity(["CCO", "CCN"])
    assert isinstance(sim, float)
    assert 0.0 <= sim <= 1.0

def test_avg_similarity_insufficient():
    sim = avg_similarity(["CCO"])
    assert sim is None