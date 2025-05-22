import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs

#sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chooseyourperfume.chem_analysis import (
    canonicalize_smiles,
    generate_fingerprint,
    assign_primary_category,
    scent_categories,
    run_pipeline
)

import subprocess

# Test 1: SMILES Canonicalization Valid Case
def test_canonicalize_smiles_valid():
    smiles = "CCO"  # Ethanol
    expected = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    assert canonicalize_smiles(smiles) == expected

# Test 2: SMILES Canonicalization Invalid Case
def test_canonicalize_smiles_invalid():
    smiles = "INVALID_SMILES"
    assert canonicalize_smiles(smiles) is None

# Test 3: Fingerprint Length
def test_generate_fingerprint_shape():
    mol = Chem.MolFromSmiles("CCO")
    fp = generate_fingerprint(mol)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    assert len(arr) == 1024

# Test 4: Primary Category Tie Case
def test_primary_category_assignment_tie():
    descriptors = "floral;woody"
    assert assign_primary_category(descriptors) is None

# Test 5: Primary Category Dominant Case
def test_primary_category_assignment_single():
    descriptors = "floral;floral;woody"
    scent_to_category = {scent: category for category, scents in scent_categories.items() for scent in scents}
    assert assign_primary_category(descriptors) == scent_to_category['floral']


def test_pipeline_runs_without_errors():
    # Simply call run_pipeline() under pytest's cwd (the project root),
    # so that load_smiles_odors("data/...") finds the CSV.
    data, tsne = run_pipeline()
    # basic sanity checks
    assert not data.empty, "Pipeline returned an empty DataFrame"
    assert tsne.shape[1] == 2, "t-SNE result must have two dimensions"



# Test 7: No Missing Primary Category After Processing
def test_no_missing_primary_category_after_processing():
    data, _ = run_pipeline()  # Explicitly run the pipeline
    assert data['primary_category'].isna().sum() == 0, "There should be no missing primary categories after processing."
