""" 
These tests check the coherence between the dictionary `scent_categories` used to classify scents 
and the available loaded data. 
"""

from chooseyourperfume.dataset import scent_categories, load_smiles_odors

# Retrieve all scent columns
df = load_smiles_odors()
scent_notes = df.columns[2:].str.strip().str.lower().tolist()

# Display all available scent notes from the dataset
print("Scent Notes Available in Dataset:")
for note in sorted(scent_notes):
    print(f"- {note}")

# Check for duplicate scent columns in the dataset loaded
if len(scent_notes) != len(set(scent_notes)):
    print("\n Duplicate scent columns found in the dataset!")
else:
    print("\n No duplicate scent columns found.")

def test_scent_categories_is_dict():
    """Ensure scent_categories is a dictionary."""
    assert isinstance(scent_categories, dict), "scent_categories should be a dictionary."

def test_unique_scent_notes():
    """Ensure each scent note appears only once across all categories."""
    all_notes = [note for notes in scent_categories.values() for note in notes]
    duplicates = set([note for note in all_notes if all_notes.count(note) > 1])
    assert not duplicates, f"Duplicate notes found in scent_categories: {duplicates}"

def test_dataset_scent_columns_in_scent_categories():
    """Verify that all scent-related columns in the dataset are covered by scent_categories."""
    df = load_smiles_odors()
    scent_columns = df.columns[2:].str.strip().str.lower().tolist()  # Skip first 2 non-scent columns
    dict_notes = [note for notes in scent_categories.values() for note in notes]
    missing_notes = set(scent_columns) - set(dict_notes)
    assert not missing_notes, f"Notes present in dataset but missing in scent_categories: {missing_notes}"

def test_scent_categories_notes_exist_in_dataset():
    """Ensure all notes in scent_categories are actually present in the dataset columns."""
    df = load_smiles_odors()
    scent_columns = df.columns[2:].str.strip().str.lower().tolist()
    dict_notes = [note for notes in scent_categories.values() for note in notes]
    unused_notes = set(dict_notes) - set(scent_columns)
    assert not unused_notes, f"Notes in scent_categories not found in dataset columns: {unused_notes}"

def test_scent_categories_are_not_empty():
    """Check that each scent category contains at least one scent note."""
    empty_categories = [category for category, notes in scent_categories.items() if not notes]
    assert not empty_categories, f"The following categories are empty: {empty_categories}"

def test_scent_categories_notes_are_strings():
    """Verify all scent notes in scent_categories are of type string."""
    invalid_notes = []
    for category, notes in scent_categories.items():
        for note in notes:
            if not isinstance(note, str):
                invalid_notes.append((category, note, type(note)))
    assert not invalid_notes, f"The following notes are not strings: {invalid_notes}"
