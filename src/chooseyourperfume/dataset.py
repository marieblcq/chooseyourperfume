import os
import pandas as pd
def load_csv(path, sep=",", encoding=None):
    """
    Load a CSV file and return a Dataframe with fallback encodings.
    Args:
        path (str): Path to CSV file.
        sep (str): Separator used in the file, defautlt is ','.
        encoding (str): File encoding.
    Returns:
        pd.DataFrame: the loaded dataset as DataFrame.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"The specified file was not found: {path}")
    try: #First try UTF-8 encoder
        df=pd.read_csv(path, sep=sep, encoding=encoding or "utf-8")
        return df
    except UnicodeDecodeError:
        # Try fallback encodings if utf8 fails
        for enc in ["ISO-8859-1", "cp1252"]:
            try:
                return pd.read_csv(path, sep=sep, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise  # re-raise the last exception if all encodings fail

# LOADING EACH DATASET 

def load_smiles_odors():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "datasets", "Multi-Labelled_Smiles_Odors_dataset.csv"))
    return load_csv(path)

def load_perfume_descriptions():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "datasets", "final_perfume_data.csv"))
    return load_csv(path)

def load_fragrantica_data():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "datasets", "fra_cleaned.csv"))
    return load_csv(path, sep=";")

def load_extended_perfume_set():
    path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "datasets", "fra_perfumes.csv"))
    return load_csv(path)
#def load_smiles_odors():
#    """Dataset 1 - Molecules in SMILES format with corresponding odors."""
#    return load_csv("data/datasets/Multi-Labelled_Smiles_Odors_dataset.csv")

#def load_perfume_descriptions():
#    """Dataset 2 - Commercial perfume descriptions, images, and notes."""
#    return load_csv("data/datasets/final_perfume_data.csv") 

#def load_fragrantica_data():
#    """Dataset 3 - General Fragrantica dataset."""
#    return load_csv("data/datasets/fra_cleaned.csv", sep=";")

#def load_extended_perfume_set():
#    """Dataset 4 - Additional structured perfume records."""
#    return load_csv("data/datasets/fra_perfumes.csv")


#DICTIONNARY 
""" 
Defines scent categories and associated descriptors used for data classification and filtering based on the dataset available.
"""
if __name__ == "__main__":
    df = load_smiles_odors()
# Retrieve all scent columns (lowercased and without extra spaces)
    scent_notes = df.columns[2:].str.strip().str.lower().tolist()

# Display all available scent notes from the dataset
    print("Scent Notes Available in Dataset:")
    for note in sorted(scent_notes):
        print(f"- {note}")

scent_categories = {
    "Floral": ['floral', 'jasmin', 'rose', 'violet', 'lily', 'hyacinth', 'lavender', 'muguet', 'chamomile', 'orangeflower', 'geranium'],
    "Fruity": ['fruity', 'apple', 'apricot', 'banana', 'berry', 'black currant', 'grape', 'grapefruit', 'melon', 'orange', 'peach', 'pear', 'pineapple', 'plum', 'raspberry', 'strawberry', 'tropical', 'cherry', 'lemon', 'juicy', 'ripe', 'fruit skin', 'sour', 'citrus', 'bergamot'],
    "Vegetal / Herbal": ['green', 'grassy', 'herbal', 'leafy', 'celery', 'cucumber', 'hay', 'hawthorn', 'weedy', 'vegetable', 'potato', 'tomato'],
    "Sweet / Gourmand": ['sweet', 'coconut', 'vanilla', 'chocolate', 'caramellic', 'buttery', 'honey', 'almond', 'milky', 'creamy', 'cocoa', 'fatty', 'malty'],
    "Nutty / Seed": ['nutty', 'hazelnut'],
    "Lactonic / Milky": ['lactonic', 'dairy'],
    "Spicy / Aromatic": ['spicy', 'clove', 'cinnamon', 'anisic', 'mint', 'coumarinic', 'camphoreous', 'aromatic', 'balsamic', 'tea', 'terpenic'],
    "Fresh / Volatile": ['fresh', 'cooling', 'ethereal', 'ozone', 'clean', 'aldehydic'],
    "Woody / Resinous": ['woody', 'sandalwood', 'cedar', 'pine', 'cortex', 'amber', 'vetiver', 'tobacco', 'warm'],
    "Animalic / Meaty": ['animal', 'meaty', 'beefy', 'musk', 'leathery', 'sweaty', 'savory', 'fishy'],
    "Smoky / Roasted / Burnt": ['smoky', 'burnt', 'roasted', 'cooked', 'popcorn', 'coffee'],
    "Fermented / Cheesy": ['cheesy', 'fermented'],
    "Sulfurous / Allium-like": ['sulfurous', 'garlic', 'onion', 'alliaceous', 'gassy', 'musty', 'cabbage', 'radish'],
    "Alcohol / Solvent": ['alcoholic', 'brandy', 'winey', 'cognac', 'solvent', 'rummy'],
    "Earthy / Mineral / Metallic": ['earthy', 'mushroom', 'metallic', 'oily', 'powdery', 'orris'],
    "Chemical / Medicinal": ['medicinal', 'ketonic', 'phenolic', 'soapy', 'sharp', 'bitter', 'pungent'],
    "Dry / Bitter / Neutral": ['dry', 'natural', 'odorless', 'waxy']
}