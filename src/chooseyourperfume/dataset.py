import os
import pandas as pd
def load_csv(path, sep=",", encoding=None):
    """

    Generic loader for CSV files.

    Args:
        path (str): File path to CSV.
        sep (str): Separator used in the file.
        encoding (str): File encoding.

    Returns:
        pd.DataFrame: the loaded dataset as DataFrame.

    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"The specified file was not found: {path}")
    try: #first try UTF-8 encoder
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

# Loading each dataset 
#test

def load_smiles_odors(): # Dataset 1 - molecules in SMILE format + corresponding odors with various descriptors
    return load_csv("data/datasets/Multi-Labelled_Smiles_Odors_dataset.csv")

def load_perfume_descriptions():  # Dataset 2 - Commercial descriptions + image + notes
    return load_csv("data/datasets/final_perfume_data.csv") 

def load_fragrantica_data():  # Dataset 3 - General Fragrantica set
    return load_csv("data/datasets/fra_cleaned.csv")

def load_extended_perfume_set():  # Dataset 4 - Additional structured perfume records
    return load_csv("data/datasets/fra_perfumes.csv")

# This block runs only if the script is executed directly
if __name__ == "__main__":
    datasets = {
        "Smiles & Odors": load_smiles_odors,
        "Perfume Descriptions": load_perfume_descriptions,
        "Fragrantica General Data": load_fragrantica_data,
        "Extended Perfumes": load_extended_perfume_set,
    }

    for name, loader in datasets.items():
        try:
            df = loader()
            print(f"\n{name} dataset loaded successfully, containing({len(df)} rows):")
            print(df.head(2)) # Display the first 2 rows for verification of the data loading
        except Exception as e:
            print(f" Error loading {name}: {e}")


#Filter
# Example loading step (adjust if already loaded)
df = load_smiles_odors()

# Filter rows where the 'fishy' category is active (i.e., value is 1)
fishy_df = df[df["fishy"] == 1]

# Display result
#print(fishy_df[["nonStereoSMILES", "descriptors"]])  # Show relevant columns

def filter_by_category(df, category):
    if category not in df.columns:
        raise ValueError(f"Category '{category}' not found in dataset.")
    return df[df[category] == 1]

# Example usage
category = "fruity"
fruity_df = filter_by_category(df, category)
#print(fruity_df[["nonStereoSMILES", "descriptors"]])


# Création de la liste des catégories (toutes les colonnes sauf les deux premières)
categories = df.columns[2:].tolist()

# Vérification
print(categories)

scent_categories = {
    "Floral": [
        'floral', 'jasmin', 'rose', 'violet', 'lily', 'hyacinth', 'lavender',
        'muguet', 'chamomile', 'orangeflower', 'geranium'
    ],
    "Fruity": [
        'fruity', 'apple', 'apricot', 'banana', 'berry', 'black currant',
        'citrus', 'grape', 'grapefruit', 'melon', 'orange', 'peach', 'pear',
        'pineapple', 'plum', 'raspberry', 'strawberry', 'tropical'
    ],
    "Vegetal / Herbal": [
        'green', 'grassy', 'herbal', 'leafy', 'celery', 'cucumber', 'hay',
        'hawthorn', 'weedy', 'vegetable'
    ],
    "Sweet / Gourmand": [
        'sweet', 'creamy', 'coconut', 'vanilla', 'chocolate', 'caramellic',
        'buttery', 'honey', 'almond', 'milky'
    ],
    "Nutty / Seed": [
        'nutty', 'hazelnut', 'almond'
    ],
    "Lactonic / Milky": [
        'lactonic', 'dairy', 'milky', 'creamy'
    ],
    "Spicy / Aromatic": [
        'spicy', 'clove', 'cinnamon', 'anisic', 'mint', 'coumarinic', 'camphoreous'
    ],
    "Fresh / Volatile": [
        'fresh', 'cooling', 'ethereal', 'ozone', 'clean', 'citrus'
    ],
    "Woody / Resinous": [
        'woody', 'sandalwood', 'cedar', 'pine', 'cortex', 'amber', 'vetiver'
    ],
    "Animalic / Meaty": [
        'animal', 'meaty', 'beefy', 'musk', 'leathery', 'sweaty', 'savory'
    ],
    "Smoky / Roasted / Burnt": [
        'smoky', 'burnt', 'roasted', 'cooked', 'popcorn', 'coffee'
    ],
    "Fermented / Cheesy": [
        'cheesy', 'fermented', 'dairy'
    ],
    "Sulfurous / Allium-like": [
        'sulfurous', 'garlic', 'onion', 'alliaceous', 'gassy', 'musty',
        'cabbage', 'radish'
    ],
    "Alcohol / Solvent": [
        'alcoholic', 'brandy', 'winey', 'cognac', 'solvent', 'rummy'
    ],
    "Earthy / Mineral / Metallic": [
        'earthy', 'mushroom', 'metallic', 'oily', 'powdery'
    ],
    "Chemical / Medicinal": [
        'medicinal', 'ketonic', 'phenolic', 'soapy', 'sharp', 'bitter'
    ],
    "Dry / Bitter / Neutral": [
        'dry', 'bitter', 'natural', 'odorless'
    ]
}
