import os 
import pandas as pd

def load_csv(path, sep=",", encoding=None):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"The specified file was not found: {path}")
    try:
        return pd.read_csv(path, sep=sep, encoding=encoding or "utf-8")
    except UnicodeDecodeError:
        for enc in ["ISO-8859-1", "cp1252"]:
            try:
                return pd.read_csv(path, sep=sep, encoding=enc)
            except UnicodeDecodeError:
                continue
        raise

def load_smiles_odors():
    return load_csv("data/datasets/Multi-Labelled_Smiles_Odors_dataset.csv")

def load_perfume_descriptions():
    return load_csv("data/datasets/final_perfume_data.csv") 

def load_fragrantica_data():
    return load_csv("data/datasets/fra_cleaned.csv", sep=";")

def load_extended_perfume_set():
    return load_csv("data/datasets/fra_perfumes.csv")

def get_scent_categories():
    df = load_smiles_odors()
    scent_notes = df.columns[2:].str.strip().str.lower().tolist()
    print("Scent Notes Available in Dataset:")
    for note in sorted(scent_notes):
        print(f"- {note}")

    return {
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
