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
    return load_csv("data/datasets/fra_cleaned.csv", sep=";")

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