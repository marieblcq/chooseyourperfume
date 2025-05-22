import os
import sys
import pandas as pd
from chooseyourperfume.dataset import load_csv


def test_load_dataset_1_smiles_odors():
    path = os.path.join("data", "datasets", "Multi-Labelled_Smiles_Odors_dataset.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "nonStereoSMILES" in df.columns
    assert "descriptors" in df.columns


def test_load_dataset_2_perfume_descriptions():
    path = os.path.join("data", "datasets", "final_perfume_data.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Name" in df.columns
    assert "Notes" in df.columns


def test_load_dataset_3_fragrantica_general():
    path = os.path.join("data", "datasets", "fra_cleaned.csv")
    df = load_csv(path=path, sep=";")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Perfume" in df.columns
    assert "Top" in df.columns


def test_load_dataset_4_fragrantica_clean():
    path = os.path.join("data", "datasets", "fra_perfumes.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Name" in df.columns
    assert "Description" in df.columns
