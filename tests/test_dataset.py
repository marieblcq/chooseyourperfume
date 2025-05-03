import os
import sys
import pandas as pd

# Ensure the src/ directory is in the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from chooseyourperfume.dataset import load_csv


def test_load_dataset_1_smiles_odors():
    path = os.path.join("data", "datasets", "Multi-Labelled_Smiles_Odors_dataset.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "SMILES" in df.columns
    assert "Odor" in df.columns


def test_load_dataset_2_perfume_descriptions():
    path = os.path.join("data", "datasets", "Perfume_Descriptions.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Name" in df.columns
    assert "Notes" in df.columns


def test_load_dataset_3_fragrantica_general():
    path = os.path.join("data", "datasets", "Fragrantica_General.csv")
    df = load_csv(path=path, sep=";")
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Perfume" in df.columns or "Name" in df.columns
    assert "Top" in df.columns


def test_load_dataset_4_fragrantica_clean():
    path = os.path.join("data", "datasets", "Fragrantica_Clean.csv")
    df = load_csv(path=path)
    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert "Name" in df.columns
    assert "Description" in df.columns
