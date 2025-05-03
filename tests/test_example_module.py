#from chooseyourperfume.example_module import hello_smiles


# Test the function
#def test_hello_smiles():
    #assert hello_smiles("C(=O)O") == "Hello, C(=O)O!", "Test failed: SMILES input"

from chooseyourperfume import load_dataset
#test data set loaded properly 
if __name__ == "__main__":
    data = load_dataset()
    print(data.head())  # show the first few rows of the dataset --> it worked