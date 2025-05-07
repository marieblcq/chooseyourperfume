from rdkit import Chem
from rdkit.Chem import Draw

smiles = "O=C(O)Cc1ccccc1"  # Benzoic acid
mol = Chem.MolFromSmiles(smiles)
img = Draw.MolToImage(mol)
img.show()  # Should open an image window if Pillow is working

import sys
import platform
import rdkit
import PIL

st.write("Python Executable:", sys.executable)
st.write("Platform:", platform.platform())
st.write("RDKit Version:", rdkit.__version__)
st.write("Pillow Version:", PIL.__version__)
