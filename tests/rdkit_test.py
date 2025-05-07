from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image

smiles = "O=C(O)Cc1ccccc1"  # Benzoic acid
mol = Chem.MolFromSmiles(smiles)

if mol:
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    img.show()  # Opens the molecule image
else:
    print("Invalid SMILES")
