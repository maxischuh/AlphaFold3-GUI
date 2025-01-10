import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

# Function to generate a 3Dmol.js viewer
def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        st.error("Invalid SMILES string.")
        return None

    # Generate 3D coordinates for the molecule
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Convert to MolBlock format
    mol_block = Chem.MolToMolBlock(mol)

    # Initialize Py3Dmol viewer
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

# Streamlit app setup
st.title("SMILES to Molecule Viewer")

# Input SMILES string
smiles_input = st.text_input("Enter a SMILES string:", "CCO")

# Render molecule button
if st.button("Render Molecule"):
    viewer = render_molecule(smiles_input)
    if viewer:
        # Generate HTML and JavaScript for the viewer
        viewer_html = viewer._make_html()
        st.components.v1.html(viewer_html, height=450)
