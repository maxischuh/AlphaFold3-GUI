import streamlit as st
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

# Function to generate a 3Dmol.js viewer
def render_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        st.error("Invalid SMILES string.")
        return

    # Generate 3D coordinates for the molecule
    mol = Chem.AddHs(mol)  # Add hydrogens
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Convert to MolBlock format
    mol_block = Chem.MolToMolBlock(mol)

    # Initialize Py3Dmol viewer
    viewer = py3Dmol.view()
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})  # Display in stick style
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
        viewer_html = viewer.render().replace('<div', '<div style="width:100%;height:400px"')
        st.components.v1.html(viewer_html, height=400)
