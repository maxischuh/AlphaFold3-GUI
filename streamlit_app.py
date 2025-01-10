import streamlit as st
import string
import json
import numpy as np


def protein_entity_ui(identifier):
    st.markdown(f"**Protein | Entity {identifier}**")
    copies = st.number_input(
        "Copies", min_value=1, max_value=16, value=1, step=1, key=f"protein_copies_{identifier}"
    )
    sequence = st.text_area(
        "Protein Sequence", height=100, key=f"protein_sequence_{identifier}"
    )
    seed = st.number_input(
        "Seed", value=42, step=1, key=f"protein_seed_{identifier}"
    )
    return copies, sequence, seed

def ligand_entity_ui(identifier):
    st.markdown(f"**Ligand | Entity {identifier}**")
    copies = st.number_input(
        "Copies", min_value=1, max_value=16, value=1, step=1, key=f"ligand_copies_{identifier}"
    )
    sequence = st.text_area(
        "Ligand SMILES", height=50, key=f"ligand_sequence_{identifier}"
    )
    seed = st.number_input(
        "Seed", value=42, step=1, key=f"ligand_seed_{identifier}"
    )
    bond_image = st.file_uploader(
        "Upload Bond Image", type=["png", "jpg", "jpeg"], key=f"bond_image_{identifier}"
    )
    bond_string = st.text_input(
        "Bond Information", key=f"bond_string_{identifier}"
    )
    return copies, sequence, seed, bond_image, bond_string

def rna_entity_ui(identifier):
    st.markdown(f"**RNA | Entity {identifier}**")
    copies = st.number_input(
        "Copies", min_value=1, max_value=16, value=1, step=1, key=f"rna_copies_{identifier}"
    )
    sequence = st.text_area(
        "RNA Sequence", height=100, key=f"rna_sequence_{identifier}"
    )
    seed = st.number_input(
        "Seed", value=42, step=1, key=f"rna_seed_{identifier}"
    )
    return copies, sequence, seed

def dna_entity_ui(identifier):
    st.markdown(f"**DNA | Entity {identifier}**")
    copies = st.number_input(
        "Copies", min_value=1, max_value=16, value=1, step=1, key=f"dna_copies_{identifier}"
    )
    sequence = st.text_area(
        "DNA Sequence", height=100, key=f"dna_sequence_{identifier}"
    )
    seed = st.number_input(
        "Seed", value=42, step=1, key=f"dna_seed_{identifier}"
    )
    return copies, sequence, seed

st.title("AlphaFold3 JSON GUI")

num_proteins = st.number_input(
    "Number of Protein Entities", min_value=0, value=1, step=1, key="num_proteins"
)
num_ligands = st.number_input(
    "Number of Ligand Entities", min_value=0, value=0, step=1, key="num_ligands"
)
num_dna = st.number_input(
    "Number of DNA Entities", min_value=0, value=0, step=1, key="num_dna"
)
num_rna = st.number_input(
    "Number of RNA Entities", min_value=0, value=0, step=1, key="num_rna"
)

def generate_ui_from_numbers(num_proteins, num_ligands, num_dna, num_rna):
    entity_outputs = []
    entity_identifiers = []

    for idx in range(int(num_proteins)):
        identifier = string.ascii_uppercase[idx]
        group = protein_entity_ui(identifier)
        entity_outputs.append(group)
        entity_identifiers.append(identifier)

    for idx in range(int(num_ligands)):
        identifier = string.ascii_uppercase[idx + int(num_proteins)]
        group = ligand_entity_ui(identifier)
        entity_outputs.append(group)
        entity_identifiers.append(identifier)

    for idx in range(int(num_dna)):
        identifier = string.ascii_uppercase[idx + int(num_proteins) + int(num_ligands)]
        group = dna_entity_ui(identifier)
        entity_outputs.append(group)
        entity_identifiers.append(identifier)

    for idx in range(int(num_rna)):
        identifier = string.ascii_uppercase[
            idx + int(num_proteins) + int(num_ligands) + int(num_dna)
        ]
        group = rna_entity_ui(identifier)
        entity_outputs.append(group)
        entity_identifiers.append(identifier)

    return entity_outputs, entity_identifiers

entity_outputs, entity_identifiers = generate_ui_from_numbers(
    num_proteins, num_ligands, num_dna, num_rna
)

def format_data(entity_outputs, entity_identifiers):
    result = {}
    for identifier, entity_group in zip(entity_identifiers, entity_outputs):
        entity_data = {}
        entity_keys = ["copies", "sequence", "seed", "bond_image", "bond_string"]
        for key, value in zip(entity_keys, entity_group):
            if key == "bond_image" and value is not None:
                entity_data[key] = value.name
            else:
                entity_data[key] = value
        result[identifier] = entity_data
    return json.dumps(result, indent=4)

if st.button("Generate JSON"):
    json_output = format_data(entity_outputs, entity_identifiers)
    st.text_area("AlphaFold3 Config", json_output, height=400)
