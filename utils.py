import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

from seqs import *


class JSONGenerator:
    """
    A class to generate JSON content for a given set of parameters.
    Attributes:
        name (str): The name of the JSON content.
        model_seeds (list): A list of model seeds.
        sequences (list): A list of sequences, each sequence should have a `to_dict` method.
        bonded_atom_pairs (list, optional): A list of bonded atom pairs. Defaults to None.
        user_ccd (dict, optional): User-defined CCD data. Defaults to None.
        dialect (str, optional): The dialect to be used. Defaults to "alphafold3".
        version (int, optional): The version of the JSON content. Defaults to 2.
    Methods:
        generate_json():
            Generates the JSON content based on the provided attributes.
        _filter_empty(data):
            Recursively filters out empty values from the provided data.
    """
    def __init__(self, name, model_seeds, sequences, bonded_atom_pairs=None, user_ccd=None, dialect="alphafold3", version=1):
        self.name = name
        self.model_seeds = model_seeds
        self.sequences = sequences
        self.bonded_atom_pairs = bonded_atom_pairs
        self.user_ccd = user_ccd
        self.dialect = dialect
        self.version = version
    
    def generate_json(self):
        json_content = {
            "name": self.name,
            "modelSeeds": self.model_seeds,
            "sequences": [self._filter_empty(seq.to_dict()) for seq in self.sequences],
            "dialect": self.dialect,
            "version": self.version,
        }

        if self.bonded_atom_pairs:
            json_content["bondedAtomPairs"] = self.bonded_atom_pairs
        if self.user_ccd:
            json_content["userCCD"] = self.user_ccd

        return json_content

    def _filter_empty(self, data):
        if isinstance(data, dict):
            return {k: self._filter_empty(v) for k, v in data.items() if v not in [None, "", [], {}]}
        elif isinstance(data, list):
            return [self._filter_empty(item) for item in data if item not in [None, "", [], {}]]
        else:
            return data


class CcdConverter:
    """
    A class to convert and manipulate chemical structures using SMILES strings and RDKit.

    Attributes:
        smiles (str): The SMILES string of the molecule.
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object.
        final_mol (rdkit.Chem.rdchem.Mol): The final RDKit molecule object after processing.
        final_labels (dict): The final atom labels after processing.
        target_atom_idx (list): The indices of target atoms.
        viewer (py3Dmol.view): The 3D viewer for molecule visualization.

    Methods:
        parse_smiles(smiles):
            Parses a SMILES string and returns an RDKit molecule object.
        
        add_hydrogens(mol):
            Adds hydrogen atoms to the molecule.
        
        generate_3d_coordinates(mol):
            Generates 3D coordinates for the molecule and optimizes its geometry.
        
        generate_atom_labels(mol):
            Generates labels for each atom in the molecule.
        
        identify_leaving_group_by_bond(mol, atom1_idx, atom2_idx):
            Identifies the leaving group by breaking a bond between two atoms.
        
        remove_atoms(mol, atom_indices):
            Removes specified atoms from the molecule.
        
        visualize_molecule(mol, atom_labels, highlight_atom_idx=None):
            Visualizes the molecule with optional highlighted atoms.
        
        first_step(smiles):
            Performs the first step of processing: parsing SMILES, adding hydrogens, generating 3D coordinates, and visualizing the molecule.
        
        second_step(atom1_idx, atom2_idx):
            Performs the second step of processing: identifying the leaving group by breaking a bond and visualizing the updated molecule.
        
        third_step(leaving_group_atoms):
            Performs the third step of processing: removing the leaving group atoms and visualizing the final molecule.
        
        fourth_step(compound_id):
            Performs the fourth step of processing: creating mmCIF fields for the final molecule.
        
        create_mmcif_fields(compound_id, mol, atom_labels):
            Creates mmCIF fields for the given molecule and atom labels.
    """
    def __init__(self):
        self.smiles = None
        self.mol = None
        self.final_mol = None
        self.final_labels = None
        self.target_atom_idx = None
        self.viewer = None

    def parse_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        Chem.SanitizeMol(mol)
        return mol

    def add_hydrogens(self, mol):
        return Chem.AddHs(mol)

    def generate_3d_coordinates(self, mol):
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.UFFOptimizeMolecule(mol)
        AllChem.AssignStereochemistryFrom3D(mol)

        return mol

    def generate_atom_labels(self, mol):
        labels = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            labels[idx] = f"{symbol}{idx}"  # Example: C13
        return labels

    def identify_leaving_group_by_bond(self, mol, atom1_idx, atom2_idx):
        # Break the bond and generate fragments
        rw_mol = Chem.RWMol(mol)
        rw_mol.RemoveBond(atom1_idx, atom2_idx)
        fragmented_mol = rw_mol.GetMol()

        # Get fragments with original atom indices
        fragments = Chem.GetMolFrags(fragmented_mol, asMols=False, sanitizeFrags=True)

        # If only one fragment, it's a cycle
        if len(fragments) == 1:
            return [], fragmented_mol

        # Identify the smaller fragment as the leaving group
        fragment_sizes = [len(fragment) for fragment in fragments]
        smaller_fragment_idx = fragment_sizes.index(min(fragment_sizes))
        leaving_group_atoms = list(fragments[smaller_fragment_idx])

        return leaving_group_atoms, fragmented_mol

    def remove_atoms(self, mol, atom_indices):
        editable_mol = Chem.EditableMol(mol)
        for idx in sorted(atom_indices, reverse=True):  # Reverse sort to handle index shifts
            editable_mol.RemoveAtom(idx)
        return editable_mol.GetMol()

    def visualize_molecule(self, mol, atom_labels, highlight_atom_idx=None):
        mol_block = Chem.MolToMolBlock(mol)
        viewer = py3Dmol.view(width=800, height=400)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {}})

        # Add labels for all atoms
        for idx, label in atom_labels.items():
            viewer.addLabel(
                label,
                {"fontSize": 10, "fontColor": "black", "backgroundColor": "rgba(255,255,255,0)"},
                {"serial": idx},
            )

        # Highlight selected atoms if provided
        if highlight_atom_idx is not None:
            for idx in highlight_atom_idx:
                viewer.setStyle(
                    {"serial": idx},
                    {"stick": {}, "sphere": {"radius": 0.5, "color": "rgba(255,0,0,0.8)"}},
                )

        viewer.zoomTo()
        return viewer._make_html()

    def first_step(self, smiles):
        self.smiles = smiles

        mol = self.parse_smiles(smiles)
        mol = self.add_hydrogens(mol)
        mol = self.generate_3d_coordinates(mol)

        self.mol = mol

        # Step 1: Visualize and choose bond to break
        atom_labels = self.generate_atom_labels(mol)
        return self.visualize_molecule(mol, atom_labels)

    def second_step(self, atom1_idx, atom2_idx):
        leaving_group_atoms, updated_mol = self.identify_leaving_group_by_bond(self.mol, atom1_idx, atom2_idx)

        if not leaving_group_atoms:
            print("No leaving group identified (cyclic structure detected).")
        else:
            print(f"Identified atoms in the leaving group: {leaving_group_atoms}")

        # Visualize updated molecule after bond removal
        updated_labels = self.generate_atom_labels(updated_mol)  # Use updated molecule for visualization
        self.final_mol = updated_mol
        self.final_labels = updated_labels
        return leaving_group_atoms, self.visualize_molecule(
            updated_mol, updated_labels, highlight_atom_idx=leaving_group_atoms
        )

    def third_step(self, leaving_group_atoms):
        final_mol = self.remove_atoms(self.final_mol, leaving_group_atoms)
        final_labels = self.generate_atom_labels(final_mol)
        self.final_mol = final_mol
        self.final_labels = final_labels
        return self.visualize_molecule(final_mol, final_labels)

    def fourth_step(self, compound_id):
        cif_lines = self.create_mmcif_fields(compound_id, self.final_mol, self.generate_atom_labels(self.final_mol))
        return str("\n".join(cif_lines))

    def create_mmcif_fields(self, compound_id, mol, atom_labels):
        atoms = mol.GetAtoms()
        bonds = mol.GetBonds()
        conformer = mol.GetConformer()

        cif_lines = [
            f"data_{compound_id}",
            "#",
            f"_chem_comp.id {compound_id}",
            f"_chem_comp.name '{compound_id}'",
            f"_chem_comp.type non-polymer",
            f"_chem_comp.formula '{compound_id}'",
            "_chem_comp.mon_nstd_parent_comp_id ?",
            "_chem_comp.pdbx_synonyms ?",
            "_chem_comp.formula_weight ?",
            "#",
        ]

        # Atom loop
        cif_lines.extend([
            "loop_",
            "_chem_comp_atom.comp_id",
            "_chem_comp_atom.atom_id",
            "_chem_comp_atom.type_symbol",
            "_chem_comp_atom.charge",
            "_chem_comp_atom.model_Cartn_x",
            "_chem_comp_atom.model_Cartn_y",
            "_chem_comp_atom.model_Cartn_z",
            "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
            "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
            "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
        ])
        for atom in atoms:
            idx = atom.GetIdx()
            pos = conformer.GetAtomPosition(idx)
            cif_lines.append(
                f"{compound_id} {atom_labels[idx]} {atom.GetSymbol()} 0 "
                f"{pos.x:.3f} {pos.y:.3f} {pos.z:.3f} "
                f"{pos.x:.3f} {pos.y:.3f} {pos.z:.3f} "
            )

        # Bond loop
        cif_lines.extend([
            "#",
            "loop_",
            "_chem_comp_bond.comp_id",
            "_chem_comp_bond.atom_id_1",
            "_chem_comp_bond.atom_id_2",
            "_chem_comp_bond.value_order",
        ])
        for bond in bonds:
            start_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType()
            bond_order = "SINGLE" if bond_type == Chem.BondType.SINGLE else \
                         "DOUBLE" if bond_type == Chem.BondType.DOUBLE else \
                         "TRIPLE" if bond_type == Chem.BondType.TRIPLE else "AROMATIC"
            cif_lines.append(
                f"{compound_id} {atom_labels[start_atom]} {atom_labels[end_atom]} {bond_order}"
            )

        cif_lines.extend([
            "#",
            "_pdbx_chem_comp_descriptor.type SMILES_CANONICAL",
            f"_pdbx_chem_comp_descriptor.descriptor '{Chem.MolToSmiles(mol)}'",
            "#",
        ])

        return cif_lines
