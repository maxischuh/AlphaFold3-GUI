import streamlit as st
import json

from seqs import *
from utils import *

def main():
    """
    Main function that sets up the AlphaFold3 GUI.
    It collects user inputs, processes them, and generates JSON for AlphaFold3.
    """
    # Setting page configuration for Streamlit
    st.set_page_config(
        page_title="AlphaFold3 GUI",
        layout="centered",
        initial_sidebar_state="auto",
    )
    
    # Title and brief explanation for the UI
    st.title("AlphaFold3 GUI")
    st.write("This app allows you to input data for various entities and generates a structured JSON compatible with AlphaFold3.")
    st.write("Use the sections below to specify your job details, entity sequences, and options for covalent binding. "
             "Covalent binding is handled step-by-step by defining leaving groups and binding sites, which streamlines the process.")

    # -- Initialize converter in session_state so it persists across re-runs --
    if "converter" not in st.session_state:
        st.session_state["converter"] = CcdConverter()

    # User input for job name and random seeds
    st.header("Job Details")
    job_name = st.text_input("Job Name", value="MyAlphaFoldJob", key="job_name")
    model_seeds_text = st.text_area("Model Seeds (comma-separated integers)", value="1, 2", key="model_seeds")
    model_seeds = [int(seed.strip()) for seed in model_seeds_text.split(",") if seed.strip().isdigit()]

    # User input for the number of entities
    st.header("Number of Entities")
    num_proteins = st.number_input("Number of Protein Entities", min_value=0, step=1, value=0, key="num_proteins")
    num_dna = st.number_input("Number of DNA Entities", min_value=0, step=1, value=0, key="num_dna")
    num_rna = st.number_input("Number of RNA Entities", min_value=0, step=1, value=0, key="num_rna")
    num_ligands = st.number_input("Number of Ligand Entities", min_value=0, step=1, value=0, key="num_ligands")

    st.write("---")

    # Global counter for unique entity identifiers: A, B, C, ...
    entity_counter = iter(range(65, 91))  # Up to 'Z'

    # ------------------------------------------------------------------------
    # Helper Functions to Generate Input Fields
    # ------------------------------------------------------------------------
    def generate_protein_inputs(count):
        """
        Collects protein-specific user input, such as sequences, modifications,
        and MSA details. Returns a list of ProteinSequence objects.
        """
        entity_data = []
        for i in range(count):
            num_copies = st.number_input(
                f"Number of Copies for Protein {i+1}",
                min_value=1, step=1, value=1,
                key=f"protein_copies_{i}"
            )
            entity_labels = [chr(next(entity_counter)) for _ in range(num_copies)]

            st.subheader(f"Protein {i+1} (labels: {', '.join(entity_labels)})")
            sequence = st.text_input(f"Protein Sequence {i+1}", key=f"protein_sequence_{i}")
            modifications = st.text_area(
                f"Modifications (JSON array) for Protein {i+1}",
                key=f"protein_modifications_{i}",
                value="[]"
            )
            try:
                modifications = json.loads(modifications)
            except json.JSONDecodeError:
                st.error(f"Invalid JSON for modifications of Protein {i+1}.")
                modifications = []

            unpaired_msa = st.text_input(f"Unpaired MSA for Protein {i+1}", key=f"protein_msa_{i}")
            unpaired_msa_path = st.text_input(f"Unpaired MSA Path for Protein {i+1}", key=f"protein_msa_path_{i}")
            paired_msa = st.text_input(f"Paired MSA for Protein {i+1}", key=f"protein_paired_msa_{i}")
            paired_msa_path = st.text_input(f"Paired MSA Path for Protein {i+1}", key=f"protein_paired_msa_path_{i}")
            templates = st.text_area(
                f"Templates (JSON array) for Protein {i+1}",
                key=f"protein_templates_{i}",
                value="[]"
            )
            try:
                templates = json.loads(templates)
            except json.JSONDecodeError:
                st.error(f"Invalid JSON for templates of Protein {i+1}.")
                templates = []

            entity_data.append(ProteinSequence(
                entity_labels,
                sequence,
                modifications,
                unpaired_msa if unpaired_msa else None,
                unpaired_msa_path if unpaired_msa_path else None,
                paired_msa if paired_msa else None,
                paired_msa_path if paired_msa_path else None,
                templates
            ))
        return entity_data

    def generate_rna_inputs(count):
        """
        Collects RNA-specific user input, such as sequences, modifications,
        and MSA details. Returns a list of RNASequence objects.
        """
        entity_data = []
        for i in range(count):
            num_copies = st.number_input(
                f"Number of Copies for RNA {i+1}",
                min_value=1, step=1, value=1,
                key=f"rna_copies_{i}"
            )
            entity_labels = [chr(next(entity_counter)) for _ in range(num_copies)]

            st.subheader(f"RNA {i+1} (labels: {', '.join(entity_labels)})")
            sequence = st.text_input(f"RNA Sequence {i+1}", key=f"rna_sequence_{i}")
            modifications = st.text_area(
                f"Modifications (JSON array) for RNA {i+1}",
                key=f"rna_modifications_{i}",
                value="[]"
            )
            try:
                modifications = json.loads(modifications)
            except json.JSONDecodeError:
                st.error(f"Invalid JSON for modifications of RNA {i+1}.")
                modifications = []

            unpaired_msa = st.text_input(f"Unpaired MSA for RNA {i+1}", key=f"rna_msa_{i}")
            unpaired_msa_path = st.text_input(f"Unpaired MSA Path for RNA {i+1}", key=f"rna_msa_path_{i}")

            entity_data.append(RNASequence(
                entity_labels,
                sequence,
                modifications,
                unpaired_msa if unpaired_msa else None,
                unpaired_msa_path if unpaired_msa_path else None
            ))
        return entity_data

    def generate_dna_inputs(count):
        """
        Collects DNA-specific user input, such as sequences and modifications.
        Returns a list of DNASequence objects.
        """
        entity_data = []
        for i in range(count):
            num_copies = st.number_input(
                f"Number of Copies for DNA {i+1}",
                min_value=1, step=1, value=1,
                key=f"dna_copies_{i}"
            )
            entity_labels = [chr(next(entity_counter)) for _ in range(num_copies)]

            st.subheader(f"DNA {i+1} (labels: {', '.join(entity_labels)})")
            sequence = st.text_input(f"DNA Sequence {i+1}", key=f"dna_sequence_{i}")
            modifications = st.text_area(
                f"Modifications (JSON array) for DNA {i+1}",
                key=f"dna_modifications_{i}",
                value="[]"
            )
            try:
                modifications = json.loads(modifications)
            except json.JSONDecodeError:
                st.error(f"Invalid JSON for modifications of DNA {i+1}.")
                modifications = []

            entity_data.append(DNASequence(
                entity_labels,
                sequence,
                modifications
            ))
        return entity_data

    # ------------------------------------------------------------------------
    # Generates user input fields for ligands, including optional covalent binding
    # ------------------------------------------------------------------------
    def generate_ligand_inputs(count):
        """
        Collects ligand-specific user input, including SMILES data and,
        optionally, covalent bond setup. Returns a list of dictionaries containing
        ligands, bond data, and user-defined CCD data.
        """
        entity_data = []

        # Initialize or get references from st.session_state so we don't lose data on each re-run
        if "ligand_bonds" not in st.session_state:
            st.session_state["ligand_bonds"] = {}   # keyed by i
        if "ligand_userccds" not in st.session_state:
            st.session_state["ligand_userccds"] = {}  # keyed by i

        # Collect all labels for potential binding
        all_entity_labels = []
        for prot_obj in proteins:
            all_entity_labels.extend([lbl for lbl in prot_obj.id])
        for rna_obj in rna:
            all_entity_labels.extend([lbl for lbl in rna_obj.id])
        for dna_obj in dna:
            all_entity_labels.extend([lbl for lbl in dna_obj.d])

        for i in range(count):
            bond = st.session_state["ligand_bonds"].get(i, None)
            userccd = st.session_state["ligand_userccds"].get(i, None)

            num_copies = st.number_input(
                f"Number of Copies for Ligand {i+1}",
                min_value=1, step=1, value=1,
                key=f"ligand_copies_{i}"
            )
            entity_labels = [chr(next(entity_counter)) for _ in range(num_copies)]
            st.subheader(f"Ligand {i+1} (labels: {', '.join(entity_labels)})")

            smiles = st.text_input(f"Ligand SMILES {i+1}", key=f"ligand_smiles_{i}")

            # Button to render the molecule visually
            if st.button(f"Render Molecule {i+1}"):
                viewer_html = st.session_state["converter"].first_step(smiles)
                if viewer_html:
                    st.components.v1.html(viewer_html, height=450)

            covalent_bond = st.checkbox(
                f"Enable Covalent Binding for Ligand {i+1}",
                value=(bond is not None),
                key=f"covalent_{i}"
            )

            # If covalent binding is enabled, additional fields for customizing bonds are displayed
            if covalent_bond:
                st.write("### Covalent Binding Options")
                with st.expander("Covalent Binding Workflow", expanded=True):
                    if (len(all_entity_labels) > 0) & (count == 1):
                        tab_define, tab_list = st.tabs(["By Indices", "Full Atom List"])

                        with tab_define:
                            st.write("#### Define Leaving Group by Atom Indices")
                            st.write("Please select the atoms between which the bond should be cleaved.")
                            atom_idx1 = st.number_input(
                                f"Ligand {i+1}: Atom Index 1",
                                min_value=0, step=1, value=0,
                                key=f"atom_idx1_{i}"
                            )
                            atom_idx2 = st.number_input(
                                f"Ligand {i+1}: Atom Index 2",
                                min_value=0, step=1, value=1,
                                key=f"atom_idx2_{i}"
                            )

                            if st.button(f"Split Bond at Indices for Ligand {i+1}"):
                                leaving_group atoms, viewer_html = st.session_state["converter"].second_step(atom_idx1, atom_idx2)
                                if viewer_html:
                                    st.components.v1.html(viewer_html, height=450)

                        with tab_list:
                            st.write("#### Define Leaving Group by Listing Atom Indices")
                            st.write("Please list all atoms which should be cleaved off.")
                            leaving_group_atoms = st.text_input(
                                f"Leaving Group Atoms (comma-separated) for Ligand {i+1}",
                                value="",
                                key=f"leaving_group_atoms_{i}"
                            )
                            leaving_group_atoms = [
                                int(idx.strip())
                                for idx in leaving_group_atoms.split(",") if idx.strip().isdigit()
                            ]
                            if st.button(f"Remove Group from Atom List {i+1}"):
                                viewer_html = st.session_state["converter"].third_step(leaving_group_atoms)
                                if viewer_html:
                                    st.components.v1.html(viewer_html, height=450)

                        st.write("#### Bond Details")
                        st.write("##### Ligand")
                        target_atom = st.selectbox(
                            f"Select Atom Index for Ligand {i+1}",
                            st.session_state["converter"].final_labels.values(),
                            key=f"ligand_atom_idx_{i}"
                        )
                        # Split out label and index
                        target_atom_label = "".join([j for j in target_atom if not j.isdigit()])
                        target_atom_idx = int("".join([j for j in target_atom if j.isdigit()]))

                        st.write("##### Protein/RNA/DNA")
                        entity_to_bind = st.selectbox(
                            f"Select Entity to Bind for Ligand {i+1}",
                            all_entity_labels,
                            key=f"entity_bind_{i}"
                        )

                        entity_seq_idx_selection = None
                        for prot in proteins:
                            if entity_to_bind == prot.id[0]:
                                entity_seq_selection = prot.sequence
                                entity_seq_idx_selection = [
                                    f"{aa}{idx+1}" for idx, aa in enumerate(entity_seq_selection)
                                ]

                        entity_residue = st.selectbox(
                            f"Select Residue of {entity_to_bind}",
                            entity_seq_idx_selection,
                        )
                        entity_atom_label = "".join([j for j in entity_residue if not j.isdigit()])
                        entity_atom_idx = int("".join([j for j in entity_residue if j.isdigit()]))

                        entity_atom_name = st.text_input(
                            f"CCD Atom Name for {entity_to_bind}",
                            key=f"entity_atom_name_{i}"
                        )

                        if st.button(f"Finalize Covalent Bond {i+1}"):
                            userccd = st.session_state["converter"].fourth_step(entity_labels[0])
                            bond = [
                                [entity_to_bind, entity_atom_idx, entity_atom_name],
                                [entity_labels[0], 1, f"{target_atom_label}{target_atom_idx}"]
                            ]
                            st.session_state["ligand_bonds"][i] = bond
                            st.session_state["ligand_userccds"][i] = userccd
                            st.success("Covalent bond data finalized.")
                    else:
                        st.warning(
                            "No protein/RNA/DNA entities defined to bind to or multiple ligands. Covalent binding disabled."
                        )
            else:
                st.session_state["ligand_bonds"][i] = None
                st.session_state["ligand_userccds"][i] = None
                bond = None
                userccd = None

            # Construct the final data structure for this ligand
            ligand = LigandSequence(
                id=entity_labels,
                ccd_codes=entity_labels,
                smiles=None if covalent_bond else smiles
            )

            entity_data.append({
                "ligand": ligand,
                "bond": bond,
                "userccd": userccd
            })

        return entity_data

    # ------------------------------------------------------------------------
    # Generate user input fields for all entity types
    # ------------------------------------------------------------------------
    st.header("Input Entity Details")
    proteins = generate_protein_inputs(num_proteins)
    rna = generate_rna_inputs(num_rna)
    dna = generate_dna_inputs(num_dna)
    ligands = generate_ligand_inputs(num_ligands)

    print(ligands)

    st.write("---")

    # ------------------------------------------------------------------------
    # Generate JSON output
    # ------------------------------------------------------------------------
    st.header("Generated JSON Output")
    if st.button("Generate JSON"):
        print(ligands)

        output_json = JSONGenerator(
            name=job_name,
            model_seeds=model_seeds,
            sequences=(
                proteins
                + rna
                + dna
                + [l["ligand"] for l in ligands]  # gather just the 'LigandSequence' parts
            ),
            bonded_atom_pairs=(
                [l["bond"] for l in ligands] if num_ligands > 0 else None
            ),
            user_ccd=(
                [l["userccd"] for l in ligands if l["userccd"] is not None]
                if num_ligands > 0 else None
            )
        ).generate_json()

        # Display JSON in the app
        st.json(output_json)

        # Provide download option for the generated JSON
        json_string = json.dumps(output_json, indent=4)
        st.download_button(
            label="Download JSON",
            data=json_string,
            file_name=f"alphafold3_config_{job_name.lower()}.json",
            mime="application/json"
        )

if __name__ == "__main__":
    main()
