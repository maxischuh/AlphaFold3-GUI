class Sequence:
    """
    Base class for all sequence types in AlphaFold 3.
    """
    def __init__(self, id, sequence_type):
        self.id = id
        self.sequence_type = sequence_type

    def to_dict(self):
        raise NotImplementedError("Subclasses must implement the to_dict method.")


class ProteinSequence(Sequence):
    def __init__(self, id, sequence, modifications=None, unpaired_msa=None, unpaired_msa_path=None, paired_msa=None, paired_msa_path=None, templates=None):
        super().__init__(id, "protein")
        self.sequence = sequence
        self.modifications = modifications or []
        self.unpaired_msa = unpaired_msa
        self.unpaired_msa_path = unpaired_msa_path
        self.paired_msa = paired_msa
        self.paired_msa_path = paired_msa_path
        self.templates = templates or []

    def to_dict(self):
        return {
            "protein": {
                "id": self.id,
                "sequence": self.sequence,
                "modifications": self.modifications,
                "unpairedMsa": self.unpaired_msa,
                "unpairedMsaPath": self.unpaired_msa_path,
                "pairedMsa": self.paired_msa,
                "pairedMsaPath": self.paired_msa_path,
                "templates": self.templates,
            }
        }


class RNASequence(Sequence):
    def __init__(self, id, sequence, modifications=None, unpaired_msa=None, unpaired_msa_path=None):
        super().__init__(id, "rna")
        self.sequence = sequence
        self.modifications = modifications or []
        self.unpaired_msa = unpaired_msa
        self.unpaired_msa_path = unpaired_msa_path

    def to_dict(self):
        return {
            "rna": {
                "id": self.id,
                "sequence": self.sequence,
                "modifications": self.modifications,
                "unpairedMsa": self.unpaired_msa,
                "unpairedMsaPath": self.unpaired_msa_path,
            }
        }


class DNASequence(Sequence):
    def __init__(self, id, sequence, modifications=None):
        super().__init__(id, "dna")
        self.sequence = sequence
        self.modifications = modifications or []

    def to_dict(self):
        return {
            "dna": {
                "id": self.id,
                "sequence": self.sequence,
                "modifications": self.modifications,
            }
        }


class LigandSequence(Sequence):
    def __init__(self, id, ccd_codes=None, smiles=None):
        super().__init__(id, "ligand")
        self.ccd_codes = ccd_codes
        self.smiles = smiles

    def to_dict(self):
        ligand_dict = {"ligand": {"id": self.id}}
        if self.ccd_codes:
            ligand_dict["ligand"]["ccdCodes"] = self.ccd_codes
        if self.smiles:
            ligand_dict["ligand"]["smiles"] = self.smiles
        return ligand_dict
