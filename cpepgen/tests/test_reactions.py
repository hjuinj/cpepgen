"""
This testfile is for function_libs/rdkit/reactions
"""

from cpepgen import chemical_linkages as reactions
from rdkit import Chem
import unittest
class test_reactions(unittest.TestCase):
    def test_make_peptide_bond_0(self):
        gly_smile = "NCC(=O)O"
        gly_mol = Chem.MolFromSmiles(gly_smile)
        # peptide Bond
        mols = reactions.make_peptide_bond(gly_mol, gly_mol)
        self.assertEqual(Chem.MolToSmiles(mols[0], True), "NCC(=O)NCC(=O)O")
        self.assertEqual(Chem.MolToSmiles(mols[1], True), "O")

    # def test_make_peptide_bond_1(self):
    #     gly_smile = "NCC(=O)O"
    #     gly_mol = Chem.MolFromSmiles(gly_smile, sanitize = False)
    #     # peptide Bond
    #     mols = reactions.make_peptide_bond(gly_mol, gly_mol)
    #     self.assertEqual(Chem.MolToSmiles(mols[0], True), "NCC(=O)NCC(=O)O")
    #     self.assertEqual(Chem.MolToSmiles(mols[1], True), "O")

    def test_make_peptide_bond_2(self):
        gly_smile = "NCC(=O)O"
        gly_mol = Chem.MolFromSmiles(gly_smile)
        # peptide Bond
        mols = reactions.make_peptide_bond(gly_mol, gly_mol)
        self.assertEqual(Chem.MolToSmiles(mols[0], True), "NCC(=O)NCC(=O)O")
        self.assertEqual(Chem.MolToSmiles(mols[1], True), "O")

    def test_make_ester_bond(self):
        gly_smile = "NCC(=O)O"
        alcohol_smile = "OCC"
        gly_mol = Chem.MolFromSmiles(gly_smile)
        alcohol_mol = Chem.MolFromSmiles(alcohol_smile)

        # ester Bond
        mols = reactions.make_ester_bond(gly_mol, alcohol_mol)
        self.assertEqual(Chem.MolToSmiles(mols[0], True), "CCOC(=O)CN")
        self.assertEqual(Chem.MolToSmiles(mols[1], True), "O")

    def test_make_disulfide_bond(self):
        cys_smile = "NC(CS)C(=O)O"
        cys_mol = Chem.MolFromSmiles(cys_smile)

        #sulfur Bond
        mols = reactions.make_disulfide_bond(cys_mol, cys_mol)
        self.assertEqual(Chem.MolToSmiles(mols[0], True), "NC(CSSCC(N)C(=O)O)C(=O)O")

if __name__ == "main":
    unittest.main()
