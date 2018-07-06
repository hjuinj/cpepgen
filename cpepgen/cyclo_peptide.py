import copy
from itertools import groupby
from rdkit import Chem
from rdkit.Chem import AllChem
from functools import reduce
from .utils import *
#from Bio.SeqUtils import seq1, seq3
"""
    TODOs:
        - scheme for residueName column names that reflects chemical modifications done to the residues (e.g. methylated/D-amino acid)
        - change 3code to 1code
        - should I keep the internal self.peptide always as RWMol()? and when user wants to then call self.peptide.GetMol()
        - support for charged residues
        - PDB block temperature factor is -nan for residues bigger than 1
        - !sort atom by residue number (done)
        - make methylate(i.e. add new terminal atoms) more general (done)
"""
class Peptide():
    """
    """
    def __init__(self, sequence , methylation_sites = [], protonation_states = {}) :
        """
        Parameters
        ----------
        sequence : str
            one-letter code of the amino acide sequence
        methylation_sites : list (optional)
            indicates which residues are methylated at backbone
            Indexing starts from 1, maximum index should be less than or equal to len(sequence)
        protonation_states : dict (optional)
            indicates which residues has non-standard(standard states has number 0) protonation states, e.g. {1:1} means res 1 has the first protonation state (which is usually the sole charged states
            Only used for registering the residue name column
            Histidine requres special assignment of which state is which
        """
        methylation_sites = set(methylation_sites)
        if len(methylation_sites) > 0  and max(methylation_sites) > len(sequence):
            raise ValueError("The number of methylation site is not valid")
        self.methylation_sites = methylation_sites
        self.protonation_states = {i : protonation_states[i] if i in protonation_states else 0 for i in range(len(sequence))}
        print(self.protonation_states)

        self.sequence = sequence

        self.numHs2AddCondition = lambda atm: (atm.GetNumImplicitHs() + atm.GetNumExplicitHs()) if (atm.GetIsAromatic() or atm.GetAtomicNum() != 6) else 0

        residues = [self.process_residue(r, idx) for idx,r in enumerate(sequence)]
        residues = [self.n_methylate(res) if (idx+1) in methylation_sites else res for idx, res in enumerate(residues)]
        self.peptide = reduce(make_peptide_bond, residues) #FIXME


    def process_residue(self, res, idx):
        """Deals with discrepencies between rdkit atom names and GROMOS atom names

        Parameters
        ----------
        res : rdkmol
            one amino acid residue

        """
        res = Chem.MolFromSequence(res, flavor = 1)
        ######## odd case of atom name differences
        if res.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName().strip() == "ILE":
            for atm in res.GetAtoms():
                if atm.GetPDBResidueInfo().GetName().strip() == "CD1" :
                    atm.GetPDBResidueInfo().SetName(self._name_column_formatter("CD", atomname = True))

        #######################################
        #########  residue name encoding into three letter codes
        #########  first letter is the one-letter code of amino acid
        ########   the second and third letter combined to give a
        ########   decimal number that encodes various properties
        # decimal number is obtained from a string of binary bits:
        # each of the bits stand for:
        # (_ _ _ _ _:currently left unused)(_:r or l residue) (_:n-methylated or not) (_ _:the state of the amino acide, usually zero
        #
        resname = self.sequence[idx].upper()
        code = ""
        #isomer
        if self.sequence[idx].isupper(): #R isomer
            code +="0"
        else: #L isomer
            code +="1"
        #N-methylation
        if idx+1 in self.methylation_sites:
            code += "1"
        else:
            code += "0"
        # protonation states
        if self.protonation_states[idx] == 0:
            code += "00"
        elif self.protonation_states[idx] == 1:
            code += "01"
        elif self.protonation_states[idx] == 2:
            code += "10"
        elif self.protonation_states[idx] == 3:
            code += "11"
        resname += str(int(code, 2)).zfill(2)

        for atm in res.GetAtoms():
            atm.GetPDBResidueInfo().SetResidueName(self._name_column_formatter(resname, residuename = True))

        return res


    def _atom_info_helper(self, atom):
        """returns an rdkit AtomInfo object for a given atom
        """
        return Chem.AtomPDBResidueInfo(atomName = "TODO",
        residueName = atom.GetPDBResidueInfo().GetResidueName(),
        residueNumber = atom.GetPDBResidueInfo().GetResidueNumber(),
        chainId = atom.GetPDBResidueInfo().GetChainId())

    def _name_column_formatter(self, name, atomname = False, residuename = False):
        """produces the correct number of spaces for the PDB format atom name column
        """
        return "{}{: <3s}".format("" if residuename is True else " ",name)

    def n_methylate(self, res = None):
        """methylate the N atom in a residues

        Parameters
        ----------
        res : rdkmol
            if is None, modification is done on the self.peptide sequence
        """
        if res is None:
            res = self.peptide
        return add_atom_at_site(res, "N", "CN")

    def gromosff_hydrogenate(self, res = None):
        """add gromos type hydrogens to the peptide
        the number of hydrogens to add to an identified atom site can be > 1
        this is not necessary for gromos as pdb2g96 can do the same thing

        Parameters
        ----------
        res : rdkmol
            if is None, modification is done on the self.peptide sequence
        """
        def atomInfo(atm):
            info = self._atom_info_helper(atm)
            name = 'H' + atm.GetPDBResidueInfo().GetName().strip()[1:] #duplicated names will be changed later
            name = self._name_column_formatter(name)
            info.SetName(name)
            return info

        def atomicNumber(atm):
            return 1

        if res is None:
            res = self.peptide

        self.peptide =  add_terminal_atom(condition = self.numHs2AddCondition,
            mol = res,
            newAtomicNumber =  atomicNumber,
            newAtomInfo = atomInfo,
            aromatic = lambda atm : False,
            bondType = 1)
        self.gromos_deduplicate_names()
        self.reorder_atoms()
        return self.peptide

    def gromos_deduplicate_names(self):
        """All duplicated atomnames per residue are appended by ascending integers (e.g. if multiple atoms have atomName "HD" in a residue they are renamed to HD1, HD2...)
        """
        atoms = groupby(self.peptide.GetAtoms(), lambda atm : (atm.GetPDBResidueInfo().GetResidueNumber(), atm.GetPDBResidueInfo().GetName().strip()))
        for key, atmgroup in atoms:
            atmlist = list(atmgroup)
            if len(atmlist) > 1:
                for i, atm in enumerate(atmlist):
                    self.peptide.GetAtomWithIdx(atm.GetIdx()).GetPDBResidueInfo().SetName(self._name_column_formatter(key[1]+str(i + 1)))

        #assertion block to ensure deduplication above do not introduce new duplicate names
        atoms = groupby(self.peptide.GetAtoms(), lambda atm : (atm.GetPDBResidueInfo().GetResidueNumber(), atm.GetPDBResidueInfo().GetName().strip()))
        for key, atmgroup in atoms:
            atmlist = list(atmgroup)
            assert len(atmlist) == 1, "New Duplicates have been introduced"
        return

    def reorder_atoms(self):
        """change index of the atoms to ensure atoms are ordered by ascending residue number
        """
        order = [(i.GetPDBResidueInfo().GetName().strip(), i.GetPDBResidueInfo().GetResidueNumber()) for i in self.peptide.GetAtoms()] # currently the atom name is not used in sorting
        order = [i[0] for i in sorted(enumerate(order), key  = lambda x : x[1][1])]
        self.peptide = Chem.RenumberAtoms(self.peptide, order)
        return

    def get_pdb(self):
        return Chem.MolToPDBBlock(self.peptide)

    def get_smiles(self):
        # return Chem.MolToSmiles(self.peptide, isomericSmiles = True, allBondsExplicit = True)
        return Chem.MolToSmiles(self.peptide, isomericSmiles = True, allBondsExplicit = True, allHsExplicit = True)


    def gen_conf(self, n = 1):
        """Generate user specified number of conformers using ETKDG
        """
        self.nConf= n
        self.peptide = Chem.RWMol(self.peptide)
        hydrogenateSites = [key for key, val in enumerate(map(self.numHs2AddCondition, self.peptide.GetAtoms())) if val == 0]

        self.peptide = Chem.AddHs(self.peptide, onlyOnAtoms = hydrogenateSites)
        # self.peptide.UpdatePropertyCache()
        print(Chem.MolToPDBBlock(self.peptide))
        AllChem.EmbedMultipleConfs(self.peptide, numConfs = n)
        self.peptide = remove_terminal_atom(self.peptide, lambda atm : atm.GetPDBResidueInfo() is None)

        return

class CycloPeptide(Peptide):
    """
    Assumptions:
        - all input mol/smiles/Peptide obj are LINEAR
    """
    def __init__(self, sequence , methylation_sites = [], protonation_states = []) :
        super(CycloPeptide, self).__init__(sequence, methylation_sites, protonation_states)
        self.peptide = make_peptide_bond(self.peptide)

"""
io.StringIO usage
https://github.com/pandegroup/openmm/blob/956df4160da8ad3bdd0629c1b5e1f29d56be671e/wrappers/python/tests/TestForceField.py#L265
bypass saving to file and retreving
"""
