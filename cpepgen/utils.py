import copy
from rdkit import Chem
from rdkit.Chem import AllChem

_aa321 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

_bondtypes = {1: Chem.BondType.SINGLE,
              1.5: Chem.BondType.AROMATIC,
              2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE,
              4: Chem.BondType.QUADRUPLE,
              5: Chem.BondType.QUINTUPLE,
              6: Chem.BondType.HEXTUPLE,
              7: Chem.BondType.ONEANDAHALF,}

def make_peptide_bond(mol, res = None, start = "N", end = "C", delete = "OXT"):
    """Performs one condesation rxn between a molecule and residue/itself
    default creates peptide bond

    Parameters
    ----------
        mol : rdkmol
            Main molecule on which reaction is performed
        res : rdkmol
            None or a single residue, when it is None, self-condensation is performed on @mol
        start : str
            atom name of one of the two atoms to which connection is established
        end : str
            atom name of the other atom to which connection is established
        delete : str
            default to hydroxy oxygen thats eliminated during condensation

    Returns
    -------
    mol : rdkmol
        modified molecule
    """
    startIdx, endIdx, deleteIdx  = -1, -1, -1
    for idx, atm in enumerate(mol.GetAtoms()):  # get the last occurence of end and delete atomname
        if atm.GetPDBResidueInfo().GetName().strip() == end:
            endIdx = idx
        elif atm.GetPDBResidueInfo().GetName().strip() == delete:
            deleteIdx = idx

    if res is not None: #residue addition
        lastResNum = -1
        for idx, atm in enumerate(mol.GetAtoms()):
            lastResNum = atm.GetPDBResidueInfo().GetResidueNumber()
        lastResNum += 1
        for idx, atm in enumerate(res.GetAtoms()): #get the last occurence of start atomname
            atm.GetPDBResidueInfo().SetResidueNumber(lastResNum)
            if atm.GetPDBResidueInfo().GetName().strip() == start:
                startIdx = idx
        startIdx += mol.GetNumAtoms()
        mol = Chem.CombineMols(mol, res)
    else: #cyclisation
        for idx, atm in enumerate(mol.GetAtoms()): #get the first occurence of start atomname
            if atm.GetPDBResidueInfo().GetName().strip() == start:
                startIdx = idx
                break

    mol = Chem.RWMol(mol)
    mol.AddBond(startIdx, endIdx, Chem.BondType.SINGLE)
    mol.RemoveAtom(deleteIdx)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

def remove_terminal_atom(mol, condition):
    """Removes atom(s) on the molecule based on expression

    Parameters
    ----------
        mol : rdkmol
            Main molecule on which deletion is performed
        condition : function
            returns boolean values for each atom in @mol, those with True are deleted

    Returns
    -------
    mol : rdkmol
        modified molecule
    """
    mol = Chem.RWMol(mol)
    todelete = [atm.GetIdx() for atm in filter(condition, copy.deepcopy(mol).GetAtoms())]
    for idx in reversed(todelete):
        mol.RemoveAtom(idx)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

#more general than add_atom_at_site()
def add_terminal_atom(mol, condition, newAtomicNumber, newAtomInfo = lambda atm : Chem.AtomMonomerInfo(), aromatic = lambda atm: False, bondType = 1): #TODO make bondType also a function?
    """create new atom(s) based on expression

    Parameters:
    ----------
    mol : rdkmol
        Main moelcule one which addition is perforemd
    condition : function
        outputs how many atoms to add at each atom site of @mol
    newAtomicNumber : function
        outputs element types to add at each atom site of @mol
    newAtomInfo : function
        outputs the rdkit AtomInfo block information add at each atom site of @mol
    aromatic : function
        True/False value for atom aromaticity
    bondType : float

    Returns
    -------
    mol : rdkmol
        modified molecule

    Example:
    --------
    residue = Chem.MolFromSequence("KK", flavor = 0)
    residue = add_terminal_atom(mol = residue,
        condition = lambda atm: atm.GetPDBResidueInfo().GetName().strip() == "N",
        newAtomName = lambda atm : "CN",
        newAtomicNumber = lambda atm: 6)
    print(Chem.MolToPDBBlock(residue))
    """
    def mapper(condition, mol):
        for idx, repeats in enumerate(map(condition, mol.GetAtoms())):
            for j in range(repeats):
                yield mol.GetAtomWithIdx(idx)


    bondType = _bondtypes[bondType]
    mol = Chem.RWMol(mol)
    for atom in mapper(condition, copy.deepcopy(mol)):
        newatom = Chem.Atom(newAtomicNumber(atom))
        newatom.SetIsAromatic(aromatic(atom))
        newatom.SetMonomerInfo(newAtomInfo(atom))
        # newatom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName = " {: <3s}".format(newAtomName(atom)),
        # residueName = atom.GetPDBResidueInfo().GetResidueName(),
        # residueNumber = atom.GetPDBResidueInfo().GetResidueNumber(),
        # chainId = atom.GetPDBResidueInfo().GetChainId(),
        # ))

        mol.AddAtom(newatom)
        mol.AddBond(atom.GetIdx(), mol.GetNumAtoms() - 1, bondType)
    # mol.UpdatePropertyCache(strict = False)
    mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()

def add_atom_at_site(mol, site, name, atomicNumber = 6, aromatic = False, bondType = 1):
    """
    ! Assumes the molecule has PDBResidueInfo()
    """
    bondType = _bondtypes[bondType]
    tmpmol = copy.deepcopy(mol)
    mol = Chem.RWMol(mol)
    for atom in tmpmol.GetAtoms():
        # print(atom.GetPDBResidueInfo().GetName())
        if atom.GetPDBResidueInfo().GetName().strip() == site:
            atom.SetIsAromatic(aromatic)
            atom.SetAtomicNum(atomicNumber)
            atom.GetPDBResidueInfo().SetName(" {: <3s}".format(name))
            mol.AddAtom(atom)
            mol.AddBond(atom.GetIdx(), mol.GetNumAtoms() - 1, bondType)
            mol.UpdatePropertyCache()
    Chem.GetSSSR(mol)
    return mol.GetMol()
