"""
This file contains several reactions for RDKIT

NOTES:
    - Check what happens with ASP ASP ester bond formation? need more specifications to do that side targets?
    - check_reaction_mols, also checking how many of the reaction sites are present in one molecule? - (more specific smart: make_aa_backbone_bond = [N][C][O:5]=[C:2][OH:1].[N:3][CH2:4][C][=O][C]>>[O:5]=[C:2][N:3][CH2:4].[OH2:1] "
    - Check in unittest also fail states
authors: Shuzhe Wang & Benjamin Schroeder
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import typing as t

#GERNEAL Functions:
def check_reaction_mols(reactants: list, reaction_smarts: str) -> bool:
    """
    This function checks if all requirement in the Educts are fullfilled for the reaction

    :param reactants: list of Chem.Mol used as reactants
    :param reaction_smarts: str defining the reaction
    :return: boolean true if everything is fullfilled
    """
    educts = [educt for educt in reaction_smarts.split(">>")[0].split(".")]
    educts_mol = [AllChem.MolFromSmarts(educt) for educt in educts]

    #TODO check Mol is None

    if (len(reactants) > len(educts)):
        raise ValueError("There are more rectants than expected for this reaction.\n Given reactants: "+str(len(reactants))+" expected reactants: "+str(len(educts)))
    elif (len(reactants) < len(educts)):
        raise ValueError("There are more "+str(len(educts))+" rectants expected for this reaction.\n But only "+str(len(reactants))+" reactants were given")
    else:
        for ind, reactant in enumerate(reactants):
            if not (reactant.HasSubstructMatch(educts_mol[ind])):
                raise ValueError("Reactant " + str(ind) + " with smile " + Chem.MolToSmiles(
                    reactants[ind]) + " has not the required reaction pattern: " + str(educts[ind]))
        return True
    return False

def do_reaction(reactants: list, reaction_smarts: str) -> t.List[Chem.Mol]:
    """
    This function is executed to perform a reaction
    :param reactants:
    :param reaction_smarts:
    :return: list[Chem.Mol] products
    """
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    except Exception as err:
        raise ValueError("Could not generate reaction object in do_reaction function:\n\t"+str("\t".join(err.args)))
    # reac1
    ps = rxn.RunReactants(tuple(reactants))
    if (ps == None):
        raise Exception("Reaction did not happen.")

    return ps[0]

def perform_generic_reaction(reactants:t.List[Chem.Mol], reaction_smarts:str):
    """
    this  wrapper is preforming an reaction with
    :param reactants:
    :param reaction_smarts:
    :return:
    """
    try:
        check_reaction_mols(reactants=reactants, reaction_smarts=reaction_smarts)
        products = do_reaction(reactants=reactants,
                               reaction_smarts=reaction_smarts)
    except Exception as err:
        raise ValueError("Could not perform Reaction.\n Failed in molecule checking Step\n\t" + str("\t".join(err.args)))
    except ValueError as err:
        raise ValueError("Could not perform Reaction.\n Failed in reaction Step:\n\t" + str("\t".join(err.args)))
    return products

#REACTIONS:

def make_custom_bond(reactant1, reactant2, reaction_smarts):
    products = perform_generic_reaction(reactants=[reactant1, reactant2], reaction_smarts=reaction_smarts)
    return products

def make_peptide_bond(carboxylGroupContaining: Chem.Mol, nitrogenGroupContaining: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smarts  = {
        "educts" : [
        "[C:1](=[O])[OH:2]",
        "[NH2,NH1,NH3+1:3]-[CH2,CH1:4]" ,
        ],
        "products" : [
        "[C:1](=[O])[N:3]-[C:4]",
        "[OH2:2]",
        ],
    }
    reaction_smarts = [".".join(reaction_smarts[i]) for i in reaction_smarts]
    reaction_smarts = ">>".join(reaction_smarts)
    # reaction_smarts = '[O:5]=[C:2][OH:1].[N:3][CH2,CH1:4]>>[O:5]=[C:2][N:3][CH2:4].[OH2:1]'

    products = perform_generic_reaction(reactants=[carboxylGroupContaining, nitrogenGroupContaining], reaction_smarts=reaction_smarts)
    # return products, reaction_smarts
    return products

def make_amide_bond(carboxylGroupContaining: Chem.Mol, nitrogenGroupContaining: Chem.Mol) -> (t.List[Chem.Mol], str):
    return make_peptide_bond(carboxylGroupContaining, nitrogenGroupContaining)

def make_ester_bond(carboxylGroupContaining: Chem.Mol, alcoholGroupContaining: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smarts = '[O:5]=[C:2][OH:1].[HO:3][CH2:4]>>[O:5]=[C:2][O:3][CH2:4].[OH2:1]'
    products = perform_generic_reaction(reactants=[carboxylGroupContaining, alcoholGroupContaining], reaction_smarts=reaction_smarts)
    # return products, reaction_smarts
    return products

def make_disulfide_bond(thiolGroupContaining1: Chem.Mol, thiolGroupContaining2: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smarts = '[C:1][SH:2].[HS:3][C:4]>>[C:1][S:2][S:3][C:4]'
    products = perform_generic_reaction(reactants=[thiolGroupContaining1, thiolGroupContaining2], reaction_smarts=reaction_smarts)
    # return products, reaction_smarts
    return products

def make_disulphide_bond(thiolGroupContaining1: Chem.Mol, thiolGroupContaining2: Chem.Mol) -> (t.List[Chem.Mol], str):
    return make_disulfide_bond(thiolGroupContaining1, thiolGroupContaining2)

a = Chem.MolFromSequence("A")
b = Chem.MolFromSequence("C")
# print(Chem.MolToPDBBlock(b))

out= make_amide_bond(a,b)
mol = out[0]
# Chem.CalcExplicitValence(mol)
mol.UpdatePropertyCache(strict = False)

print(Chem.MolToPDBBlock(mol), Chem.MolToSmiles(mol))
