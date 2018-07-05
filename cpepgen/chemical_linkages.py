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
def check_reaction_mols(reactants: list, reaction_smart: str) -> bool:
    """
    This function checks if all requirement in the Educts are fullfilled for the reaction

    :param reactants: list of Chem.Mol used as reactants
    :param reaction_smart: str defining the reaction
    :return: boolean true if everything is fullfilled
    """
    educts = [educt for educt in reaction_smart.split(">>")[0].split(".")]
    educts_mol = [AllChem.MolFromSmarts(educt) for educt in educts]

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

def do_reaction(reactants: list, reaction_smart: str) -> t.List[Chem.Mol]:
    """
    This function is executed to perform a reaction
    :param reactants:
    :param reaction_smart:
    :return: list[Chem.Mol] products
    """
    try:
        rxn = AllChem.ReactionFromSmarts(reaction_smart)
    except Exception as err:
        raise ValueError("Could not generate reaction object in do_reaction function:\n\t"+str("\t".join(err.args)))
    # reac1
    ps = rxn.RunReactants(tuple(reactants))
    if (ps == None):
        raise Exception("Reaction did not happen.")

    return ps[0]

def generic_reaction(reactants:t.List[Chem.Mol], reaction_smart:str):
    """
    this  wrapper is preforming an reaction with
    :param reactants:
    :param reaction_smart:
    :return:
    """
    try:
        check_reaction_mols(reactants=reactants, reaction_smart=reaction_smart)
        products = do_reaction(reactants=reactants,
                               reaction_smart=reaction_smart)
    except Exception as err:
        raise ValueError("Could not perform Reaction.\n Failed in reaction Step:\n\t" + str("\t".join(err.args)))
    except ValueError as err:
        raise ValueError("Could not perform Reaction.\n Failed in molecule checking Step\n\t" + str("\t".join(err.args)))
    return products

#REACTIONS:
def make_peptide_bond(caboxylGroupContaining: Chem.Mol, nitrogenGroupContaining: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smart = '[O:5]=[C:2][OH:1].[N:3][CH2:4]>>[O:5]=[C:2][N:3][CH2:4].[OH2:1]'
    products = generic_reaction(reactants=[caboxylGroupContaining, nitrogenGroupContaining], reaction_smart=reaction_smart)
    return products, reaction_smart

def make_ester_bond(caboxylGroupContaining: Chem.Mol, alcoholGroupContaining: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smart = '[O:5]=[C:2][OH:1].[HO:3][CH2:4]>>[O:5]=[C:2][O:3][CH2:4].[OH2:1]'
    products = generic_reaction(reactants=[caboxylGroupContaining, alcoholGroupContaining], reaction_smart=reaction_smart)
    return products, reaction_smart

def make_sulfur_bond(tihiolGroupContaining1: Chem.Mol, tihiolGroupContaining2: Chem.Mol) -> (t.List[Chem.Mol], str):
    reaction_smart = '[C:1][SH:2].[HS:3][C:4]>>[C:1][S:2][S:3][C:4]'
    products = generic_reaction(reactants=[tihiolGroupContaining1, tihiolGroupContaining2], reaction_smart=reaction_smart)
    return products, reaction_smart
