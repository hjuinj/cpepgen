
#############################################################################
# Example Usage Below
#############################################################################

# omphalotinA = CycloPeptide("WVIVVGVIGVIG", [2,4,5,6,7,8,9,11,12]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = CycloPeptide("wVIvVGVIGVIG", [2,4,5,6,7,8,9,11,12]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = CycloPeptide("AAAAAAAAAAAA") #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = CycloPeptide("LALGPLALGP", [2,4,7,9]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = CycloPeptide("LALpTLTLpT", [2,5,7,10]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = Peptide("WWW")
# omphalotinA = CycloPeptide("LHLpALHLpA", [2,5,7,10]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# omphalotinA = CycloPeptide("LALpALTLpA", [2,5,7,10]) #Trp - meVal - Ile - meVal - meVal - meGly - meVal - meIle - meGly - Val - Ile - Gly  9 methlyations
# print(Chem.MolToPDBBlock(omphalotinA.peptide))
#
# print(Chem.MolToSmiles(omphalotinA.peptide, isomericSmiles = True, allBondsExplicit = True))

# omphalotinA.peptide= omphalotinA.gromosff_hydrogenate()
# print(Chem.MolToPDBBlock(omphalotinA.peptide))
#
# omphalotinA.gen_conf(1)
# print(Chem.MolToPDBBlock(omphalotinA.peptide))
# Chem.MolToPDBFile(omphalotinA.peptide, "init.pdb")
# for i in range(3):
#     Chem.MolToPDBFile( omphalotinA.peptide, "{}{}".format("./tmp", i), confId = i)

"""
MDAnalysis.core.AtomGroup.AtomGroup(sorted(x.atoms, key = lambda x : x.resid))
"""

# tmp = CycloPeptide("MEA", methylation_sites = [1,2,3])
# print(Chem.MolToPDBBlock(omphalotinA.peptide))
#((omphalotinA.peptide))
