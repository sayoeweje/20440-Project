import posixpath
import sys
import string
import re

file = open(sys.argv[1],'r').readlines()
out = open(sys.argv[1].replace(".pdb","")+"_polypeptidePairs",'w')
A=""
B=""
Achain = ""
Bchain = ""

def isResidue(residue):
    if residue in ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","PRO","GLY",
    "ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP", "G", "C", "A", "U", "DG", "DC", "DA", "DT"]:
        return(True)
    else:
        return(False)

for line in file:
    if line[0:4]=="ATOM" and A!="":
        residuename = line[17:20];residuename=residuename.lstrip().strip()
        chainname = line[21];chainname=chainname.lstrip().strip()
        residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
        if A==residuename+residueposition:
            continue
        if isResidue(residuename) and residueposition.isdigit():
            B = residuename+residueposition
            Bchain = chainname
            posA_index = re.search(r"\d", A)
            posB_index = re.search(r"\d", B)
            posA = int(A[posA_index.start():])
            posB = int(B[posB_index.start():])
            if Achain==Bchain and posB == posA + 1:
                out.write(A+Achain+"\t"+B+Bchain+"\tSCORE\n")
    if line[0:4]=="ATOM":
        residuename = line[17:20];residuename=residuename.lstrip().strip()
        chainname = line[21];chainname=chainname.lstrip().strip()
        residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
        if isResidue(residuename) and residueposition.isdigit():
            A = residuename+residueposition
            Achain = chainname

out.close()
        
