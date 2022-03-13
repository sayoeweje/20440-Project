#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scipy as sp
from pylab import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import networkx as nx
import seaborn as sns
import pandas as pd  
import collections
from IPython.display import clear_output
import random
from Bio.PDB.PDBParser import PDBParser
from mpl_toolkits.mplot3d import Axes3D


# In[ ]:


parser = PDBParser()
structure_id = "Cas9"
filename = "4oo8.pdb"
structure = parser.get_structure(structure_id, filename)

for key, value in structure.header.items():
    print(key, ' : ', value)


# In[ ]:


model = structure[0]
protein = model["A"]
guide = model["B"]
target_DNA = model["C"]


# In[ ]:


Cas9_graph = nx.Graph()
for i in range(len(protein)):
    try:
        residue = protein[i+3]
        ca = residue["CA"].get_vector()
        Cas9_graph.add_node(residue.get_resname()+str(i), pos=(ca[0], ca[1]))
    except:
        print("Cas9 Residue "+ str(residue) + ": This residue appears to be missing")

fig = figure(num=1, figsize=(10, 10), dpi=80)
ax  = fig.add_subplot(111)
pos = nx.get_node_attributes(Cas9_graph,'pos')
g = nx.draw_networkx(Cas9_graph, pos, node_size = 5, with_labels = True, label_size = 2)
ax.set_axis_off()


# In[ ]:


for i in range(len(protein)):
    try:
        residue = protein[i+3]
        print("Cas9 Residue "+ str(i+3) + ": " + str(protein[i+3]))
        ca = residue["CA"].get_vector()
        print("\tAlpha carbon location: " + str(ca))
        print(ca[0])
        print(residue.get_resname())
    except:
        print("Cas9 Residue "+ str(residue) + ": This residue appears to be missing")


# In[ ]:


for j in range(len(guide)):
    try:
        nucleotide = guide[j+1]
        print("sgRNA Base "+ str(j+1) + ": " + str(guide[j+1]))
        print("\tPhosphate Location: "+nucleotide["P"].get_vector())
    except:
        print("sgRNA Base "+ str(nucleotide) + ": This residue appears to be missing")


# In[ ]:


for nucleotide in range(len(target_DNA)):
    try:
        print("Target DNA Base Pair "+ str(nucleotide) + ": " + str(target_DNA[nucleotide]))
    except:
        print("Target DNA Base Pair "+ str(nucleotide) + ": This residue appears to be missing")

