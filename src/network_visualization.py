import scipy as sp
from pylab import *
import numpy as np
import matplotlib as mpl
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
import matplotlib.patches as mpatches

parser = PDBParser() # Initiates PDBParser object to read .pdb file
structure_id = "Cas9"
filename = "4oo8.pdb" # PDB file containing Cas9 + sgRNA, target DNA crystal structure
structure = parser.get_structure(structure_id, filename)

model = structure[0] # Get structural elements of pdb
protein = model["A"] # Selects Chain A in pdb file (Cas 9 protein)
guide = model["B"] # Selects Chain B in pdb file (sgRNA)
target_DNA = model["C"] # Selects Chain C in pdb file (Target DNA)

graph = nx.Graph() 
locations = [] # Array to hold 3D positions of each residue in structure
for i in range(len(protein)):
    try:
        residue = protein[i+3] # Cas9 starts with a Lys3 (Res 1 and 2 undefined, 1-indexed)
        
        # Gets 3D location of residue based on alpha carbon position
        ca = residue["CA"].get_vector()

        residue_info = [residue.get_resname()+str(i+3), ca[0], ca[1], ca[2]] 
        locations.append(residue_info) # Add residue to location array

        # Creates node for residue in 2-D network, sets color
        graph.add_node(residue.get_resname()+str(i+3), pos=(ca[0], ca[1]), color='gray')
    except:
        # Exception for any residues that are missing from structure
        print("Cas9 Residue "+ str(residue) + ": This residue appears to be missing")
        
for j in range(len(guide)):
    try:
        nucleotide = guide[j+1] 

        # Gets 3D location of base based on phosphate position
        p = nucleotide["P"].get_vector() 

        nucleotide_info = [nucleotide.get_resname()+str(j+1), p[0], p[1], p[2]]
        locations.append(nucleotide_info) # Add base to location array

        # Creates node for base in 2-D network, sets color
        graph.add_node(nucleotide.get_resname()+str(j+1), pos=(p[0], p[1]), color='orange')
    except:
        # Exception for any bases that are missing from structure
        print("sgRNA Base "+ str(nucleotide) + ": This base appears to be missing")

for k in range(len(target_DNA)):
    try:
        nucleotide = target_DNA[k+1]

         # Gets 3D location of base based on phosphate position
        p = nucleotide["P"].get_vector()

        nucleotide_info = [nucleotide.get_resname()+str(k+1), p[0], p[1], p[2]]
        locations.append(nucleotide_info) # Add base to location array

        #  Creates node for base in 2-D network, sets color
        graph.add_node(nucleotide.get_resname()+str(k+1), pos=(p[0], p[1]), color='green')
    except:
        # Exception for any bases that are missing from structure
        print("Target DNA Base Pair "+ str(nucleotide) + ": This base appears to be missing")

# Calculates distances between all residues
for x in range(len(locations)):
    residue = locations[x]
    for y in range(len(locations)):
        residue_test = locations[y]
        distance = np.sqrt((float(residue_test[1])-float(residue[1]))**2
                           + (float(residue_test[2])-float(residue[2]))**2 
                           + (float(residue_test[3])-float(residue[3]))**2) 
        if distance < 8.5 and abs(x-y) > 1: # 8.5 Angstorms is threshold distance for 
            # contribution to hydrophobic packing
            graph.add_edge(residue[0], residue_test[0], weight=distance)

# Plots network

fig = figure(num=1, figsize=(10, 15), dpi=300)
ax  = fig.add_subplot(111)
pos = nx.get_node_attributes(graph,'pos')
color = nx.get_node_attributes(graph, 'color')
cmap = plt.cm.coolwarm
edges, weights = zip(*nx.get_edge_attributes(graph, 'weight').items())
g = nx.draw_networkx(graph, pos, node_size = 20, with_labels = False, 
                     node_color = color.values(), edge_color=weights, edge_cmap=cmap)

# Assigns color/legend values to nodes
handles, labels = ax.get_legend_handles_labels() 
handles.append(mpatches.Patch(color='grey',label='Cas9 residues'))
handles.append(mpatches.Patch(color='orange',label='sgRNA'))
handles.append(mpatches.Patch(color='green',label='Target DNA strand'))
ax.set_axis_off()
ax.set_title("2D Representation of Cas9 Structural Network")
ax.legend(handles=handles, loc='best')
ax.margins(x=0, y=0)

# Normalizes weights of edges in range of distances under packing threshold
norm = mpl.colors.Normalize(vmin=0,vmax=8.5)

fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
             orientation='horizontal', location='bottom', 
             label='Distance between Cas9 residues (Angstroms)', pad=0)

plt.savefig("4oo8_network.pdf", bbox_inches = 'tight',
    pad_inches = 0.5)





