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
structure_id = "Cas"
filename = "4oo8_Cas.pdb" # PDB file containing Cas9 + sgRNA, target DNA crystal structure
structure = parser.get_structure(structure_id, filename)

model = structure[0] # Get structural elements of pdb
protein = model["A"] # Selects Chain A in pdb file (Cas 9 protein)

NetworkScores = {}
with open("FinalSum_Cas") as file:
    for line in file:
 
        (key, value) = line.split()

        NetworkScores[key] = float(value)

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
        networkscore_residue = NetworkScores.get(residue.get_resname()+str(i+3))
        residue_name = residue.get_resname()+str(i+3)
        graph.add_node(residue_name, pos=(ca[0], ca[1]), node_color = networkscore_residue, 
        node_size = 10*(1.271323607873502+networkscore_residue))
    except:
        # Exception for any residues that are missing from structure
        print("Cas9 Residue "+ str(residue) + ": This residue appears to be missing")

# Plots network
fig = figure(num=1, figsize=(10, 15), dpi=300)
ax  = fig.add_subplot(111)
pos = nx.get_node_attributes(graph,'pos')

node_colors_dic = nx.get_node_attributes(graph, 'node_color')
node_size_dic = nx.get_node_attributes(graph, 'node_size')
node_names = list(node_colors_dic.keys())

for i in range(len(node_names)):
    try:
        graph.add_edge(node_names[i], node_names[i+1])
    except:
        print("End of node list")

node_colors_array = list(node_colors_dic.values())
node_size_array = list(node_size_dic.values())

g = nx.draw_networkx(graph, pos, node_size = node_size_array,
                     node_color = node_colors_array, with_labels = False, 
                                    cmap = plt.cm.coolwarm)

norm = mpl.colors.Normalize(vmin=min(list(NetworkScores.values())),vmax=max(list(NetworkScores.values())))

cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=plt.cm.coolwarm), ax=ax, 
             orientation='horizontal', location='bottom', pad=0)

cbar.ax.tick_params(labelsize=20)
cbar.set_label(label='Network Score', size = 20)
ax.set_axis_off()
ax.set_title("2D Representation of Cas9 Structural Network")
ax.margins(x=0, y=0)

plt.savefig("4oo8_network_"+structure_id+".pdf", bbox_inches = 'tight',
    pad_inches = 0.5)