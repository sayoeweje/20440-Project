 # Overview

This repository contains the necessary scripts and data to conduct a structure-based network analysis of Streptococcus pyogenes Cas9 (SpCas9) when bound to a single guide RNA (sgRNA) and its target DNA. The structure is derived from a pdb file (4oo8.pdb) that contains a crystal structure reported in Nishimasu et al (2014).

The analysis is inspired by the work of Gaiha et al (2019), who first conceived this method. In this network-based protein representation, amino acid residues are defined as nodes, and sums of the energies of the non-covalent interactions between residues are defined as weighted edges. To construct the residue interaction network, both an energetic network and a centroid network were created. In the energetic network, the bond energies of all non-covalent interactions between residues (van der Waals interactions, hydrogen bonds, water-bridged bonds, salt bridges, disulfide bonds, pi-pi interactions, pi-cation interactions and metal coordinated bonds) were calculated and combined into a summed weighted edge. Only interactions between terminal side chain atoms were considered. The centroid network, originally created to account for the contributions of hydrophobic packing to protein folding, contained edges present only if the distance between two residues was ≤ 8.5 Å (excluding immediately neighboring residues in the protein chain). The edges from the energetic and centroid networks were summed, and the following centrality metrics were calculated: 

- Second Order Degree Centrality (SD): The number of second-order interactions between a residue and residues in other higher order structures (i.e. secondary structures, averages results from analysis using classically defined secondary structures using Stride and inferred secondary structures using a random-walk approach). 
- Summed Edge Betweenness (SEB): The frequency that a node’s edges are used as the shortest path between all nodes in the network.
- Residue Ligand Proximity (RP): The euclidean distance of a residue’s center of mass to the center of mass of the protein’s ligand.

where the final network score = SD + SEB - RP (each value individually normalized). A residue's network score is indicative of its relative importance to the protein network. To determine the extent to which the assigned network scores are predictive of the true importance of individual residues, a linear correlation analysis is conducted between this data and in vitro activity data of a library of ~8500 Cas9 mutants (derived from Spencer et Zhang 2017). 

Citations:

Nishimasu, H., Ran, F. A., Hsu, P. D., Konermann, S., Shehata, S. I., Dohmae, N., Ishitani, 
R., Zhang, F., & Nureki, O. (2014). Crystal structure of Cas9 in complex with guide RNA and 
target DNA. Cell, 156(5), 935–949.

Gaiha, G. D., Rossin, E. J., Urbach, J., Landeros, C., Collins, D. R., Nwonu, C., ... & Walker, B. D. (2019). Structural topology defines protective CD8+ T cell epitopes in the HIV proteome. Science, 364(6439), 480-484.

Spencer, J. M., & Zhang, X. (2017). Deep mutational scanning of S. pyogenes Cas9 reveals important functional domains. Scientific reports, 7(1), 1-14.






# src

## Network Analysis Pipeline

### I. Initialize directory structure.

	1. Create a directory in your chosen location called "NetworkAnalysis".

	> mkdir NetworkAnalysis
	> cd NetworkAnalysis

	2. Populate this directory with the dependencies provided in this Git repo. You will need:
		- 4oo8.pdb
		- AAarea
		- allatoms_radii
		- atomicWeights
		- buildEdges.sh
		- calculateFinalScore.R
		- makeWaterBonds.py
		- modularity_analysis_centroid.R
		- modularity_analysis_energetics.R
		- pdb2edgeSC.py
		- pdb2polypeptide.py
		- terminalAtoms
		- uniqueAtoms

	3. For each structure you'd like to analyze, create a new subdirectory named for the PDB code, and include the relevant PDB file within that directory (here we've provided 4oo8 as an example)). 

	> mkdir 4oo8
	> mv 4oo8.pdb 4oo8
	

### II. Build the network.

	1. Run the buildEdges.sh script from within the NetworkAnalysis folder.

	> sh buildEdges.sh 4oo8/4oo8.pdb

	2. You will now see a number of files within the 1BTL/ directory, which correspond to the intermediate and final network files. The final file to be analyzed is 1BTL_net.

### III. Analyze the network.

There are 3 steps to analyzing the network: the energetics analysis, the centroid analysis, and the final analysis. Note that most of the analysis (energetics and centroid) is already conducted by buildEdges.sh. The final analysis needs to be run separately using the calculateFinalScore.R script. Run this from within the 1BTL/ directory, with the following usage:

> cd 1BTL
> R --vanilla --args 1BTL/1BTL/ < ../calculateFinalScore.R 

The output is named "FinalSum" and is a list of nodes (amino acids) with their scores as follows:

head FinalSum

GLY156 -1.12139304234615
HIS158 -1.64504045031281
ALA150 -1.9452762424746
GLU147 -1.52075794536626
LEU139 1.15367020700574
LEU162 2.19848342502133
ALA280 -0.961502115446914
ILE279 -0.00840942137073508
ARG164 1.77030331167432
GLU171 -0.735560986652963
These are the final network scores.

## Network Visualization

network_visualization.py can be used to create a plot of a network. Prior to running, move the pdb file and a file containing the network scores for each residue (i.e. the FinalSum file generated following Network Analysis) into the src folder. The output will be a pdf file.

## Linear Regression Analysis

mutational_network_analysis.py can be used to conduct a linear regression analysis between the residue network scores and in vitro mutability data derived from the Spencer et Zhang data set. Prior to running, move AminoAcids_Domains_Mutability.xlsx, AminoAcids_NetworkScores.xlsx, Mutability_NetworkScores.xlsx, Spencer_et_al_2017_Cas9_mutagenesis.xlsx, and the appropriate FinalSum file into the src folder. 

# Data

Excel files: Contains sheets with network scores, mutability values, and both for use in linear regression analysis. 

Network scores: FinalSum files contain the network scores for residues following specific analyses of target bound Cas9 RNP. Centrality metrics used during analysis are specified by the file name.

PDB files: 4oo8 (pdb file for target bound Cas9 RNP)

# Folder Structure
data: Contains needed excel files, network score files (generated in Network Analysis), and pdb file.
figures: Contains linear regression plots, network plots.
NetworkAnalysis: Contains scripts from Gaiha et al necessary to conduct network analysis, along with raw data following analysis of target DNA bound Cas9 RNP. 
src: Contains scripts for linear regression analysis and network visualization.
stride: software used in network analysis to assign secondary structures. 

# Installation

## Dependencies for Network Analysis 

	1. Download PHENIX (Python-based Hierarchical ENvironment for Integrated Xtallography) for protonation.
		- Download binaries from https://www.phenix-online.org/download/ (after requesting a password via email).
		- Follow installation instructions and note the final location; if you are using a Mac, it will likely be installed in your Applications directory. 
		- Make sure that line 18 in the file buildEdges.sh reflects this path. 

	2. Download STRIDE for secondary structure assignment.
		- Download binaries from http://webclu.bio.wzw.tum.de/stride/install.html
		- Create a directory called "Stride" one level above your NetworkAnalysis folder and copy any Stride-related files into this directory.

## Dependencies for Network Visualization and Linear Correlation Analysis

To run the code, the user needs to have an IDE that can run python scripts. The
user also must also have 4008.pdb in the same folder as the python script.

The packages necessary are:
Bio version 1.3.3
ipython version 8.1.1
matplotlib version 3.5.1
networkx version 2.6.3
numpy version 1.21.2