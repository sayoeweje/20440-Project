import numpy as np
import pandas as pd
from statistics import mean
from statistics import pstdev
from scipy import stats
import csv
import re
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

mutability_scores = pd.read_excel(io='Spencer_et_al_2017_Cas9_mutagenesis.xlsx',
                   sheet_name='Mutability Scores')

mutability_scores_list = mutability_scores['Mutability Score'].tolist()
mean_mutability_score = mean(mutability_scores_list)
stdev_mutability_score = pstdev(mutability_scores_list)

for index, row in mutability_scores.iterrows():
    normalized_mutability_score = (row[1] - mean_mutability_score)/stdev_mutability_score
    mutability_scores = mutability_scores.replace(row[1], normalized_mutability_score)

df1 = pd.read_excel(io='Spencer_et_al_2017_Cas9_mutagenesis.xlsx',
                   sheet_name='All Count Data')

df1.drop_duplicates(subset=['AA Position'], inplace=True)

s1 = pd.Series(df1['AA Position'], name='AA Position')
s2 = pd.Series(df1['Domain'], name='Domain')

aa_domains = pd.concat([s1,s2], axis=1)

aa_domains_mutability = pd.merge(aa_domains, mutability_scores, on=['AA Position'])

aa_domains_mutability = aa_domains_mutability.dropna()

aa_domains_mutability.to_excel('AminoAcids_Domains_Mutability.xlsx')

with open('FinalSum') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    aa_network_scores = pd.DataFrame(csv_reader, columns=['AA Position','Network Score'])

aa_network_scores = aa_network_scores[aa_network_scores['Network Score'] != 'NA']

aa_network_scores['Network Score'] = aa_network_scores['Network Score'].astype(float)
aa_network_scores = aa_network_scores.sort_values(by=['Network Score'])

for aa in aa_network_scores['AA Position']:
    aa_new = int(re.sub('\D', '', aa))
    aa_network_scores = aa_network_scores.replace(to_replace=aa, value=aa_new)


aa_network_scores = aa_network_scores.sort_values(by=['AA Position'])

aa_network_scores.to_excel('AminoAcids_NetworkScores.xlsx')

aa_domains_network_mutability = pd.merge(aa_network_scores, aa_domains_mutability, on=['AA Position'])
columns_swap = ['AA Position','Domain','Network Score','Mutability Score']
aa_domains_network_mutability=aa_domains_network_mutability.reindex(columns=columns_swap)

aa_domains_network_mutability.to_excel('Mutability_NetworkScores.xlsx')

fig = plt.figure(figsize = (15,10))
x_vals = aa_domains_network_mutability.loc[:,'Network Score']
y_vals = aa_domains_network_mutability.loc[:,'Mutability Score']
plt.scatter(x_vals, y_vals, s=10, color='dodgerblue')
m, b = np.polyfit(x_vals, y_vals, 1)
plt.plot(x_vals, m*x_vals + b, color='orange')
plt.xlabel('Network Score', fontsize = 15)
plt.ylabel('Mutational Tolerance', fontsize = 15)
plt.title('Structural Network Analysis of Cas9 Mutational Tolerance', fontsize = 20)
model_x = np.array(x_vals).reshape((-1,1))
model_y = np.array(y_vals)
correlation, pvalue = stats.spearmanr(model_x, model_y)
fig_text = 'Correlation = {} \n p-value = {}'.format(correlation,pvalue)
plt.text(6.5,5, fig_text, fontsize = 15)
plt.show()