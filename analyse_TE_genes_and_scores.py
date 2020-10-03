### Katharina Sielemann ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1 ###

#/FULL/PATH/TO/LIST/OF/TE/GENES.txt (contains one gene per line)

#imports
import plotly.graph_objs as go
from plotly.offline import plot
import scipy.stats, random
import numpy as np

#read files
datei = open("/FULL/PATH/TO/LIST/OF/TE/GENES.txt", "r")
te_genes = datei.readlines()
datei.close()

datei = open("/FULL/PATH/TO/QUOD/OUTPUT/gene_dispensability_scores.csv", "r")
ds_results = datei.readlines()
datei.close()

te_genes = [value.strip() for value in te_genes]
#extract genes and scores (TE vs. non-TE)
TE_scores = []
non_TE_scores = []
for line in ds_results:
	line = line.strip().split(",")
	gene = line[0].split(".")[1]
	if gene in te_genes:
		TE_scores.append(float(line[1]))
	else:
		non_TE_scores.append(float(line[1]))


max_number_1 = max([value for value in list(TE_scores) if value != np.inf])
max_number_2 = max([value for value in list(non_TE_scores) if value != np.inf])
max_number = max(max_number_1, max_number_2)

filtered_non_TE_scores = []
for number in non_TE_scores:
	if number >= max_number:
		value = max_number
		filtered_non_TE_scores.append(value)
	else:
		filtered_non_TE_scores.append(number)

filtered_TE_scores = []
for number in TE_scores:
	if number >= max_number:
		value = max_number
		filtered_TE_scores.append(value)
	else:
		filtered_TE_scores.append(number)

#test for significance: Welch's t-test
print("mean of 'TE' scores: " + str(np.mean(filtered_TE_scores)) + "n=" + str(len(filtered_TE_scores)))
print("mean of 'non-TE' scores: " + str(np.mean(filtered_non_TE_scores)) + "n=" + str(len(filtered_non_TE_scores)))
print(scipy.stats.mannwhitneyu(filtered_TE_scores, filtered_non_TE_scores))


