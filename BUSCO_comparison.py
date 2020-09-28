### Katharina Sielemann ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1 ###

#/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/: 1.tsv = result of BUSCO run using reference 1 (e.g. embryophyta)
#/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/: 2.tsv = result of BUSCO run using reference 2 (e.g. chlorophyta)
#/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/: 3.tsv = result of BUSCO run using reference 3 (e.g. brasssicales)

#/PATH/TO/QUOD/OUTPUT/gene_dispensability_scores.csv

#imports
import plotly.graph_objs as go
from plotly.offline import plot
import scipy.stats, random
import numpy as np
import pandas as pd
import dabest

#read files
datei = open("/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/1.tsv", "r")
busco_results = datei.readlines()
datei.close()

datei = open("/PATH/TO/QUOD/OUTPUT/gene_dispensability_scores.csv", "r")
ds_results = datei.readlines()
datei.close()

#extract genes and scores (BUSCO vs. non-BUSCO)
all_genesandscores = {}
all_genes = []
all_scores = []
for line in ds_results:
	elements = line.strip().split(",")
	gene = elements[0].split(".")[1]
	all_genes.append(gene)
	all_scores.append(float(elements[1]))
	all_genesandscores[gene] = float(elements[1])
    
busco_results = busco_results[5:]
busco_genes = []
busco_scores = []
for line in busco_results:
	try:
	  line = line.strip().split("\t")
	  if line[1] == "Complete":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Duplicated":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Fragmented":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])  
	except KeyError:
		pass
		
all_genesandscores2 = {}
all_genes2 = []
all_scores2 = []
for index, gene in enumerate(all_genes):
	if gene not in busco_genes: #non-BUSCO genes
		all_genes2.append(gene)
		all_scores2.append(float(all_scores[index]))
		all_genesandscores2[gene] = float(all_scores[index])


max_number_1 = max([value for value in list(all_scores2) if value != np.inf])
max_number_2 = max([value for value in list(busco_scores) if value != np.inf])
max_number = max(max_number_1, max_number_2)
Enon_BUSCO_plot_data = []
for number in sorted(list(all_scores2)):
	if number >= max_number:
		value = max_number
		Enon_BUSCO_plot_data.append(value)
	else:
		Enon_BUSCO_plot_data.append(number)
EBUSCO_plot_data = []
for number in sorted(list(busco_scores)):
	if number >= max_number:
		value = max_number
		EBUSCO_plot_data.append(value)
	else:
		EBUSCO_plot_data.append(number)		


datei = open("/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/2.tsv", "r")
busco_results = datei.readlines()
datei.close()

datei = open("/PATH/TO/QUOD/OUTPUT/gene_dispensability_scores.csv", "r")
ds_results = datei.readlines()
datei.close()

#extract genes and scores (BUSCO vs. non-BUSCO)
all_genesandscores = {}
all_genes = []
all_scores = []
for line in ds_results:
	elements = line.strip().split(",")
	gene = elements[0].split(".")[1]
	all_genes.append(gene)
	all_scores.append(float(elements[1]))
	all_genesandscores[gene] = float(elements[1])
    
busco_results = busco_results[5:]
busco_genes = []
busco_scores = []
for line in busco_results:
	try:
	  line = line.strip().split("\t")
	  if line[1] == "Complete":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Duplicated":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Fragmented":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])  
	except KeyError:
		pass
		
all_genesandscores2 = {}
all_genes2 = []
all_scores2 = []
for index, gene in enumerate(all_genes):
	if gene not in busco_genes: #non-BUSCO genes
		all_genes2.append(gene)
		all_scores2.append(float(all_scores[index]))
		all_genesandscores2[gene] = float(all_scores[index])


max_number_1 = max([value for value in list(all_scores2) if value != np.inf])
max_number_2 = max([value for value in list(busco_scores) if value != np.inf])
max_number = max(max_number_1, max_number_2)
Cnon_BUSCO_plot_data = []
for number in sorted(list(all_scores2)):
	if number >= max_number:
		value = max_number
		Cnon_BUSCO_plot_data.append(value)
	else:
		Cnon_BUSCO_plot_data.append(number)
CBUSCO_plot_data = []
for number in sorted(list(busco_scores)):
	if number >= max_number:
		value = max_number
		CBUSCO_plot_data.append(value)
	else:
		CBUSCO_plot_data.append(number)	


datei = open("/PATH/TO/BUSCO/OUTPUT/FULL/TABLE/3.tsv", "r")
busco_results = datei.readlines()
datei.close()

datei = open("/PATH/TO/QUOD/OUTPUT/gene_dispensability_scores.csv", "r")
ds_results = datei.readlines()
datei.close()

#extract genes and scores (BUSCO vs. non-BUSCO)
all_genesandscores = {}
all_genes = []
all_scores = []
for line in ds_results:
	elements = line.strip().split(",")
	gene = elements[0].split(".")[1]
	all_genes.append(gene)
	all_scores.append(float(elements[1]))
	all_genesandscores[gene] = float(elements[1])
    
busco_results = busco_results[5:]
busco_genes = []
busco_scores = []
for line in busco_results:
	try:
	  line = line.strip().split("\t")
	  if line[1] == "Complete":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Duplicated":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])
	  elif line[1] == "Fragmented":
		  busco_genes.append(line[2])
		  busco_scores.append(all_genesandscores[line[2]])  
	except KeyError:
		pass
		
all_genesandscores2 = {}
all_genes2 = []
all_scores2 = []
for index, gene in enumerate(all_genes):
	if gene not in busco_genes: #non-BUSCO genes
		all_genes2.append(gene)
		all_scores2.append(float(all_scores[index]))
		all_genesandscores2[gene] = float(all_scores[index])


max_number_1 = max([value for value in list(all_scores2) if value != np.inf])
max_number_2 = max([value for value in list(busco_scores) if value != np.inf])
max_number = max(max_number_1, max_number_2)
Bnon_BUSCO_plot_data = []
for number in sorted(list(all_scores2)):
	if number >= max_number:
		value = max_number
		Bnon_BUSCO_plot_data.append(value)
	else:
		Bnon_BUSCO_plot_data.append(number)
BBUSCO_plot_data = []
for number in sorted(list(busco_scores)):
	if number >= max_number:
		value = max_number
		BBUSCO_plot_data.append(value)
	else:
		BBUSCO_plot_data.append(number)	
		

#Levene's test to test for equal variances: Is the variance larger for one distribution (e.g. non-BUSCOs)?
#output: test statistic, p-value
import numpy as np
print("BUSCOs, non-BUSCOs (variance)")
print("1")
print(np.var(EBUSCO_plot_data))
print(np.var(Enon_BUSCO_plot_data))
print(scipy.stats.levene(EBUSCO_plot_data,Enon_BUSCO_plot_data))

print("2")
print(np.var(CBUSCO_plot_data))
print(np.var(Cnon_BUSCO_plot_data))
print(scipy.stats.levene(CBUSCO_plot_data,Cnon_BUSCO_plot_data))

print("3")
print(np.var(BBUSCO_plot_data))
print(np.var(Bnon_BUSCO_plot_data))
print(scipy.stats.levene(BBUSCO_plot_data,Bnon_BUSCO_plot_data))

		
#dabest
dict_data = {"BUSCO (2)":pd.Series(CBUSCO_plot_data), "non-BUSCO (2)":pd.Series(Cnon_BUSCO_plot_data),
			"BUSCO (1)":pd.Series(EBUSCO_plot_data), "non-BUSCO (1)":pd.Series(Enon_BUSCO_plot_data),
			"BUSCO (3)":pd.Series(BBUSCO_plot_data), "non-BUSCO (3)":pd.Series(Bnon_BUSCO_plot_data)}
df = pd.DataFrame(dict_data)


multi = dabest.load(df,idx=(("BUSCO (2)","non-BUSCO (2)"),
							("BUSCO (1)","non-BUSCO (1)"),
							("BUSCO (3)","non-BUSCO (3)")))

(multi.mean_diff.statistical_tests).to_csv("/FULL/PATH/TO/OUTPUT/DIRECTORY/dabest_BUSCO_comparison.csv")




