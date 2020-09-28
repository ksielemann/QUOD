### Katharina Sielemann ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1 ###

#prior to this analysis: run QUOD.py for the (I) whole dataset including all accessions and (II) the replicate dataset of the same accession


#imports
import os, glob, sys
from argparse import ArgumentParser
import scipy.stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


__usage__ = """
python3 variance_in_repl_test.py
--all_covs <FULL_PATH_TO_accession_coverage.txt_FILE_FOR_ALL_ACCESSIONS> (in QUOD output folder of whole dataset)
--replicate_scores <FULL_PATH_TO_gene_dispensability_scores.csv_FILE_FOR_REPLICATE_DATASET> (in QUOD output folder of replicate dataset)
--output_dir <FULL_PATH_TO_OUPUT_FOLDER>
REQUIREMENTS: os, glob, sys, argparse, scipy, pandas, numpy, matplotlib
			"""

input_parameters = ArgumentParser()
input_parameters.add_argument("--all_covs", dest="all_covs")
input_parameters.add_argument("--replicate_scores", dest="repl_scores")
input_parameters.add_argument("--output_dir", "--output", "--out", "--o", dest="output_directory")

if "--help" in sys.argv or "-h" in sys.argv:
    print(__usage__)
    sys.exit(1)
    
args = input_parameters.parse_args()
if args.all_covs is None:
	print("\n'--all_covs' was not set'")
	print(__usage__)
	sys.exit(1)	
elif args.repl_scores is None:
	print("\n'--replicate_scores' was not set'")
	print(__usage__)
	sys.exit(1)
elif args.output_directory is None:
	print("\n'--output_dir' was not set'")
	print(__usage__)
	sys.exit(1)
else:
	#iterative (100x) random (n=14) selection of datasets
	input_matrix = pd.read_csv(args.all_covs, sep="\t")
	input_matrix = input_matrix.set_index("gene")
	if os.path.isdir(args.output_directory + "iterative_random_sets_scores/") == False:
		os.makedirs(args.output_directory + "iterative_random_sets_scores/")
	for i in range(100):
		print("iteration " + str(i+1))
		cov_matrix = input_matrix.sample(14, axis="columns").copy()
		N = len(list(cov_matrix))
		cX = cov_matrix.div(cov_matrix.mean(axis="index"), axis="columns")
		ds = 1/((cX.sum(axis="columns")).div(N))
		ds.to_csv(args.output_directory + "iterative_random_sets_scores/" + str(i+1) + "_ds.csv", header=False)


	datei = open(args.repl_scores,"r")
	lines = datei.readlines()
	datei.close()
	all_values = []
	for line in lines:
		line = line.strip().split(",")
		if float(line[1]) != np.inf:
			all_values.append(float(line[1]))
	max_number = max(all_values)

	#take average of iterative random sets
	score_files = glob.glob(args.output_directory + "iterative_random_sets_scores/*_ds.csv")
	variances = []
	stdeviations = []
	for scorefile in score_files:
		datei = open(scorefile,"r")
		ds = datei.readlines()
		datei.close()
		scores = []
		for line in sorted(list(ds)):
			line = line.strip().split(",")
			if float(line[1]) >= max_number:
				value = max_number
				scores.append(float(value))
			else:
				scores.append(float(line[1]))					
		variances.append(np.var(scores))
		stdeviations.append(np.std(scores))

	#Col-0 replicates (n=14)
	datei = open(args.repl_scores,"r")
	ds = datei.readlines()
	datei.close()	

	all_values = []
	for line in ds:
		line = line.strip().split(",")
		if line[1] != np.inf:
			all_values.append(line[1])
	max_number = max(all_values)

	scores_col = []
	for line in sorted(list(ds)):
		line = line.strip().split(",")
		if line[1] >= max_number:
			value = max_number
			scores_col.append(float(value))
		else:
			scores_col.append(float(line[1]))


	#boxplot
	my_dict = {'randomly,\niterative\nselected datasets': variances, 'replicates': np.var(scores_col)}
	medianprops = dict(linestyle='-', linewidth=1.5, color='royalblue')
	meanprops = dict(linestyle='dashed', linewidth=1.5, color='royalblue')
	plt.figure(figsize=(6,7))
	data = [variances, np.var(scores_col)]
	plt.ylim(ymin=0, ymax = max(variances))
	plt.xlim(xmin=0.8, xmax=2.2)

	plt.scatter(2, np.var(scores_col), color='grey', s=4, alpha=1)
	x = np.random.normal(1, 0.04, len(variances))
	plt.scatter(x, variances, color='grey', s=1, alpha=0.4)   
	plt.boxplot(my_dict['randomly,\niterative\nselected datasets'], showmeans=True, meanline=True, showfliers=False, medianprops=medianprops, meanprops=meanprops, widths=0.3)
	plt.xticks([1,2],["randomly,\niterative\nselected datasets", "replicates"])
	plt.ylabel('variance of the dispensability score (ds)', fontsize=16)
	plt.tick_params(axis='both', which='both', labelsize=14)
	plt.savefig(args.output_directory + "variance_in_replicates.png")
	plt.close()


	#calculate mean scores of iterative randomly selected datasets
	score_files = glob.glob(args.output_directory + "iterative_random_sets_scores/*_ds.csv")
	gene_scores = {}
	for scorefile in score_files:
		datei = open(scorefile,"r")
		ds = datei.readlines()
		datei.close()
		for zeile in ds:
			line = zeile.strip().split(",")
			#print(line[0])
			if line[0] not in gene_scores.keys():
				if float(line[1]) <= float(max_number):
					gene_scores[line[0]] = [float(line[1])]
				else:
					gene_scores[line[0]] = [float(max_number)]
			else:
				if float(line[1]) <= float(max_number):
					gene_scores[line[0]] += [float(line[1])]
				else:
					gene_scores[line[0]] += [float(max_number)]		


	scores_randomsets = []
	for gene in gene_scores.keys():
		scores_randomsets.append(np.mean(gene_scores[gene]))	


	#calculate differences
	print("\nreplicates (var,std):")
	print(np.var(scores_col))
	print(np.std(scores_col))
	print("\niterative randomly selected datasets (var,std):")
	print(np.var(scores_randomsets))
	print(np.std(scores_randomsets))

	print("\n")
	print(scipy.stats.levene(scores_randomsets, scores))

