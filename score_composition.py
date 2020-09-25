### Katharina Sielemann ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1 ###

#imports
import os, glob, sys
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt

__usage__ = """
python3 score_composition.py
--gene_file <FULL_PATH_TO_GENE_FILE> (one gene per line (gene name from gff file used for QUOD))
--score_file <FULL_PATH_TO_gene_dispensability_scores_FILE> (QUOD output file 'gene_dispensability_scores.csv')
--accession_cov_file <FULL_PATH_TO_accession_coverage_FILE> (QUOD output file 'accession_coverage.txt')
--output_dir <FULL_PATH_TO_OUPUT_FOLDER>
--visualize (optional argument)
REQUIREMENTS: os, glob, sys, argparse, numpy (optional: matplotlib for visualization)
			"""

input_parameters = ArgumentParser()
input_parameters.add_argument("--gene_file", dest="gene_file")
input_parameters.add_argument("--score_file", dest="score_file")
input_parameters.add_argument("--accession_cov_file", dest="cov_matrix")
input_parameters.add_argument("--output_dir", "--output", "--out", "--o", dest="output_directory")
input_parameters.add_argument("--visualize", action='store_true')

if "--help" in sys.argv or "-h" in sys.argv:
    print(__usage__)
    sys.exit(1)
    
args = input_parameters.parse_args()
if args.gene_file is None:
	print("\n'--gene_file' was not set'")
	print(__usage__)
	sys.exit(1)
elif args.score_file is None:
	print("\n'--score_file' was not set'")
	print(__usage__)
	sys.exit(1)
elif args.cov_matrix is None:
	print("\n'--accession_cov_file' was not set'")
	print(__usage__)
	sys.exit(1)
elif args.output_directory is None:
	print("\n'--output_dir' was not set'")
	print(__usage__)
	sys.exit(1)	
else:
	#extracting read coverage of the selected genes
	print("\nextracting read coverage of the selected genes...")
	datei = open(args.gene_file, "r")
	gene_lines = datei.readlines()
	datei.close()
	genes = [x.strip().rstrip() for x in gene_lines]
	#coverages per gene
	datei = open(args.cov_matrix, "r")
	cov_lines = datei.readlines()
	datei.close()
	prefix, accessions = cov_lines[0].split("\t",1)
	dict_gene_covs = {}
	cov_lines = cov_lines[1:]
	for line in cov_lines:
		covs = line.strip().split("\t")[1:]
		gene = line.strip().split("\t")[0]
		dict_gene_covs[gene] = covs
	#score per gene
	datei = open(args.score_file, "r")
	score_lines = datei.readlines()
	datei.close()
	dict_gene_score = {}
	for line in score_lines:
		line = line.strip().split(",")
		dict_gene_score[line[0]] = line[1]
	#table	
	out = ["gene_name\tds_score\ttotal_mean_cov\tmean_cov_of_lowest_10%\tmean_cov_of_highest_10%\tnumber_of_zero_cov_accessions\t" + accessions]
	for gene in genes:
		covs_floats = [float(x) for x in dict_gene_covs[gene]]
		covs_floats = sorted(covs_floats)
		lowest_10_list = covs_floats[:int(len(covs_floats) * .1)]
		highest_10_list = covs_floats[int(len(covs_floats) * .9):]
		total_mean = sum(covs_floats)/len(covs_floats)
		mean_lowest = sum(lowest_10_list)/len(lowest_10_list)
		mean_highest = sum(highest_10_list)/len(highest_10_list)
		out.append(gene + "\t" + dict_gene_score[gene] + "\t" + str(total_mean) + "\t" + str(mean_lowest) + "\t" + str(mean_highest) + "\t" + str(len([x for x in covs_floats if x==0])) + "\t" + "\t".join(dict_gene_covs[gene]) + "\n")
	datei = open(args.output_directory + "score_composition.tsv","w")
	datei.write(("").join(out))
	datei.close()

if args.visualize == True:
	print("\nvisualizing...")
	#plot distribution of the coverage in all samples; one box plot per gene
	if os.path.isdir(args.output_directory + "plots/") == False:
		os.makedirs(args.output_directory + "plots/")
	for gene in genes:
		print(gene)
		my_dict = {'coverage': [float(x) for x in dict_gene_covs[gene]]}
		medianprops = dict(linestyle='-', linewidth=1.5, color='royalblue')
		meanprops = dict(linestyle='dashed', linewidth=1.5, color='royalblue')
		plt.figure(figsize=(6,7))
		data = [[float(x) for x in dict_gene_covs[gene]]]
		plt.ylim(ymin=0, ymax = max([float(x) for x in dict_gene_covs[gene]]))
		plt.xlim(xmin=0.9, xmax=1.1)
		y = [float(x) for x in dict_gene_covs[gene]]
		x = np.random.normal(1, 0.0125, len(y))
		plt.scatter(x, y, color='grey', s=1, alpha=0.4)   
		plt.boxplot(my_dict.values(), showmeans=True, meanline=True, showfliers=False, medianprops=medianprops, meanprops=meanprops, widths=0.125, vert=True)
		plt.xticks([1],[gene])
		plt.ylabel('coverage', fontsize=16)
		plt.tick_params(axis='both', which='both', labelsize=14)
		plt.plot([np.mean([float(x) for x in dict_gene_covs[gene]]),np.mean([float(x) for x in dict_gene_covs[gene]])],[0.9,1.1], '--', linewidth=1, color='royalblue', alpha=0.5)
		plt.savefig(args.output_directory + "plots/" + gene + "_coverage_distribution_boxplot.png")
		plt.close()
	
print("\n---finished---")
