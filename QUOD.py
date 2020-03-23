### Katharina Frey ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1.0 ###

#imports
import os, glob, sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

__usage__ = """
python3 calculate_dispensability_score.py
--input_dir <FULL_PATH_TO_FOLDER_INPUT_BAM_FILES> (file names = sample names)
--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILES> (optional argument)
--gff <FULL_PATH_TO_REFERENCE_ANNOTATION_FILE>
--output_dir <FULL_PATH_TO_OUPUT_FOLDER>
--min_cov_per_genome <INTEGER> (default = 10, optional argument)
--visualize (optional argument)
REQUIREMENTS: os, glob, sys, argparse, pandas, numpy (optional: matplotlib for visualization)
			"""

input_parameters = ArgumentParser()
input_parameters.add_argument("--input_dir", "--input", "--in", "--i", dest="input_directory", required=True)
input_parameters.add_argument("--bam_is_sorted", action='store_true')
input_parameters.add_argument("--gff", dest="gff_file", required=True)
input_parameters.add_argument("--output_dir", "--output", "--out", "--o", dest="output_directory", required=True)
input_parameters.add_argument("--min_cov_per_genome", dest="min_cov", default=10)
input_parameters.add_argument("--visualize", action='store_true')

if "--help" in sys.argv or "-h" in sys.argv:
    print(__usage__)
    sys.exit(1)
    
organelles = ["chrc", "chrm", "chloroplast", "mitochondria"]

args = input_parameters.parse_args()
if args.input_directory is None:
	print("\n'--input_dir' was not set'")
	print(__usage__)	
elif args.output_directory is None:
	print("\n'--output_dir' was not set'")
	print(__usage__)
elif args.gff_file is None:
	print("\n'--gff' was not set'")
	print(__usage__)	
else:
	#calculating read coverage depth per position
	print("\ncalculating read coverage depth per position...")
	if os.path.isdir(args.output_directory + "cov_files/") == False:
		os.makedirs(args.output_directory + "cov_files/")
	bam_files = glob.glob(args.input_directory + "*.bam*")
	for bamfile in bam_files:
		if bamfile.endswith(".bam.gz"):
			name = bamfile[(len(args.input_directory)):-7]
		else:	
			name = bamfile[(len(args.input_directory)):-4]
		samtools = "samtools"
		bedtools = "genomeCoverageBed"
		if os.path.isfile(args.output_directory + "cov_files/" + name + ".cov") == False:
			if args.bam_is_sorted == True:
				sorted_bam_file = bamfile
			else:
				sorted_bam_file =  bamfile + "_sorted.bam"
				cmd = samtools + " sort -m 5000000000 --threads 8 " + bamfile + " > " + sorted_bam_file
				os.popen( cmd )
			os.system(bedtools + " -d -split -ibam " + sorted_bam_file + " > " + args.output_directory + "cov_files/" + name + ".cov")	
		else:
			print(name + ".cov already exists. Continue with further steps...")	

	#calculating read coverage depth per gene
	print("\ncalculating read coverage depth per gene...")
	if os.path.isdir(args.output_directory + "cov_rep_files/") == False:
		os.makedirs(args.output_directory + "cov_rep_files/")
	def get_gene_positions_form_gff( gff_file ):
		gene_positions = []
		with open( gff_file, "r" ) as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if parts[2] == "gene":
						chromosome = parts[0].lower()
						if chromosome in organelles:
							pass	
						elif chromosome == "random" or chromosome.endswith("un"):
							pass
						else:	
							gene_positions.append( { 'chr': chromosome, 'start': int( parts[3] ), 'end': int( parts[4] ), 'ID': parts[-1].split(';')[0].replace( "ID=", "" ) } )
				line = f.readline()
		return gene_positions
	
	def load_coverage( cov_file ):
		coverage = {}
		with open( cov_file, "r" ) as f:
			line = f.readline()
			prev_chr = line.split('\t')[0].lower()
			if prev_chr in organelles:
				pass		
			cov = []
			while line:
				parts = line.strip().split('\t')
				if parts[0].lower() != prev_chr:
					coverage.update( { prev_chr: cov } )
					prev_chr = parts[0].lower()
					if prev_chr in organelles:
						pass	
					cov = []
				cov.append( float( parts[2] ) )
				line = f.readline()
			coverage.update( { prev_chr: cov } )
		return coverage	
		
	def calculate_cov_per_gene( gene_positions, cov ):
		for idx, gene in enumerate( gene_positions ):
			gene_positions[ idx ].update( { 'cov': np.mean( cov[ gene['chr'] ][ gene['start']:gene['end'] ] ) } )
		return gene_positions
		
	gff_pos = get_gene_positions_form_gff(args.gff_file)
	cov_files = glob.glob(args.output_directory + "cov_files/*.cov")
	for covfile in cov_files:
		name = covfile[(len(args.output_directory + "cov_files/")):-4]
		if os.path.isfile(args.output_directory + "cov_rep_files/" + name + "_report.txt") == False:
			coverage = load_coverage(covfile)
			gene_positions = calculate_cov_per_gene(gff_pos, coverage)
			datei = open(args.output_directory + "cov_rep_files/" + name + "_report.txt","w")
			for gene in gene_positions:
				datei.write( "\t".join( map( str, [ gene["ID"], gene['cov'] ] ) ) + '\n' )
			datei.close()
		else:
			print(name + "_report.txt already exists. Continue with further steps...")	
	
	#identify samples with average coverage higher than cutoff (default=10)
	print("\nextracting samples with average coverage higher than cutoff...")
	cov_cutoff = int(args.min_cov)
	if os.path.isfile(args.output_directory + "sample_cov.txt") == False:
		cov_rep_files = glob.glob(args.output_directory + "cov_rep_files/*.txt")
		av_cov = []
		for rep_file in cov_rep_files:
			name = rep_file[(len(args.output_directory + "cov_rep_files/")):-11]
			datei = open(rep_file,"r")
			sample = datei.readlines()
			datei.close()
			covs = []
			for line in sample:
				line = line.strip().split("\t")
				covs.append(float(line[1]))	
			if np.mean(covs) >= cov_cutoff:	
				av_cov.append(name + "\t" + str(np.mean(covs)))
		datei = open(args.output_directory + "sample_cov.txt", "w")
		datei.write("\n".join(av_cov))
		datei.close()			
	else:
		print("sample_cov.txt already exists. Continue with further steps...")	
			
	#construct matrix
	print("\nconstructing matrix...")
	sample_cov = pd.read_csv(args.output_directory + "sample_cov.txt", header = None, sep="\t")
	all_samples = sample_cov[0].to_list()
	if len(all_samples) <= 1:
		print("no samples left after applying coverage cutoff...")
		sys.exit(0)
	dictionary_matrix = {}
	for sample in all_samples:
		filename = (args.output_directory + "cov_rep_files/" + sample + "_report.txt")
		rep_df = pd.read_csv(filename, header = None, sep="\t")
		dictionary_matrix[sample] = rep_df[1].to_list()
		if "gene" not in dictionary_matrix.keys():
			dictionary_matrix["gene"] = rep_df[0].to_list()
	matrix = pd.DataFrame.from_dict(dictionary_matrix)
	matrix = matrix.set_index("gene")
	matrix.to_csv(args.output_directory + "accession_coverage.txt", sep="\t")
	
	#calculating gene dispensability scores
	print("\ncalculating gene dispensability scores...")
	cov_matrix = pd.read_csv(args.output_directory + "accession_coverage.txt", sep="\t")
	N = len(list(cov_matrix)[1:])
	cov_matrix = cov_matrix.set_index("gene")
	cX = cov_matrix.div(cov_matrix.mean(axis="index"), axis="columns")
	ds = (cX.sum(axis="columns")).div(N)
	ds.to_csv(args.output_directory + "gene_dispensability_scores.csv", header=False)

if args.visualize == True:
	print("\nvisualizing...")
	
	#histogram
	plt.yscale("log")
	plt.xlabel("dispensability score (ds)")
	plt.ylabel("number of genes")
	highscores_in_bin2 = []
	for number in sorted(list(ds)):
		if number >= 2:
			value = 2
			highscores_in_bin2.append(value)
		else:
			highscores_in_bin2.append(number)	
	plt.xlim(xmin=0, xmax = 2)
	cm = plt.cm.cool
	n, bins, patches = plt.hist(highscores_in_bin2, bins = np.arange(min(highscores_in_bin2), max(highscores_in_bin2) + 0.02, 0.02))
	for i, p in enumerate(patches):
		plt.setp(p, 'facecolor', cm((i)/(len(bins))))
	plt.savefig(args.output_directory + "score_distribution_hist.png")
	plt.close()
	
	#violin plot
	def drawViolinPlot(xlabel, xticks, xticklabels, ylabel, bandwidth):
		axis.set_xlabel(xlabel);
		axis.set_xticks(xticks);
		axis.set_xticklabels(xticklabels);
		axis.set_ylabel(ylabel);
		axis.violinplot(sequences, showmeans=True, showmedians=False, bw_method=bandwidth);
	sequences = highscores_in_bin2
	figure, axis = plt.subplots(1, 1);
	plt.subplots_adjust(hspace=1);
	bandwidth = None;
	drawViolinPlot("all genes", np.arange(len(highscores_in_bin2)+1), ('all genes'), "dispensability score (ds)", bandwidth);
	bandwidth = 0.3;
	plt.savefig(args.output_directory + "score_distribution_violin.png")
	plt.close()
	
print("\n---finished---")
