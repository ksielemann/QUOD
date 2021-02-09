[![DOI](https://zenodo.org/badge/249387105.svg)](https://zenodo.org/badge/latestdoi/249387105)

# QUOD


## Reference-based QUantification Of gene-Dispensability
<img src="https://github.com/ksielemann/QUOD/blob/master/QUOD_logo.png" alt="drawing" width="500"/>


### Background:

Dispensability of genes in a phylogenetic lineage, e.g. a species, genus, or higher-level clade, is gaining relevance as most genome sequencing projects move to a pangenome level.  
Most analyses classify genes as core genes, which are present in (almost) all investigated individual genomes, and dispensable genes, which only occur in a single or a few investigated genomes. The binary classification as ‘core’ or ‘dispensable’ is often based on arbitrary cutoffs of presence/absence in the analysed genomes.  
Instead of classifying a gene as core or dispensable, QUOD assigns a dispensability score to each gene. Hence, QUOD facilitates the identification of candidate dispensable genes which often underlie lineage-specific adaptation to varying environmental conditions.  


### Concept:

<p align="center">
<img src="https://github.com/ksielemann/QUOD/blob/master/QUOD_concept.png" alt="drawing" width="750"/>
</p>

Illustration of the QUOD method using an artificial dataset. On the left side, genes are classified as ‘core’ (dark blue) or ‘dispensable’ (dark grey) according to a cutoff. On the right side, gene dispensability is quantified according to a dispensability score based on the normalised coverage in a read mapping. Coloring of genes (right side) indicates different dispensability scores. Extremely rare genes can be easily detected using QUOD.


### Usage:

Genomic reads of the investigated genomes are mapped to the corresponding reference genome sequence. The resulting BAM files of these mappings are then subjected to QUOD. In addition, the corresponding annotation file (in GFF format) has to be provided. '--out' is used to specify the output folder. A gene-specific dispensability score for each gene based on normalised coverage in read mappings is calculated [1].

Optionally, genomes with an average coverage below a given cutoff (default=10) can be discarded and excluded from further analyses using the '--min_cov_per_genome' parameter. The optional parameter '--bam_is_sorted' prevents extra sorting of the BAM files in case these files are already sorted. '--visualize' activates the visualization module which constructs a histogram and a box plot of the dispensability score distribution.  

~~~
python3 QUOD.py
  
--in STR        full path to folder with input bam files

--gff STR       full path to reference annotation file
  
--out STR       full path to output folder

  
Optional:
  
--min_cov_per_genome INT      default = 10
  
--bam_is_sorted               prevents extra sorting of BAM files
  
--visualize
~~~

- REQUIREMENTS: pandas, numpy (optional: matplotlib for visualization)


#### Usage example (test set):
The test dataset for QUOD comprises genomic reads of four randomly selected accessions of the *Arabidopsis thaliana* Nordborg set. The reads were received from the Sequence Read Archive (SRA) and mapped against the AthNd1_v2c reference genome sequence [2] using bowtie2 [3]. To reduce the size of the files, the first MB of NdCChr1 was extracted. All BAM files provided here are already sorted and should be used as input for QUOD. The test dataset including all relevant files can be downloaded from 'PUB – Publications at Bielefeld University' (doi: <https://doi.org/10.4119/unibi/2946079>).  

~~~
python3 QUOD.py --in /input_bams_testset/ --bam_is_sorted --gff AthNd1_v2c_chr1_1mb.gff3 --out /output_QUOD/ --visualize
~~~


### Output:
A dispensability score is calculated for each gene. Optionally, the results can be visualized as a colored histogram and a box plot.

<p align="center">
<img src="https://github.com/ksielemann/QUOD/blob/master/score_distribution.png" alt="drawing" width="400"/>
</p>


### Score composition of selected genes:
For further investigation of the score composition of selected genes-of-interest, the script 'score_composition.py' can be used:
~~~
python3 score_composition.py
--gene_file              <FULL_PATH_TO_GENE_FILE> (one gene per line (gene name from gff file used for QUOD))
--score_file             <FULL_PATH_TO_gene_dispensability_scores_FILE> (QUOD output file 'gene_dispensability_scores.csv')
--accession_cov_file     <FULL_PATH_TO_accession_coverage_FILE> (QUOD output file 'accession_coverage.txt')
--output_dir             <FULL_PATH_TO_OUPUT_FOLDER>

Optional:
--visualize
~~~
As output, a table including  
(I) the dispensability score,  
(II) the average coverage of all investigated genome sequences,  
(III) the average coverage of the accessions with the highest and (IV) lowest 10 % of all coverage values, respectively,  
(V) the number of accessions with zero coverage and  
(VI) the coverage for each accession, separately, is provided.  
Further, the coverage distribution for each gene can be visualized in a box plot.  


### Statistical analyses:
The following scripts, performed to evaluate QUOD statistically, are provided:  

(A) variance_in_repl_test.py          (Investigation of the variance within replicate datasets of the same accession in comparison to iteratively, randomly selected subsets of the whole dataset.)  
(B) analyse_TE_genes_and_scores.py    (Welch's t-test to compare dispensability scores of TE genes and non-TE genes.)  
(C) mixed_linear_effects.py           (Mixed linear modelling to compare the dispensability score with either gene length, exon number or distance to closest TE gene.)  
(D) BUSCO_comparison.py               (Welch's t-test to compare dispensability scores of BUSCO genes and non-BUSCO genes. Levene's test to investigate the variation of the dispensability scores of BUSCO and non-BUSCO genes.)  

More details for these analyses can be found in the Method section of our publication.  


### How to cite:
Sielemann, K., Weisshaar, B. & Pucker, B. Reference-based QUantification Of gene Dispensability (QUOD). Plant Methods 17, 18 (2021). <https://doi.org/10.1186/s13007-021-00718-5> <br/><br/><br/>  


##### References:  
<sup>[1] Pucker B, Brockington SF. Genome-wide analyses supported by RNA-Seq reveal non-canonical splice sites in plant genomes. BMC Genomics 19, 980 (2018). <https://doi.org/10.1186/s12864-018-5360-z></sup>  

<sup>[2] Pucker B, et al. A chromosome-level sequence assembly reveals the structure of the Arabidopsis thaliana Nd-1 genome and its gene set. PloS one 14.5 (2019): e0216233. <https://doi.org/10.1371/journal.pone.0216233></sup>  

<sup>[3] Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359. <https://doi.org/10.1038/nmeth.1923></sup>
