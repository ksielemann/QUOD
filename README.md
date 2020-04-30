# QUOD
Reference-based QUantification Of gene-Dispensability
![alt text](https://github.com/KatharinaFrey/QUOD/blob/master/QUOD_logo.png)


### Usage:

  python3 QUOD.py
  --in STR*				full path to folder with input bam files
  --gff STR* 			full path to reference annotation file
  --out STR*			full path to output folder

  Optional:
  --min_cov_per_genome INT* 	default = 10
  --bam_is_sorted 			prevents extra sorting of BAM files
  --visualize


REQUIREMENTS: pandas, numpy (optional: matplotlib for visualization)


![Illustration of the QUOD method using an artificial dataset. On the left side, genes are classified as ‘core’ or ‘dispensable’ according to a cutoff. On the right side, gene dispensability is quantified according to a dispensability score based on the normalised coverage in a read mapping (a-f: investigated genomes). Coloring of genes (right side) indicates different dispensability scores. Extremely rare genes can be easily detected using QUOD.](https://github.com/KatharinaFrey/QUOD/blob/master/QUOD_concept.png)


### Reference:
Frey K., Weisshaar B., Pucker B.; Reference-based QUantification Of gene Dispensability (QUOD); bioRxiv 2020.04.28.065714; doi: <https://doi.org/10.1101/2020.04.28.065714>
