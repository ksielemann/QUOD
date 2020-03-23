# QUOD
Reference-based QUantification Of gene-Dispensability
![alt text](https://github.com/KatharinaFrey/QUOD/blob/master/QUOD_logo.jpg)


USAGE:\n
python3 QUOD.py

--input_dir <FULL_PATH_TO_FOLDER_INPUT_BAM_FILES> (file names = sample names)

--bam_is_sorted <PREVENTS_EXTRA_SORTING_OF_BAM_FILES> (optional argument)

--gff <FULL_PATH_TO_REFERENCE_ANNOTATION_FILE>

--output_dir <FULL_PATH_TO_OUPUT_FOLDER>

--min_cov_per_genome <INTEGER> (default = 10, optional argument)

--visualize (optional argument)

REQUIREMENTS: os, glob, sys, argparse, pandas, numpy (optional: matplotlib for visualization)
