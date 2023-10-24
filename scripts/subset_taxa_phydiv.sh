#!/bin/bash
set -e

# Replace hardcodes
IQTREE=/home/frederickjaya/Downloads/iqtree-2.2.2.7-Linux/bin/iqtree2
TREE=/home/frederickjaya/GitHub/gtdb_trees/2306_global/unscaled_shrunk0.05/output.tree
LOCI_DIR=/home/frederickjaya/Dropbox/gtdb/01_data/gtdb_r207_full_concat

# Modify the below line to select k taxa 
for k in {25,50,100,250,500,1000,2500,5000}; do
	# Subset k taxa with the maximum phylogenetic diversity
	mkdir -p $k
	$IQTREE -te $TREE -k $k -pre $k/unscaled_shrunk0.05_pd${k}
	# Parse taxa names
	sed '0,/The optimal PD set has/d' $k/unscaled_shrunk0.05_pd${k}.pda | sed '/^$/,$d' > ${k}/unscaled_shrunk0.05_pd${k}.txt
	# Subset taxa from each loci (this is done in nf)  
	#mkdir -p $k/loci
	#for l in $LOCI_DIR/*.faa; do
	#	echo $l
	#	faSomeRecords.py --fasta $l --list ${k}/unscaled_shrunk0.05_pd${k}.txt --outfile $k/loci/`basename $l`
	#done
	for i in $k/unscaled_shrunk0.05_pd${k}.pda; do
		# Parse subtree to separate tree file
		grep -A1 "Corresponding sub-tree" $i | tail -1 > \
			$k/`basename $i .pda`.tree
		done
done
