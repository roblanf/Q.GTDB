# Q.GTDB  
Q matrices for GTDB

# Installation  

Download repo:  
```bash
git clone https://github.com/roblanf/Q.GTDB.git
```  

Install conda environment:  
```bash
cd Q.GTDB/
conda env create -f env.yml
conda activate qgtdb
```

Install faSomeRecords and AliStat:

```bash
cd bin/ 
# faSomeRecords
wget https://raw.githubusercontent.com/santiagosnchez/faSomeRecords/master/faSomeRecords.py
chmod +x faSomeRecords.py

# AliStat
wget https://github.com/thomaskf/AliStat/archive/refs/tags/v1.14.tar.gz
tar -xzvf v1.14.tar.gz
rm v1.14.tar.gz
cd AliStat-1.14/
make
mv alistat ../
cd ../
rm -r AliStat-1.14/
```

# Analysis  

## 1. Input data  
This section is an overview of the input files required. The following sections
will include information on how to generate these depending on the analysis.  

### 1a. Per-locus protein alignments  
The pipeline has been developed so that the Q-matrices are trained on the 
whole, 120 locus data set. Each locus is its own file, with all taxa included.
All locus files (*.faa) should be in the same directory, for example:  
```bash
r207_loci/
├── gtdb_r207_bac120_PF00380.20.faa
├── gtdb_r207_bac120_PF00410.20.faa
├── ...
└── gtdb_r207_bac120_TIGR03953.faa
```  

### 1b. Phylogeny  
This is used to remove fast-evolving sequences (i.e. long branches) that may 
bias the training. For example, the r207 tree (made by GTDB using the filtered
and concatenated alignment using FastTree) `data/gtdb_r207_bac120_unscaled.decorated.tree`

### 1c. Taxa list
The pipeline requires an input text file listing all the species that should be
included for training. The text file should have one species name (i.e. header
in the *.faa) per-line. For example, `data/r207_nitrospinota.taxa` consists of
all 62 sequences that belong to the Nitrospinota phylum as of version r207.  

## 2. Training: group-specific  

### 2a. Preparing the taxa list  
Prepare the species list for a single taxonomic rank (i.e. the Nitrospinota
phylum) using a taxonomy .tsv such as `gtdb_r214_selected_genomes.bacteria.family.tsv`.  

### 2b. Removing long branches  
Exclude fast-evolving species (i.e. long branches) that may bias training.  

Subset the group-specific tree from the full reference tree:  
```bash
scripts/get_subtree.py [ref_tree]
# output: pruned.tree
```  

Output the initial total tree length and distribution of branch lengths:  
```bash
scripts/tree_length.py pruned.tree
```

The tree length will be output to stdout and the histogram to
`branch_length_histogram.png`.  

Remove long branches with treeshrink and output branch lengths again:  
```bash
run_treeshrink -t pruned.tree
scripts/tree_length.py pruned_treeshrink/output.tree
```  

Identify the number of taxa removed:  
```bash
grep -oP "\t" pruned_treeshrink/output.txt | wc -l
```

Create new taxa list with dropped taxa:  
```bash
# Assuming that these are all genbank accessions that start with G..
grep -oP "G\d+" pruned_treeshrink/output.tree > treeshrunk.taxa
```

### 2c. Assigning loci for training and testing  
Subset relevant taxa first:
```bash
mkdir -p 00_subset_taxa
for loc in r207_loci/*.faa; do
	faSomeRecords.py --fasta $loc --list treeshrunk.taxa --outfile 00_subset_taxa/${loc}
```  

### 2d. Assigning loci for training and testing  

Randomly assign 20 loci for testing, and 100 for training:  
```
scripts/get_subtree.py 00_subset_taxa
```
