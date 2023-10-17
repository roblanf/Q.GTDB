# Q.GTDB  
Q matrices for GTDB

# Installation  

Download repo:  
```
git clone https://github.com/roblanf/Q.GTDB.git
```  

Install conda environment:  
```
cd Q.GTDB/
conda env create -f env.yml
conda activate qgtdb
```

Install faSomeRecords and AliStat:

```
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
```
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
```
scripts/get_subtree.py [ref_tree]
# output: pruned.tree
```  

Output the initial total tree length and distribution of branch lengths:  
```
scripts/tree_length.py pruned.tree
```

The tree length will be output to stdout and the histogram to
`branch_length_histogram.png`.  

Remove long branches with treeshrink and output branch lengths again:  
```
run_treeshrink -t pruned.tree
scripts/tree_length.py pruned_treeshrink/output.tree
```  
