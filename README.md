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

### 1.1. Per-locus protein alignments  
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

### 1.2. Phylogeny  
This is used to remove fast-evolving sequences (i.e. long branches) that may 
bias the training. For example, the r207 tree (made by GTDB using the filtered
and concatenated alignment using FastTree) `data/gtdb_r207_bac120_unscaled.decorated.tree`

### 1.3. Taxa list
The pipeline requires an input text file listing all the species that should be
included for training. The text file should have one species name (i.e. header
in the *.faa) per-line. For example, `data/r207_nitrospinota.taxa` consists of
all 62 sequences that belong to the Nitrospinota phylum as of version r207.  

## 2. Group-specific Q-matrix  

### 2.1. Subset taxa  

**Preparing the taxa list**  

Prepare the species list for a single taxonomic rank (i.e. the Nitrospinota
phylum) using a taxonomy .tsv such as `gtdb_r214_selected_genomes.bacteria.family.tsv`.

**Removing long branches**  

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

**Subset taxa**  
Create new taxa list with dropped taxa:  
```bash
# Assuming that these are all genbank accessions that start with G..
grep -oP "G\d+" pruned_treeshrink/output.tree > treeshrunk.taxa
```

Subset relevant taxa first:
```bash
mkdir -p 00_subset_taxa
for loc in r207_loci/*.faa; do
	faSomeRecords.py --fasta $loc --list treeshrunk.taxa --outfile 00_subset_taxa/${loc}
done
```  

### 2.2. Assigning loci for training and testing  

Randomly assign 20 loci for testing, and 100 for training:  
```bash
scripts/assign_loci_random.sh 00_subset_taxa
```

### 2.3. Initial model selection  
Identify the best models for each locus to determine the starting models for
training:
```bash
mkdir -p 01_model_selection
iqtree2 -S 00_subset_taxa/training_loci -T 4 -pre 01_model_selection/training_loci
iqtree2 -S 00_subset_taxa/testing_loci -T 4 -pre 01_model_selection/testing_loci
cat 01_model_selection/*.best_scheme > 01_model_selection/combined.best_scheme
```  

This scripts counts the frequency of best-fitting models across all loci.
Then, it selects the models in the top 90% (default cut-off) as the starting
models `-mset` for training:  
```bash
scripts/count_top_models.py 01_model_selection/combined.best_scheme > 01_model_selection/starting_models.txt
```

**Example output**  

> [('Q.yeast', 73), ('LG', 27), ('Q.insect', 8), ('Q.pfam', 7), ('HIVw', 1), ('mtZOA', 1)]  
> Total loci: 117  
> Cut-off: 0.9  
> Selecting most frequent models up to 106 loci.  
> Starting models: Q.yeast,LG,Q.insect  

### 2.4. Training  

**Training modes**  

| Mode              | Name    | Fixed topology? | Linked branch lengths? |
| ----------------- | ------- | --------------- | ---------------------- |
| Unconstrained     | uncon   | No              | No                     |
| Semi-constrained  | semicon | No              | Yes                    |
| Fully-constrained | fullcon | Yes             | Yes                    |

**Fully constrained**  
```bash
scripts/estimate_q.py \ 
        --mode fullcon \
        --loci 00_subset_taxa/training_loci/ \
        -mset Q.yeast,LG,Q.insect \
        -te pruned_treeshrink/output.tree \
        -T 4 -v
```

**Semi-constrained**  
```bash
scripts/estimate_q.py \ 
        --mode semicon \
        --loci 00_subset_taxa/training_loci/ \
        -mset Q.yeast,LG,Q.insect \
        -T 4 -v
```

**Unconstrained**  
```bash
scripts/estimate_q.py \ 
        --mode uncon \
        --loci 00_subset_taxa/training_loci/ \
        -mset Q.yeast,LG,Q.insect \
        -T 4 -v
```

### 2.5. Testing  
Run ModelFinder on the test loci, using the best starting models identified,
and the estimated Q-matrices.

For example:  
```bash
mkdir -p 03_testing
iqtree2 -T 4 -S 00_subset_loci/testing_loci/ -m MF -pre 03
_testing/test_loci_mf -mset LG,Q.pfam,Q.insect,02_uncon/Q.uncon_i2,02_semicon/Q.s
emicon_i2,02_fullcon/Q.fullcon_i2
```

## 3. Global Q matrix  

### 3.1. Subset taxa 


