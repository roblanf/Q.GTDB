# Q.GTDB  
Q matrices for GTDB

## Installation  

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

## Analysis  

### 1. Input data  

#### 1a. Per-locus protein alignments  
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

#### 1b. Taxa  
Next, the pipeline requires an input text file listing all the species that
should be included for training. One species name (i.e. header in the *.faa)
files per line. For example, `data/r207_nitrospinota.taxa` consists of all 62
sequences that belong to the Nitrospinota phylum as of version r207.  
