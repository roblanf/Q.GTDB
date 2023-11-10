analysis="Q.GTDB_sub250"

echo "Setting up loci for "$analysis > log.txt

# 3. just get the taxa we want from the loci

echo "Subsetting alignments" >> log.txt

mkdir -p loci
for taxon_list in taxon_lists/*.txt; do
  
    base_name=$(basename "$taxon_list" .txt)
  
    for loc in ../alignments/*.faa; do
        filename=$(basename $loc)
        new_filename="${base_name}_${filename}"

        faSomeRecords $loc $taxon_list loci/${new_filename}
        
    done
done

echo "" >> log.txt
echo "Sub-alignments created for each gene and taxon list." >> log.txt
alignment_count=$(find loci/ -name "*.faa" | wc -l)
echo "A total of $alignment_count alignments were created." >> log.txt

# 4. Split the loci between training and testing

echo "splitting alignments into testing and training" >> log.txt


test_set=$(ls loci | sort -R | tail -1000)

echo "Test alignments: " >> log.txt
echo $test_set >> log.txt
mkdir -p testing_loci
for file in $test_set; do
    mv "loci/$file" testing_loci/
done

# now we need to de-duplicate the remaining loci. No training comes from duplicates...
# we then remove any files with <10 sequences
# we also remove any sequences that are all gaps

mkdir loci_deduped
mkdir loci_clean
for file in loci/*; do
    fname=$(basename $file)
    seqkit rmdup --by-seq -o loci_deduped/$fname $file
    seqkit grep -w 0 -svrp "^-+$" loci_deduped/$fname > loci_clean/$fname # remove sequences that are all gaps
    nseq=$(grep -c '>' loci_clean/$fname)  # -c will count the number of matches
    echo "$nseq sequences left"
    if [ "$nseq" -lt 10 ]; then
        rm loci_clean/$fname
    fi    
done

# clean up!
rm -rf loci_deduped

mkdir -p training_loci_1k
mkdir -p training_loci_5k
mkdir -p training_loci_10k
training_set1k=$(ls loci_clean | sort -R | tail -1000)
training_set5k=$(ls loci_clean | sort -R | tail -5000)
training_set10k=$(ls loci_clean | sort -R | tail -10000)


for file in $training_set1k; do 
    cp "loci_clean/$file" training_loci_1k/
done

for file in $training_set5k; do 
    cp "loci_clean/$file" training_loci_5k/
done

for file in $training_set10k; do 
    cp "loci_clean/$file" training_loci_10k/
done

# check! 

echo "" >> log.txt
echo "Number of testing loci : " >> log.txt
ls testing_loci/ | wc -l >> log.txt

echo "" >> log.txt
echo "Number of training loci : " >> log.txt
alignment_count=$(find loci/ -name "*.faa" | wc -l)

n1k=$(find training_loci_1k/ -name "*.faa" | wc -l)
n5k=$(find training_loci_5k/ -name "*.faa" | wc -l)
n10k=$(find training_loci_10k/ -name "*.faa" | wc -l)
echo "training_loci_1k: $n1k" >> log.txt
echo "training_loci_5k: $n5k" >> log.txt
echo "training_loci_10k: $n10k" >> log.txt