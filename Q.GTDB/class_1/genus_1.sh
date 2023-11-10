analysis="genus_1"

# 1. set up
mkdir $analysis
cd $analysis

echo ""
echo "Setting up analysis for "$analysis > log.txt

cp ../$analysis.txt .

# 2. get the subtree (produces $analysis.tree)

echo ""
echo "Getting subtree for "$analysis".txt taxon list" >> log.txt

../scripts/get_subtree.py ../r207_original_clean.tree $analysis.txt 


# 3. just get the taxa we want from the loci

echo ""
echo "Subsetting alignments" >> log.txt

mkdir -p loci
for loc in ../alignments/*.faa; do
    filename=$(basename $loc)
    faSomeRecords $loc $analysis.txt loci/${filename}
done

# remove the full alignments, we don't want them
rm loci/gtdb_r207_bac120_full.faa 
rm loci/gtdb_r207_bac120_concatenated.faa 


# 4. Split the loci between training and testing


echo ""
echo "splitting alignments into testing and training" >> log.txt

cd loci
test_set=$(ls | sort -R | tail -20)

echo ""
echo "Test alignments: " >> log.txt
echo $test_set >> log.txt

mkdir -p training_loci
mkdir -p testing_loci


mv $test_set testing_loci
mv *.faa training_loci
cd ..

# check! 

echo ""
echo "Number of training loci : " >> log.txt
ls loci/training_loci/ | wc -l >> log.txt

echo ""
echo "Number of testing loci : " >> log.txt
ls loci/testing_loci/ | wc -l >> log.txt

# 5. Estimate the models

echo ""
echo "Estimating initial models with IQ-TREE2" >> log.txt

mkdir 02_fullcon # first make the output directory
model_set="LG,Q.pfam,Q.insect,Q.yeast"
iqtree2 -T 100 -p loci/training_loci -m MFP -cmax 8 -mset $model_set -te $analysis.tree -pre 02_fullcon/iteration_1

# get the list of models, and save it to models.txt
grep '^ *[^ ]\+:' 02_fullcon/iteration_1.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr > 02_fullcon/models.txt

echo ""
echo "List of models best fit to training loci: " >> log.txt
cat 02_fullcon/models.txt >> log.txt

# now we set the initial model to LG, because it's usually best for these data
# and for whatever reason, if it's Q.pfam, there's a bug in IQ-TREE that means it doesn't work!

initial_model="LG"
echo ""
echo "Initial Model will be set to" >> log.txt
echo $initial_model >> log.txt

# 7. Estimate the Q matrix

echo ""
echo "Estimating Q matrix with IQ-TREE2" >> log.txt

iqtree2 -T 100 -S loci/training_loci -p 02_fullcon/iteration_1.best_scheme.nex -te 02_fullcon/iteration_1.treefile --init-model $initial_model --model-joint GTR20+FO -pre 02_fullcon/iteration_1.GTR20


