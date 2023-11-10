analysis="Q.GTDB_sub"
training="training_loci_5k"
name="5k"

# 5. Estimate the models
log=$training"_log.txt"

echo "Log file for "$training > $log
echo "Estimating initial models with IQ-TREE2 for" $training >> $log

mkdir 02_fullcon_$training # first make the output directory
model_set="LG,Q.pfam,Q.insect,Q.yeast"
iqtree2 -T 128 -st AA -S $training -cmax 4 -mset $model_set -pre 02_fullcon_$training/iteration_1

# get the list of models, and save it to models.txt
grep '^ *[^ ]\+:' 02_fullcon_$training/iteration_1.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr > 02_fullcon_$training/models.txt

echo "List of models best fit to training loci: " >> $log
cat 02_fullcon_$training/models.txt >> $log

# now we get the init model as the first model in that list

initial_model=$(awk 'NR==1 {print $2}' 02_fullcon_$training/models.txt)
echo "Initial Model will be set to" >> $log
echo $initial_model >> $log

# 7. Estimate the Q matrix

echo "Estimating Q matrix with IQ-TREE2" >> $log
echo "Using LG as the initial model due to IQ-TREE bug!" >> $log

iqtree2 -T 128 -S 02_fullcon_$training/iteration_1.best_scheme.nex -te 02_fullcon_$training/iteration_1.treefile --init-model LG --model-joint GTR20+FO -pre 02_fullcon_$training/iteration_1.GTR20


# extract the matrix
new_model="Q.GTDB_sub_"$name
echo "" >> $log
echo "## Model Testing ##" >> $log
echo "Extracting model and saving to "$new_model >> $log

grep -A 21 "can be used as input for IQ-TREE" 02_fullcon_$training/iteration_1.GTR20.iqtree | tail -n20 > $new_model

cat $new_model >> $log


# # test it on test loci
echo "" >> $log
echo "Testing model on test loci..." >> $log
mkdir 03_testing_$training
iqtree2 -T 248 -S testing_loci/ -st AA -m MF -mset LG,Q.pfam,Q.insect,Q.yeast,$new_model -pre 03_testing_$training/test_loci_mf
echo "Frequency table of best-fit models" >> $log
grep '^ *[^ ]\+:' 03_testing_$training/test_loci_mf.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr >> $log
