analysis="Q.GTDB_sub"
training="training_loci_5k"

# 5. Estimate the models

echo "Estimating initial models with IQ-TREE2 for $training" >> log.txt

mkdir 02_fullcon_$training # first make the output directory
model_set="LG,Q.pfam,Q.insect,Q.yeast"
iqtree2 -T 200 -st AA -S $training -cmax 4 -mset $model_set -pre 02_fullcon_$training/iteration_1

# get the list of models, and save it to models.txt
grep '^ *[^ ]\+:' 02_fullcon_$training/iteration_1.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr > 02_fullcon_$training/models.txt

echo "List of models best fit to training loci: " >> log.txt
cat 02_fullcon_$training/models.txt >> log.txt

# now we get the init model as the first model in that list

initial_model=$(awk 'NR==1 {print $2}' 02_fullcon_$training/models.txt)
echo "Initial Model will be set to" >> log.txt
echo $initial_model >> log.txt

# 7. Estimate the Q matrix

echo "Estimating Q matrix with IQ-TREE2" >> log.txt

iqtree2 -T 200 -S 02_fullcon_$training/iteration_1.best_scheme.nex -te 02_fullcon_$training/iteration_1.treefile --init-model LG --model-joint GTR20+FO -pre 02_fullcon_$training/iteration_1.GTR20
