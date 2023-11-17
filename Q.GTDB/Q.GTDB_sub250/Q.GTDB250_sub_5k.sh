analysis="Q.GTDB_sub250"
training="training_loci_5k"
name="5k"
new_model=$analysis"_"$name
log=$training"_log.txt"

echo "Log file for "$training > $log
echo "Estimating initial models with IQ-TREE2 for" $training >> $log

mkdir 02_fullcon_$training # first make the output directory
model_set="LG,Q.pfam,Q.insect,Q.yeast"
iqtree2 -T 128 -st AA -S $training -cmax 8 -mset $model_set -pre 02_fullcon_$training/iteration_1

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


# test it on the full and the reduced datasets
# we don't need to do LG since we did that for phylum_1...

# small alignment
raxml-ng --msa ../concat_alignments/gtdb_r207_bac120_concatenated.faa --model PROTGTR{$new_model}+G --threads 16 --force perf_threads --tree r207_original_clean.tree --evaluate --lh-epsilon 0.1  --prefix 03_testing_$training/$new_model"_G_reduced_aln"

# big alignment
raxml-ng --msa ../concat_alignments/gtdb_r207_bac120_full.faa --model PROTGTR{$new_model}+G --threads 16 --force perf_threads --tree r207_original_clean.tree --evaluate --lh-epsilon 0.1  --prefix 03_testing_$training/$new_model"_G_full_aln"


# finally we make a little table of the likelihoods etc. for the analyses, and compare it to LG

# Specify your list of filenames
f1=03_testing_$training/$new_model"_G_reduced_aln.raxml.log"
f2="../phylum_1/03_testing/LGG.raxml.log"
f3=03_testing_$training/$new_model"_G_full_aln.raxml.log"
f4="../phylum_1/03_testing/LGG_full.raxml.log"
declare -a filenames=($f1 $f2 $f3 $f4)

# Create a new file to store the table
log_file="03_testing_"$training"/log_table.txt"

# Print the header of the table
echo -e "Likelihood\tAIC\tTime\tFilename" > "$log_file"

# Loop through each specified file
for file in "${filenames[@]}"; do
    # Check if the file exists
    if [[ -f "$file" ]]; then
        # Extract the required values using grep and awk
        likelihood=$(grep 'Final LogLikelihood' "$file" | awk '{print $3}')
        aic=$(grep 'AIC score' "$file" | awk '{print $3}')
        time=$(grep 'Elapsed time' "$file" | awk '{print $3}')
        filename=$(basename "$file")

        # Append the extracted values to the log file
        echo -e "$likelihood\t$aic\t$time\t$filename" >> "$log_file"
    else
        echo "The file $file does not exist." >> "$log_file"
    fi
done

# Display the table
cat "$log_file"

echo "" >> $log
echo "Comparing fit of new model on full alignment and r207 tree" >> $log

cat $log_file >> $log

