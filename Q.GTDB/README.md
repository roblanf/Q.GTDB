# Q.GTDB  

The aim for this matrix is to estimate a single amino acid matrix for the GTDB dataset.

The dataset comprises 120 genes, each for ~85K taxa. The main challenges are that this is just too many taxa. So we need to come up with sensible ways of estimating matrices from a subset of the taxa. I will use a lot of the suggestions from https://github.com/roblanf/Q.GTDB, with some modifications.

## Getting set up

Install and activate the conda environment like so:

```{bash}
conda env create -f env.yml
conda activate qgtdb
```

## Input data

For the purposes of this work, the 120 gene alignments are in:

```
/alignments
```

these aren't on GitHub - they're too big.

The pre-estimated tree for these data (estimated in FastTree using an LG model) is:

```
gtdb_r207_bac120_unscaled.decorated.tree
```

there is also a 'clean' version of this tree, without the decorations:

```
r207_original_clean.tree
```

finally, there's a big taxonomy table here:

```
gtdb_r207_bac120_curation_taxonomy.tsv
```

This file is a bit annoying. It looks like this:

```
G009834515	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__Pseudomonas_E sp009834515
G900187605	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__Pseudomonas_E sp900187605
G013433315	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__Pseudomonas_E crudilactis
G018138145	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__Pseudomonas_E koreensis_A
```

So we first process it into a taxonomy file that's tab delimited for easier processing:

```
awk 'BEGIN { FS=OFS="\t" } { gsub(";", "\t", $2) } 1' input.tsv > output.tsv
```

now it looks like this:

```
G009834515	d__Bacteria	p__Proteobacteria	c__Gammaproteobacteria	o__Pseudomonadales	f__Pseudomonadaceae	g__Pseudomonas_E	s__Pseudomonas_E sp009834515
G900187605	d__Bacteria	p__Proteobacteria	c__Gammaproteobacteria	o__Pseudomonadales	f__Pseudomonadaceae	g__Pseudomonas_E	s__Pseudomonas_E sp900187605
G013433315	d__Bacteria	p__Proteobacteria	c__Gammaproteobacteria	o__Pseudomonadales	f__Pseudomonadaceae	g__Pseudomonas_E	s__Pseudomonas_E crudilactis
G018138145	d__Bacteria	p__Proteobacteria	c__Gammaproteobacteria	o__Pseudomonadales	f__Pseudomonadaceae	g__Pseudomonas_E	s__Pseudomonas_E koreensis_A
```

## Taxonomy summary

in R...

```{R}
library(tidyverse)
tax <- read_delim("taxonomy.tsv", delim="\t", col_names = FALSE)
names(tax) <- c("id", "clade", "phylum", "class", "order", "family", "genus", "species")

unique_counts <- tax %>%
    gather(column, value) %>%                # Convert data to long format
    group_by(column) %>%                     # Group by column names
    summarise(unique_entries = n_distinct(value)) %>% # Count unique values
    arrange(unique_entries)
```

Gives:

```
> unique_counts
# A tibble: 8 × 2
  column  unique_entries
  <chr>            <int>
1 clade                1
2 phylum             169
3 class              428
4 order             1460
5 family            3650
6 genus            15342
7 id               62291
8 species          62291
```

### Assessing the distribution of gaps

Since we want to estimate substitutino models, we don't want to include genomes with too many gaps. Let's take a look at that first.

```
gtdb_r207_bac120_full.faa
```

I can then calculate gap proportions like this (for sure there are quicker and smarter ways using actual software, but I was about to go home for the evening, so this sufficed):

```{bash}
awk '/^>/ {if (seq) {print id, seq}; id=substr($1, 2); seq=""; next} {seq=seq$0} END {print id, seq}' gtdb_r207_bac120_full.faa | 
while read -r id seq; do
    total_length=$(echo -n "$seq" | wc -c)
    gap_count=$(echo -n "$seq" | sed 's/[^-]//g' | wc -c)
    gap_proportion=$(echo "scale=2; $gap_count / $total_length" | bc -l)
    echo -e "$id\t$gap_proportion"
done > gaps.tsv
```

We can look at the distribution:

```{r}
gaps <- read_delim("gaps.tsv", delim="\t", col_names = c("id", "proportion_gaps"))
ggplot(gaps, aes(x = proportion_gaps)) + geom_histogram(bins=80)
```

Let's add the gaps to the big dataset:

```{r}
tax <- tax %>%
  inner_join(gaps, by = "id")
```

Let's look at the proportion of sequences with fewer gaps than each threshold:

```{r}
results <- map_dfr(seq(0, 1, by = 0.01), function(threshold) {
  proportion_below_threshold <- mean(tax$proportion_gaps < threshold)
  tibble(threshold = threshold, proportion_below = proportion_below_threshold)
})

head(results, 40)

ggplot(results, aes(x = threshold, y = proportion_below)) +
  geom_point() +
  geom_line()
```

This shows that we lose half the data by removing everything with >10% gaps. Let's check how different thresholds would affect the number of phyla and other ranks we'd have.

```{r}

thresholds <- seq(100, 0, by = -10)

for (thresh in thresholds) {
  col_name <- paste0("unique_entries_", thresh)
  
  temp_counts <- tax %>%
    filter(proportion_gaps <= (thresh / 100)) %>%
    gather(column, value, -id, -proportion_gaps) %>%
    group_by(column) %>%
    summarise(!!col_name := n_distinct(value))
  
  unique_counts <- left_join(unique_counts, temp_counts, by = "column")
}

percentage_retained <- unique_counts[,3:12] / unique_counts$unique_entries

rownames(percentage_retained) <- unique_counts$column

```

This shows us what proportion of each rank (phylum, class, etc) is retained with different gap cutoffs:

```
> percentage_retained
        unique_entries_100 unique_entries_90 unique_entries_80 unique_entries_70 unique_entries_60 unique_entries_50 unique_entries_40 unique_entries_30 unique_entries_20 unique_entries_10
clade                    1                 1                 1         1.0000000         1.0000000         1.0000000         1.0000000         1.0000000         1.0000000         1.0000000
phylum                   1                 1                 1         1.0000000         1.0000000         0.9881657         0.9822485         0.9644970         0.8579882         0.5798817
class                    1                 1                 1         1.0000000         1.0000000         0.9953271         0.9836449         0.9485981         0.8247664         0.5490654
order                    1                 1                 1         1.0000000         1.0000000         0.9952055         0.9712329         0.9082192         0.7664384         0.4910959
family                   1                 1                 1         0.9997260         0.9986301         0.9901370         0.9602740         0.8660274         0.7216438         0.4646575
genus                    1                 1                 1         0.9997393         0.9986312         0.9866380         0.9431626         0.8473472         0.7088385         0.4561987
id                      NA                NA                NA                NA                NA                NA                NA                NA                NA                NA
species                  1                 1                 1         0.9998234         0.9983465         0.9822286         0.9328153         0.8389013         0.6966335         0.4836814
```

So a sensible cutoff here is to only keep taxa with <20% gaps - this means we keep ~70% of the species, and more than 80% of the phyla and classes are still represented.

### Assessing long branches

We  also need to assess if any taxa are on crazy long branches. This can happen if the genomes are poorly sequenced or aligned, and can be a big problem for estimating the Q matrix. So let's use TreeShrink for this.

```
run_treeshrink.py -t r207_original_clean.tree
```

Now let's look at what happened, by analysing the tree before and after:

```
scripts/tree_length.py r207_original_clean.tree
mv branch_length_histogram.png before_pruning.png
scripts/tree_length.py r207_original_clean_treeshrink/output.tree
mv branch_length_histogram.png after_pruning.png
```

Before

![](before_pruning.png)

After

![](after_pruning.png)

The total tree length dropped from 6908.397739999869 to 6899.0687299998635. 

Treeshrink removed 18 branches, by removing the following taxa:

```
G002325765      
G001730065      
G002682985      
G000508245      
G000277795      
G000319365      
G000281235      
G000477415      
G000179035      
G002135175      
G000238995 
G000200735       
G001645765      
G902712995      
G000008205      
G000815065      
G900660545      
G001951355      
G000218525
```

We can look at them like this

```{r}
treeshrink <- c("G002325765","G001730065","G002682985","G000508245","G000277795","G000319365","G000281235","G000477415","G000179035","G002135175","G000238995","G000200735","G001645765","G902712995","G000008205","G000815065","G900660545","G001951355","G000218525")

tax %>% filter(id %in% treeshrink)
```

```
> tax %>% filter(id %in% treeshrink)
# A tibble: 19 × 9
   id         clade       phylum            class                  order                 family                  genus               species                           proportion_gaps
   <chr>      <chr>       <chr>             <chr>                  <chr>                 <chr>                   <chr>               <chr>                                       <dbl>
 1 G002682985 d__Bacteria p__Proteobacteria c__Gammaproteobacteria o__Enterobacterales_A f__Enterobacteriaceae_A g__Tremblaya        s__Tremblaya phenacola                       0.62
 2 G001730065 d__Bacteria p__Proteobacteria c__Alphaproteobacteria o__Rickettsiales      f__AB1-6                g__AB1-6            s__AB1-6 sp001730065                         0.52
 3 G002325765 d__Bacteria p__Proteobacteria c__Alphaproteobacteria o__Rickettsiales      f__UBA1459              g__UBA1459          s__UBA1459 sp002325765                       0.66
 4 G900660545 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Metamycoplasmataceae g__Mycoplasmopsis_A s__Mycoplasmopsis_A cynos                    0.34
 5 G000008205 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Metamycoplasmataceae g__Mesomycoplasma   s__Mesomycoplasma hyopneumoniae              0.34
 6 G000815065 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Metamycoplasmataceae g__Mesomycoplasma   s__Mesomycoplasma flocculare                 0.33
 7 G000218525 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Metamycoplasmataceae g__Mesomycoplasma   s__Mesomycoplasma ovipneumoniae_A            0.34
 8 G001951355 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Metamycoplasmataceae g__Mesomycoplasma   s__Mesomycoplasma ovipneumoniae              0.34
 9 G000277795 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A wenyonii                  0.56
10 G000508245 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A ovis                      0.57
11 G000319365 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A haemominutum              0.56
12 G000281235 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A haemolamae                0.57
13 G000179035 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A suis                      0.54
14 G000477415 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_A  s__Eperythrozoon_A parvum                    0.54
15 G000200735 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_B  s__Eperythrozoon_B haemofelis                0.44
16 G000238995 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_B  s__Eperythrozoon_B haemocanis                0.47
17 G001645765 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_B  s__Eperythrozoon_B haemobos                  0.44
18 G902712995 d__Bacteria p__Firmicutes     c__Bacilli             o__Mycoplasmatales    f__Mycoplasmoidaceae    g__Eperythrozoon_B  s__Eperythrozoon_B haemohominis              0.51
```

These are all species with a LOT of gaps, highlighting the issues with these kinds of taxa!


### Correlation between branch lengths and gaps

Let's check to see the scale of the problem. 

First we get the terminal branch lengths, then add them to the main tibble

```{r}
t <- read.tree("r207_original_clean.tree")
tip_names <- t$tip.label
terminal_branch_lengths <- t$edge.length[t$edge[, 2] <= length(tip_names)]
tips <- tibble(id = tip_names, branch_length = terminal_branch_lengths)

# add to main tibble
tax <- left_join(tax, tips, by = "id")
```

Now we can plot it out...

```{r}
ggplot(tax, aes(x = proportion_gaps, y = branch_length)) +
    geom_boxplot(aes(group = proportion_gaps), size = 0.1, alpha = 0.1) +
    geom_point(data = filter(tax, id %in% treeshrink), 1, colour = 'red')
```

![](bl_plot.png)

We can see that the treeshrink taxa all have a lot of gaps, and a mild but not alarming tendency for more gaps to be associated with longer branch lengths.

## Subset taxa


### Subset by gaps

First we'll remove taxa based on our gap threshold

```{r}
subset = tax %>% 
            filter(proportion_gaps<0.20)
```

### Subset by treeshrink

Now we can remove the treeshrink taxa

```{r}
subset = subset %>% 
            filter(!id %in% treeshrink)
```

Actually this didn't remove any more taxa - they were already removed based on gaps, as expected from the plot above.

Now let's see what we've got:

```{r}
unique_counts <- subset %>%
    gather(column, value) %>%                # Convert data to long format
    group_by(column) %>%                     # Group by column names
    summarise(unique_entries = n_distinct(value)) %>% # Count unique values
    arrange(unique_entries)
```

```
> unique_counts
# A tibble: 9 × 2
  column          unique_entries
  <chr>                    <int>
1 clade                        1
2 proportion_gaps             19
3 phylum                     140
4 class                      343
5 order                     1090
6 family                    2570
7 genus                    10603
8 id                       42385
9 species                  42385
```


### One random genome per phylum

Let's start by making three Q matrices, each by selecting a single random ID from every phylum:

```{r}
sample_phylum_id <- function(data) {
  data %>%
    group_by(phylum) %>%
    sample_n(1) %>%
    ungroup() %>%
    pull(id)
}

phylum_1 <- sample_phylum_id(subset)
phylum_2 <- sample_phylum_id(subset)
phylum_3 <- sample_phylum_id(subset)

```

Then lets see how different these are:

```{r}
common_1_2 <- length(intersect(phylum_1, phylum_2))
common_1_3 <- length(intersect(phylum_1, phylum_3))
common_2_3 <- length(intersect(phylum_2, phylum_3))
```

They have 38, 41, and 44 IDs in common, respectively, so the majority (about 100) of the IDs are different between each pair. This is good. If we get very similar matrices with these three lists, that will indicate that the details of the taxon selection didn't matter much.

Let's write out these lists for future usage...

```
writeLines(phylum_1, "phylum_1.txt")
writeLines(phylum_2, "phylum_2.txt")
writeLines(phylum_3, "phylum_3.txt")
```

### One random genome per class

We can also do one random genome per class, like so:


```{r}
sample_class_id <- function(data) {
  data %>%
    group_by(class) %>%
    sample_n(1) %>%
    ungroup() %>%
    pull(id)
}

class_1 <- sample_class_id(subset)
class_2 <- sample_class_id(subset)
class_3 <- sample_class_id(subset)

```

Then lets see how different these are:

```{r}
common_1_2 <- length(intersect(class_1, class_2))
common_1_3 <- length(intersect(class_1, class_3))
common_2_3 <- length(intersect(class_2, class_3))
```

They have ~120 IDs in common, respectively, so the majority (about 200) of the IDs are different between each pair. This is good. 

Let's write out these lists for future usage...

```
writeLines(class_1, "class_1.txt")
writeLines(class_2, "class_2.txt")
writeLines(class_3, "class_3.txt")
```



### One random genome per order

We can also do one random genome per order, like so:


```{r}
sample_order_id <- function(data) {
  data %>%
    group_by(order) %>%
    sample_n(1) %>%
    ungroup() %>%
    pull(id)
}

order_1 <- sample_order_id(subset)
order_2 <- sample_order_id(subset)
order_3 <- sample_order_id(subset)

```

Then lets see how different these are:

```{r}
common_1_2 <- length(intersect(order_1, order_2))
common_1_3 <- length(intersect(order_1, order_3))
common_2_3 <- length(intersect(order_2, order_3))
```

They have ~120 IDs in common, respectively, so the majority (about 200) of the IDs are different between each pair. This is good. 

Let's write out these lists for future usage...

```
writeLines(order_1, "order_1.txt")
writeLines(order_2, "order_2.txt")
writeLines(order_3, "order_3.txt")
```


## Estimate Q Matrices

Now we estimate a model for each of the taxon lists.

I'll walk through the first one in great detail. The rest are just bash scripts based on the first one!

### Q.phylum_1

#### 1. Make a folder and cd to it

```{bash}
mkdir phylum_1
cd phylum_1
cp ../phylum_1.txt .
```

#### 2. Get the subtree

```{bash}
../scripts/get_subtree.py r207_original_clean.tree phylum_1.txt 
```

#### 3. just get the taxa we want from the loci

```{bash}
mkdir -p loci
for loc in ../alignments/*.faa; do
    filename=$(basename $loc)
    faSomeRecords $loc phylum_1.txt loci/${filename}
done

# remove the full alignment, we don't want that!
rm loci/gtdb_r207_bac120_full.faa 
```

#### 4. Split the loci between training and testing

Here we choose at random 20 loci for testsing, which leaves 100 for training.

```{bash}
cd loci
mkdir -p training_loci
mkdir -p testing_loci

test_set=$(ls | sort -R | tail -20)

mv $test_set testing_loci
mv *.faa training_loci
cd ..
```

#### 5. Estimate the models

Now we look through all models and estimate the best model for each locus, using the original r207 sub-tree as our tree.

Here we should be careful to set the number of threads equal to (or less than, if necessary) the number of training loci.

This is just a way to keep IQ-TREE efficient (seems silly, and it is, we're working on it!!!).

The following command broken down:

* `-T 100` 100 threads (I have 100 loci, this will allocate 1 thread to each locus)
* `-p loci/training_loci` this points to my folder of 100 training loci, each in its own `.faa` file
* `-m MFP` this calls modelfinder on every locus
* `-cmax 8` this allows up to 8 free-rate categories per locus (higher is better, but takes longer)
* `-te phylum_1.tree` sets the tree to a pre-estimated tree for these loci (the same tree across all loci)
* `-pre 02_fullcon/iteration_1` sets the output directory

```{bash}
mkdir 02_fullcon # first make the output directory
iqtree2 -T 100 -p loci/training_loci -m MFP -cmax 8 -te phylum_1.tree -pre 02_fullcon/iteration_1
```

The most important output from this is the file `iteration_1.best_scheme.nex`, which has the best model for each locus in it. 

If you scroll through it, you'll see them like this:

```
  charpartition mymodels =
    LG+F+I+R7: gtdb_r207_bac120_PF00466.21.faa,
    Q.pfam+R6: gtdb_r207_bac120_PF02576.18.faa,
    LG+F+I+R8: gtdb_r207_bac120_TIGR00006.faa,
    LG+F+I+R8: gtdb_r207_bac120_TIGR00019.faa,
    LG+I+R8: gtdb_r207_bac120_TIGR00020.faa,
    rtREV+F+I+R6: gtdb_r207_bac120_TIGR00029.faa,
    LG+F+I+R8: gtdb_r207_bac120_TIGR00054.faa,
    Q.pfam+R8: gtdb_r207_bac120_TIGR00059.faa,
    LG+R6: gtdb_r207_bac120_TIGR00061.faa,
...
```

Let's check the models from that analysis:

We can do this with grep and awk...

```{bash}
grep '^ *[^ ]\+:' 02_fullcon/iteration_1.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr
```

This gives us:

```
     61 LG
     24 Q.pfam
      7 Q.yeast
      6 Q.insect
      1 WAG
      1 rtREV
```

Importantly, this tells us that we should use the LG model as the initial model for our analysis, since it's the best fit in 61% of the loci. This means that starting with LG model parameters is likely to be our best bet at getting an even better model.

#### 7. Estimate the Q matrix

Now we estimate the first iteration fo the matrix

```{bash}
iqtree2 -T 100 -S loci/training_loci -p 02_fullcon/iteration_1.best_scheme.nex -te 02_fullcon/iteration_1.treefile --init-model LG --model-joint GTR20+FO -pre 02_fullcon/iteration_1.GTR20
```

Once this is done we can do more iterations...



### Q.class_1

Let's do the same for the class_1.txt list, but this time with a single bash script...


```{bash}

analysis="class_1"

# 1. set up
mkdir $analysis
cd $analysis

echo "Setting up analysis for "$analysis > log.txt

cp ../$analysis.txt .

# 2. get the subtree (produces $analysis.tree)

echo "Getting subtree for "$analysis".txt taxon list" >> log.txt

../scripts/get_subtree.py ../r207_original_clean.tree $analysis.txt 


# 3. just get the taxa we want from the loci

echo "Subsetting alignments" >> log.txt

mkdir -p loci
for loc in ../alignments/*.faa; do
    filename=$(basename $loc)
    faSomeRecords $loc $analysis.txt loci/${filename}
done

# remove the full alignment, we don't want that!
rm loci/gtdb_r207_bac120_full.faa 


# 4. Split the loci between training and testing


echo "splitting alignments into testing and training" >> log.txt

cd loci
test_set=$(ls | sort -R | tail -20)

echo "Test alignments: " >> log.txt
echo $test_set >> log.txt

mkdir -p training_loci
mkdir -p testing_loci


mv $test_set testing_loci
mv *.faa training_loci
cd ..

# check! 

echo "Number of training loci : " >> log.txt
ls loci/training_loci/ | wc -l >> log.txt

echo "Number of testing loci : " >> log.txt
ls loci/testing_loci/ | wc -l >> log.txt

# 5. Estimate the models

echo "Estimating initial models with IQ-TREE2" >> log.txt

mkdir 02_fullcon # first make the output directory
iqtree2 -T 100 -p loci/training_loci -m MFP -cmax 8 -te $analysis.tree -pre 02_fullcon/iteration_1

# get the list of models, and save it to models.txt
grep '^ *[^ ]\+:' 02_fullcon/iteration_1.best_scheme.nex | awk -F: '{print $1}' | awk '{print $NF}' | cut -d'+' -f1 | sort | uniq -c | sort -nr > 02_fullcon/models.txt

echo "List of models best fit to training loci: " >> log.txt
cat 02_fullcon/models.txt >> log.txt

# now we get the init model as the first model in that list

initial_model=$(awk 'NR==1 {print $2}' 02_fullcon/models.txt)
echo "Initial Model will be set to" >> log.txt
echo $initial_model >> log.txt

# 7. Estimate the Q matrix

echo "Estimating Q matrix with IQ-TREE2" >> log.txt

iqtree2 -T 100 -S loci/training_loci -p 02_fullcon/iteration_1.best_scheme.nex -te 02_fullcon/iteration_1.treefile --init-model $initial_model --model-joint GTR20+FO -pre 02_fullcon/iteration_1.GTR20

```
