# Q.GTDB  

The aim for this matrix is to estimate a single amino acid matrix for the GTDB dataset.

The dataset comprises 120 genes, each for ~85K taxa. The main challenges are that this is just too many taxa. So we need to come up with sensible ways of estimating matrices from a subset of the taxa. I will use a lot of the suggestions from https://github.com/roblanf/Q.GTDB, with some modifications.


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

## Subset taxa

### Remove taxa with lots of gaps

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
ggplot(gaps, aes(x = proportion_gaps)) + geom_histogram(bins=100)
```

Let's look at the proportion of sequences with fewer gaps than each threshold:

```{r}
results <- map_dfr(seq(0, 1, by = 0.01), function(threshold) {
  proportion_below_threshold <- mean(gaps$proportion_gaps < threshold)
  tibble(threshold = threshold, proportion_below = proportion_below_threshold)
})

head(results, 20)
```

This gives:

```

```

So a sensible cutoff here is to only keep taxa with <15% gaps - this means we only throw out x% of the dataset.

### Subset the remaining taxa

There are two obvious ways to subset the taxa - by taxonomy or by phylogenetic diversity.

Let's create a few lists of different subsets first.

First let's 
