#!/usr/bin/env python3

import sys
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np

# Read newick
# tree_file = "/home/fredjaya/GitHub/gtdb_trees/data/trees/gtdb_r207_bac120_unscaled.decorated.tree"
tree_file = sys.argv[1]
tree = Phylo.read(tree_file, "newick")

# Output total branch length
print(tree.total_branch_length())

# Save all branch lengths to list and plot histogram
bl = [clade.branch_length for clade in tree.find_clades() if clade.branch_length]

log_bin_edges = np.logspace(np.log10(min(bl)), np.log10(max(bl)), num=50)

plt.hist(bl, bins=log_bin_edges)
plt.xscale("log")
plt.xlabel("Branch length")
plt.ylabel("Count")
text = f"Total number of branches: {len(bl)}"
plt.text(0.7, 0.9, f"{len(bl)} branches", transform=plt.gca().transAxes, fontsize=10)

plt.savefig("branch_length_histogram.png")

# plt.show()
