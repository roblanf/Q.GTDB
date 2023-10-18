#!/usr/bin/env python3

import re
import sys
import math


def get_top_models(scheme_file):
    """
    Input: .best_scheme file output by iqtree2 -S all_loci -m MF
    Output: top_models.txt
        - Number of best-fitting Q-matrices for all loci
        - Total loci
        - Cut-off (i.e. get models that appear in the top X% of loci)
            default: 0.9
        - Comma-delimited string of inferred starting models
    """
    schemes = []
    model_counts = {}

    with open(scheme_file, "r") as scheme_file:
        # Parse best_scheme file
        for line in scheme_file:
            line = line.strip()
            line = re.sub(" =", "", line)
            model, loc_name = line.split(", ")
            # Remove all appended "+F"
            model = re.sub("F$", "", model)
            schemes.append([line, loc_name])

            # Count frequency of models
            model_counts[model] = model_counts.get(model, 0) + 1

    model_counts = sorted(model_counts.items(), key=lambda x: x[1], reverse=True)
    print(model_counts)

    """
    Identify the models that appear in the top X% of loci.
    Default = 90%
    """
    nloci = len(schemes)
    cutoff = 0.9
    maxloci = math.ceil(nloci * cutoff)

    print(
        f"Total loci: {nloci}\nCut-off: {cutoff}\nSelecting most frequent models up to {maxloci} loci."
    )

    cumulative_count = 0
    starting_models = []

    for i in model_counts:
        if cumulative_count <= maxloci:
            cumulative_count += i[1]
            starting_models.append(i[0])

    print(f"Starting models: {','.join(starting_models)}")
    return


if __name__ == "__main__":
    best_scheme_file = sys.argv[1]
    get_top_models(best_scheme_file)
