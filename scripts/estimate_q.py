#!/usr/bin/env python3

import argparse
import subprocess
import os
import re
import sys
import math
from pearsons import read_qmatrix, pearsons_correlation


def run_command(cmd):
    """
    Wrapper for bash commands to deal with stderr and exit_codes
    """
    if args.verbose:
        print(cmd, "\n")
    process = subprocess.Popen(
        cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    exit_code = process.returncode
    if exit_code > 0:
        sys.exit(stderr)
    return stdout, stderr, exit_code


def grep_iqtree(string, filename):
    index = [idx for idx, l in enumerate(lines) if string in l][0]
    with open(filename, "w") as f:
        # TODO: softcode for single concat tree
        f.writelines([l.lstrip() for l in lines[index + 2 : index + 22]])


def read_qmatrix(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    # Convert lines to floats
    lines = [list(map(float, line.strip().split())) for line in lines]

    n = len(lines)
    matrix = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j, value in enumerate(lines[i]):
            matrix[i][j] = value
            matrix[j][i] = value
    return matrix


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", type=str, help="Pick one: uncon, semicon, fullcon", required=True
    )
    parser.add_argument(
        "-s", "--loci", type=str, help="Path to concatenated loci", required=True
    )
    parser.add_argument(
        "-mset",
        "--starting_models",
        type=str,
        help='e.g. "LG,Q.yeast,Q.inset"',
        required=True,
    )
    parser.add_argument(
        "-te",
        "--constraint_tree",
        help="Estimate Q with a fixed topology.",
        required=False,
    )
    parser.add_argument(
        "-T", "--threads", help="Threads to use (-T)", default=int(1), required=False
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print commands", required=False
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = cli()
    if args.mode not in ["uncon", "semicon", "fullcon"]:
        sys.exit("--mode must be one of uncon, semicon, or fullcon")
    iqtree_binary = "~/Downloads/iqtree-2.2.2.7-Linux/bin/iqtree2"
    seed = 1
    pearsons_rho = 0
    rho_cutoff = 0.999
    iteration = 0  # pro tip: change this to resume from i
    models = args.starting_models

    while pearsons_rho < rho_cutoff:
        iteration += 1
        if iteration > 1:
            models = (
                args.starting_models + f",02_{args.mode}/Q.{args.mode}_i{iteration-1}"
            )

        run_command(f"mkdir -p 02_{args.mode}")

        if args.mode == "uncon":
            run_command(
                f"{iqtree_binary} --seed {seed} -T {args.threads} -S {args.loci} -m MFP -mset {models} -cmax 8 -pre 02_uncon/i{iteration}"
            )
            run_command(
                f"{iqtree_binary} -seed {seed} -T {args.threads} -S {args.loci} -te 02_{args.mode}/i1.treefile --init-model LG --model-joint GTR20+FO -pre 02_{args.mode}/i{iteration}.GTR20"
            )
        elif args.mode == "fullcon":
            run_command(
                f"{iqtree_binary} --seed {seed} -T {args.threads} -p {args.loci} -m MFP -mset {models} -cmax 8 -te {args.constraint_tree} -pre 02_fullcon/i{iteration}"
            )
        elif args.mode == "semicon":
            run_command(
                f"{iqtree_binary} --seed {seed} -T {args.threads} -p {args.loci} -m MFP -mset {models} -cmax 8 -pre 02_semicon/i{iteration}"
            )

        if args.mode in ["semicon", "fullcon"]:
            run_command(f"sed -i 's/, //' 02_{args.mode}/i{iteration}.best_scheme.nex")
            run_command(
                f"{iqtree_binary} -seed {seed} -T {args.threads} -S {args.loci} -p 02_{args.mode}/i{iteration}.best_scheme.nex -te 02_{args.mode}/i1.treefile --init-model LG --model-joint GTR20+FO -pre 02_{args.mode}/i{iteration}.GTR20"
            )

        with open(f"02_{args.mode}/i{iteration}.GTR20.iqtree") as iq:
            """
            Parse .iqtree file to get new Q.matrix and calculate the
            correlation between the previous one
            """
            lines = iq.readlines()
            # Write Q-matrix
            grep_iqtree(
                "can be used as input for IQ-TREE",
                f"02_{args.mode}/Q.{args.mode}_i{iteration}",
            )

        if iteration > 1:
            qold = read_qmatrix(f"02_{args.mode}/Q.{args.mode}_i{iteration-1}")
            qnew = read_qmatrix(f"02_{args.mode}/Q.{args.mode}_i{iteration}")
            rho = pearsons_correlation(qold, qnew)
            print(f"rho={rho}")
            if rho > rho_cutoff:
                print("Very similar to previous Q, stop.")
            pearsons_rho = rho
