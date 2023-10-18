#!/bin/bash

# Path to directory with all loci
DIR=$1

mkdir -p training_loci
mkdir -p testing_loci

test_set=$(ls $DIR | sort -R | tail -20)

mv $test_set testing_loci
mv ${DIR}/*.faa training_loci
