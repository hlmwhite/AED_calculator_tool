#!/bin/bash

# wrapper tool for running the AED calculator script on the command line

# example usage
# 'AED_wrapper.sh test.genome.gff stringtie.test.gff'

gffcmp_string=$1
ann_gff=$2
ev_gff=$3

    annIDs=$(printf "$gffcmp_string" | awk '{print $4}' | tr ',' '\n')
    evIDs=$(printf "$gffcmp_string" | awk '{print $5}'| tr ',' '\n')
    for annID in $annIDs
    do
        for evID in $evIDs
        do
            echo "Calculating AED for annotation ID: $annID and evidence ID: $evID"
            python3 example_calculator_v2.py -A $ann_gff -E $ev_gff -a $annID -e $evID 
        done
    done