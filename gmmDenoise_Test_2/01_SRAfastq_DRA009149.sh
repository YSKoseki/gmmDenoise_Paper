#!/bin/bash

# 01_SRAfastq_DRA009149.sh
# Last updated on 2022.11.9 by YK
# A shell script to download FASTQs from NCBI SRA

# Set variables: directory path and SRA BioProject ID
DIR="./01-SRAfastq_DRA009149"
BIO=PRJDB8911

# Make the output directory
mkdir -pv "${DIR}/FASTQ"

# Download multiple SRAs and convert them into gzipped FASTQs
grabseqs sra \
  ${BIO} \
  -t 6 \
  -m "../metadata.csv" \
  -o "${DIR}/fastq"
