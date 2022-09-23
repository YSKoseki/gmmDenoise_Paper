#!/bin/bash

# 01_SRAfastq_DRA005106.sh
# Last updated on 2022.2.27 by YK
# A shell script to download FASTQs from NCBI SRA

# Set variables: directory path and SRA BioProject ID
DIR="./01-SRAfastq_DRA005106"
BIO=PRJDB5158

# Make the output directory
mkdir -pv "${DIR}/FASTQ"

# Download multiple SRAs and convert them into gzipped FASTQs
grabseqs sra \
  ${BIO} \
  -t 10 \
  -m "../metadata.csv" \
  -o "${DIR}/fastq"
