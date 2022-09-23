#!/bin/bash

# 01_SRAfastq_DRA010703.sh
# Last updated on 2022.4.15 by YK
# A shell script to download FASTQs from NCBI SRA

# Set variables: directory path and SRA BioProject ID
DIR="./01-SRAfastq_DRA010703"
BIO=PRJDB10433

# Make the output directory
mkdir -pv "${DIR}/FASTQ"

# Download multiple SRAs and convert them into gzipped FASTQs
grabseqs sra \
  ${BIO} \
  -t 10 \
  -m "../metadata.csv" \
  -o "${DIR}/fastq"
