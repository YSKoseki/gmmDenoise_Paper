#!/bin/bash

# 01_SRAfastq_DRA009149_v2.sh
# Last updated on 2022.11.13 by YK
# A shell script to download FASTQs from NCBI SRA

# Set variables: directory path and SRA BioProject ID
DIR="./01-SRAfastq_DRA009149"
BIO=PRJDB8911

# Make the output and temp directory
mkdir -pv "${DIR}/FASTQ"
mkdir -pv "${DIR}/FASTQ/temp"

# Use SRA Run Selector to download the accession list as "SRR_Acc_List.txt"
# Make sure that the list is in the current directory
 
# Download SRAs
prefetch \
  --output-directory "${DIR}/FASTQ/temp" \
  --option-file SRR_Acc_list_test.txt

# Convert the SRAs into FASTQs
for file in ${DIR}/FASTQ/temp/*/*.sra
do
  fasterq-dump \
    --split-files \
    --outdir "${DIR}/FASTQ" \
    --progress \
    --threads 12 \
    ${file}
done

# Compress the FASTQs
for file in ${DIR}/FASTQ/*.fastq
do
  gzip \
    --keep \
    --verbose \
    --stdout \
    ${file} > ${file}.gz
done

# Clean-up the temp directory
rm -r "${DIR}/FASTQ/temp"
