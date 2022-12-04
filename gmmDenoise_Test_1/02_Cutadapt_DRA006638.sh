#!/bin/bash

# 02_Cutadapt_DRA006638.sh
# Last updated on 2022.9.27 by YK
# A shell script to trim the pair of primers from paired-end sequencing reads

# Set variables: directory paths, FASTQ suffixes, and primer sequences
INDIR="./01-SRAfastq_DRA006638/FASTQ"
OUTDIR="./02-Cutadapt_DRA006638"
FWSUF=_1.fastq.gz
RVSUF=_2.fastq.gz
FWSUF2=_trimmed_R1.fastq.gz # for trimmed forward-FASTQs
RVSUF2=_trimmed_R2.fastq.gz # for trimmed reverse-FASTQs
FWADAPT=^NNNNNNCCGGTTGCATATATGGACCTATTACNNN
RVADAPT=^NNNNNNGCTATTRTAGTCTGGTAACGCAAG

# Make output directory and create files of command log and stat table
mkdir -pv "${OUTDIR}/trimmedFASTQ" 
: > "${OUTDIR}/primtrim.txt"
: > "${OUTDIR}/primtrim.tsv"

# Write header (tab-separated elements) to the stat table
echo -e "sample\ttotRP\tfiltdRP\tpctRP\ttotBP\tfiltdBP\tpctBP" > "${OUTDIR}/primtrim.tsv"

# Find forward-read FASTQs, return their file names, and sort them 
FASTQR1=`find ${INDIR} -type f -name "*${FWSUF}" | sed 's!^.*/!!' | sort`

# Primer trimming by cutadapt
echo "${FASTQR1}" | while read line
do
  FW=`echo ${line}`
  RV=`echo ${FW} | sed "s/${FWSUF}/${RVSUF}/g"`
  FW2=`echo ${FW} | sed "s/${FWSUF}/${FWSUF2}/g"`
  RV2=`echo ${RV} | sed "s/${RVSUF}/${RVSUF2}/g"`
  RUN=`echo ${line} | cut -f 1 -d "_"`
  # Run 'cutadapt'
  OUTPUT=`cutadapt --discard-untrimmed \
    -g ${FWADAPT} -G ${RVADAPT} \
    -o "${OUTDIR}/trimmedFASTQ/${FW2}" -p "${OUTDIR}/trimmedFASTQ/${RV2}" \
    "${INDIR}/${FW}" "${INDIR}/${RV}"`
  # Save logs
  echo -e "${RUN}\n\n${OUTPUT}\n\n\n" >> "${OUTDIR}/primtrim.txt" 2>&1
  # Extract the selected stats and write them to the stat table
  if echo "${OUTPUT}" | grep 'No reads processed!' >/dev/null; then
    echo -e "${RUN}\t\t\t\t\t\t" >> "${OUTDIR}/${ID}.tsv" 2>&1
  else
    TOTRP=`echo "${OUTPUT}" | grep 'Total read pairs' | sed -e 's/[^0-9]//g'`
    FILTDRP=`echo "${OUTPUT}" | grep 'Pairs written' | cut -f 2 -d '(' | sed -e 's/[^0-9]//g'`
    PCTRP=`echo "${OUTPUT}" | grep 'Pairs written' | cut -f 3 -d '(' | tr -d ')'`
    TOTBP=`echo "${OUTPUT}" | grep 'Total basepairs' | sed -e 's/[^0-9]//g'`
    FILTDBP=`echo "${OUTPUT}" | grep 'Total written' | cut -f 2 -d '(' | sed -e 's/[^0-9]//g'`
    PCTBP=`echo "${OUTPUT}" | grep 'Total written' | cut -f 3 -d '(' | tr -d ')'`
    echo -e "${RUN}\t${TOTRP}\t${FILTDRP}\t${PCTRP}\t${TOTBP}\t${FILTDBP}\t${PCTBP}" >> "${OUTDIR}/primtrim.tsv" 2>&1
  fi
done
