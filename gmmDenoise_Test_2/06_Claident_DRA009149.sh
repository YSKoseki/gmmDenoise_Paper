#!/bin/bash

# 06_Claident_DRA009149.sh
# Last updated on 2022.11.9 by YK
# A shell script to run the program Claident to assign taxonomy to sequences 
# For the commands and their options, see the document on https://github.com/astanabe/MetabarcodingTextbook
# clidentseq 0.9.2022.01.26

# Set variables
# Number of CPU cores to use
NCPU=6
## Path of input FASTA
INFILE="./05-Merger_DRA009149/02-ASV_seqs_merged.fa"
## Path of output directory
OUTDIR="./06-Claident_DRA009149"

# The wrapper function
  runclaident() {
    local CPU=$1
    local INFILE=$2
    local OUTDIR=$3
    mkdir -p "${OUTDIR}"
#   Construct cache DB to reduce runtime of 'clidentseq' command
    clmakecachedb \
      --blastdb=animals_mt_species \
      --numthreads=${CPU} \
      "${INFILE}" \
      "${OUTDIR}/01-animals_mt_species_cached"
#   The first axonomic assignment based on the default 'QCauto' method
##  Retrieve neighborhood sequences from the cache DB
    clidentseq \
      --blastdb="${OUTDIR}/01-animals_mt_species_cached" \
      --numthreads=${CPU} \
      "${INFILE}" \
      "${OUTDIR}/02-QCauto_GBID.dat"
##  Assign taxonomy to the retrieved sequences, 
##    dealing with inconsistency among matched sequences by the set of arguments:
##    'maxpopposer' the maximum permissible proportion of opposer sequences (%)
##    'minsoratio' the minimum permissible ratio (X:1) of supporter(X) to opposer(1) sequences
    classigntax \
      --taxdb=animals_mt_species \
      --maxpopposer=0.05 \
      --minsoratio=19 \
      "${OUTDIR}/02-QCauto_GBID.dat" \
      "${OUTDIR}/03-QCauto_MaxPOP05_MinSOR19.tsv"
#   The second taxonomic assignment based on a 'Nearest Neighborhood' method
##  Retrieve top-1 (and tie) sequences whose percent-identity to query is 99% or higher
    clidentseq \
      blastn -task megablast -word_size 16 end \
      --method="1,99%" \
      --blastdb="${OUTDIR}/01-animals_mt_species_cached" \
      --numthreads=${CPU} \
      "${INFILE}" \
      "${OUTDIR}/04-NN99_GBID.dat"
##  Assign taxonomy to the retrieved sequences,
##    reducing the required minimum number of neighborhood sequences to 1 by 'minnsupporter' and 
##    dealing with inconsistency among matched sequences by:
##    'maxpopposer' the maximum permissible proportion of opposer sequences (%)
##    'minsoratio' the minimum permissible ratio (X:1) of supporter(X) to opposer(1) sequences
    classigntax \
      --taxdb=animals_mt_species \
      --minnsupporter=1 \
      "${OUTDIR}/04-NN99_GBID.dat" \
      "${OUTDIR}/05-NN99.tsv"
#   Merge the two taxonomic assignments based on strict consensus LCA rule, giving preference to the first  
    clmergeassign \
      --priority=descend \
      "${OUTDIR}/03-QCauto_MaxPOP05_MinSOR19.tsv" \
      "${OUTDIR}/05-NN99.tsv" \
      "${OUTDIR}/06-QCauto_MaxPOP05_MinSOR19+NN99.tsv"
}

# Run the wrapper
runclaident ${NCPU} ${INFILE} ${OUTDIR}
