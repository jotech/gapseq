#!/usr/bin/env bash
set -euo pipefail

INPUT="${1:-proteins.fasta}"
PREFIX="${2:-protein_family}"
THREADS="${3:-1}"

echo "Input:   $INPUT"
echo "Prefix:  $PREFIX"
echo "Threads: $THREADS"

# ------------------------------------------------------------
# 1. Confirm that sequence identifiers are unique
# ------------------------------------------------------------
seqkit seq -n "$INPUT" |
    sed 's/[[:space:]].*$//' |
    sort |
    uniq -d > "${PREFIX}.duplicate_ids.txt"

if [[ -s "${PREFIX}.duplicate_ids.txt" ]]; then
    echo "ERROR: Duplicate sequence identifiers were found:"
    head -20 "${PREFIX}.duplicate_ids.txt"
    echo "Rename duplicate identifiers before running the pipeline."
    exit 1
fi

# ------------------------------------------------------------
# 2. Multiple sequence alignment
# ------------------------------------------------------------
mafft \
    --auto \
    --thread "$THREADS" \
    "$INPUT" \
    > "${PREFIX}.aligned.fasta"

# ------------------------------------------------------------
# 3. Trim poorly aligned columns
# ------------------------------------------------------------
clipkit \
    "${PREFIX}.aligned.fasta" \
    -m smart-gap \
    -o "${PREFIX}.trimmed.fasta"

# ------------------------------------------------------------
# 4. Infer a maximum-likelihood tree with IQ-TREE 3
#
# -s          : input alignment
# -st AA      : amino-acid alignment
# -m MFP      : select the best-fit protein model with ModelFinder
# -B 1000     : 1,000 ultrafast bootstrap replicates
# --alrt 1000 : 1,000 SH-aLRT replicates
# -T          : number of CPU threads
# --prefix    : prefix for all IQ-TREE output files
# ------------------------------------------------------------
iqtree3 \
    -s "${PREFIX}.trimmed.fasta" \
    -st AA \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$THREADS" \
    --prefix "${PREFIX}.iqtree"

echo
echo "Finished."
echo
echo "Main outputs:"
echo "  Alignment:       ${PREFIX}.aligned.fasta"
echo "  Trimmed MSA:     ${PREFIX}.trimmed.fasta"
echo "  ML tree:         ${PREFIX}.iqtree.treefile"
echo "  Consensus tree:  ${PREFIX}.iqtree.contree"
echo "  IQ-TREE report:  ${PREFIX}.iqtree.iqtree"
echo "  IQ-TREE log:     ${PREFIX}.iqtree.log"


