#!/bin/bash

#SBATCH --job-name=concat_loops_all_chr
#SBATCH --output=1.3_concat_loops_%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g

# Concatenates the per-chromosome output of 1.2 into one loop table per
# condition, renaming the columns to fithic's so 1.4 can join on coordinates.
# Unlike 1.1 and 1.2 this does the work itself rather than generating SLURM jobs.
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (workingDir defaults to $PWD)
#
#   sbatch 1.3_concat_loops.sh (preferable)
#   bash   1.3_concat_loops.sh
#
# Override any INPUT VARIABLE below on the command line:
#
#   SUBSETS="condA condB" FORCE=1 bash 1.3_concat_loops.sh
#
#   sbatch --export=ALL,SUBSETS="condA condB",FORCE=1 1.3_concat_loops.sh
# ---------------------------------------------------------------------------
#
# Input  (from 1.2): <resultsRoot>/filtered_loops_per_chr_fdr<fdr>/<condition>/
#                        filtered_loops.fdr<fdr>.<chrom>.txt
#                    columns: chr  left  right  size
#
# Output (for 1.4): <resultsRoot>/filtered_loops_<resKb>kb_fdr<fdr>/<condition>/
#                        <condition>.coords.fdr<fdr>.txt
#                    columns: chr1  fragmentMid1  chr2  fragmentMid2  size

source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================
# THIS MUST MATCH WHAT 1.2 WAS RUN ON.

# Conditions to concatenate. 
read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}"

# Chromosomes to concatenate, in output order. Mouse autosomes by default; use
# chr1..chr22 for human.
read -r -a chroms <<< "${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19}"

# Output root, shared by steps 1.1 - 1.4.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"

# This repo. Defaults to the directory you run from.
workingDir="${WORKING_DIR:-$(pwd)}"

resolution="${RESOLUTION:-10000}"
fdrThreshold="${FDR_THRESHOLD:-0.01}"

# Set FORCE=1 to overwrite an existing output file instead of skipping it.
FORCE="${FORCE:-0}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
curr_date=$(date +"%y%m%d")
resKb=$((resolution / 1000))
inputDir="${resultsRoot}/filtered_loops_per_chr_fdr${fdrThreshold}"
outputDir="${resultsRoot}/filtered_loops_${resKb}kb_fdr${fdrThreshold}"

mkdir -p ${outputDir}

if [ ! -d "${inputDir}" ]; then
    echo "Error: 1.2 output directory not found: ${inputDir}"
    echo "  Run 1.2_filter_loops.sh first, or set RESULTS_DIR= to where it wrote."
    exit 1
fi

echo "=================================================="
echo "Concatenating filtered loops across chromosomes"
echo "=================================================="
echo "Input directory:  ${inputDir}"
echo "Output directory: ${outputDir}"
echo "Chromosomes:      ${#chroms[@]} (${chroms[0]}..${chroms[$((${#chroms[@]}-1))]})"
echo ""

for subset in ${subsets[@]}; do
    dir="${inputDir}/${subset}"
    outSubDir="${outputDir}/${subset}"
    outFile="${outSubDir}/${subset}.coords.fdr${fdrThreshold}.txt"

    if [ ! -d "${dir}" ]; then
        echo "Warning: No 1.2 output for ${subset}, skipping..."
        continue
    fi

    if [ -f "${outFile}" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "Skipping ${subset}: output already exists (FORCE=1 to overwrite)"
        echo "  ${outFile}"
        echo ""
        continue
    fi

    mkdir -p ${outSubDir}
    echo "Processing: ${subset}"

    # Write to a temp file so an interrupted run cannot leave a half-written
    # table that the FORCE check above would then skip over.
    tmpFile="${outFile}.tmp.$$"
    printf 'chr1\tfragmentMid1\tchr2\tfragmentMid2\tsize\n' > ${tmpFile}

    nFound=0
    nMissing=0
    for chrom in ${chroms[@]}; do
        chromFile="${dir}/filtered_loops.fdr${fdrThreshold}.${chrom}.txt"
        if [ ! -f "${chromFile}" ]; then
            echo "  Warning: missing ${chrom} (${chromFile})"
            nMissing=$((nMissing + 1))
            continue
        fi
        # 1.2 writes: chr  left  right  size  (with a header line).
        # fithic wants both ends named, so chr is emitted twice: loops here are
        # cis by construction, since 1.2 filters chr1 == chr2 == --chromosome.
        tail -n +2 ${chromFile} | awk 'NF{print $1"\t"$2"\t"$1"\t"$3"\t"$4}' >> ${tmpFile}
        nFound=$((nFound + 1))
    done

    if [ ${nFound} -eq 0 ]; then
        echo "  Error: no chromosome files found for ${subset}, nothing written"
        rm -f ${tmpFile}
        echo ""
        continue
    fi

    mv ${tmpFile} ${outFile}
    nLoops=$(( $(wc -l < ${outFile}) - 1 ))
    echo "  ✓ ${nFound} chromosomes, ${nMissing} missing, ${nLoops} loops"
    echo "  ✓ Wrote: ${outFile}"
    echo ""
done

echo "=================================================="
echo "Concatenated tables in: ${outputDir}"
echo "Next step:"
echo "  bash ${workingDir}/1.4_intersect_replicates.sh"
echo "=================================================="
