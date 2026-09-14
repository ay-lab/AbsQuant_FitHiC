#!/bin/bash

#SBATCH --job-name=intersect_per_replicate
#SBATCH --output=1.4_intersect_replicates_%j.out
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64g

# Intersects the concatenated loop set from 1.3 with a list of fithic
# significances files -- one per biological replicate. Like 1.3, this runs the
# work itself rather than generating jobs.
#
# The question it answers: of the high-confidence loops called on the combined
# sample, which are also called by fithic in each individual replicate?
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (workingDir defaults to $PWD)
#
#   sbatch 1.4_intersect_replicates.sh (preferable; bash loads a full fithic
#                                       table per replicate, GBs of RAM)
#   bash   1.4_intersect_replicates.sh
#
# Override any INPUT VARIABLE below on the command line:
#
#   FITHIC_FILES="/path/rep1.txt /path/rep2.txt.gz" bash 1.4_intersect_replicates.sh
#
#   sbatch --export=ALL,FITHIC_FILES="/path/rep1.txt /path/rep2.txt.gz" \
#     1.4_intersect_replicates.sh
# ---------------------------------------------------------------------------
#
# Input  (from 1.3): <resultsRoot>/filtered_loops_<resKb>kb_fdr<fdr>/<condition>/
#                        <condition>.coords.fdr<fdr>.txt
#        (fithic):   one significances.txt[.gz] per replicate
#
# Output: <resultsRoot>/fithic_filtered_loops_bioreplicates_<resKb>kb_fdr<fdr>/
#             loops/<condition>/<replicate>.<...>.fdr<fdr>.significances.txt
#
source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================
# THIS MUST MATCH WHAT 1.3 WAS RUN ON.
# Conditions to process.
read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}"

# Which fithic files to intersect against, in order of precedence:
#   1. FITHIC_FILES - a space-separated list, used for every condition
#   2. the fithicFiles array below, if you fill it in
#   3. otherwise derived from fithicDir + replicateNamesList + fithicTemplate
fithicFiles=()
if [ -n "${FITHIC_FILES:-}" ]; then
    read -r -a fithicFiles <<< "${FITHIC_FILES}"
fi

# Biological replicates, used only when deriving paths (case 3 above).
read -r -a replicateNamesList <<< "${REPLICATE_NAMES:-pTh17-1 npTh17 Treg Th1 Th2 Th0}"

# Root of per-replicate fithic output
fithicDir="${FITHIC_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"

# Filename of the fithic call table.
fithicTemplate="${FITHIC_TEMPLATE:-L20000.U3000000.p2.b200.spline_pass2.res${RESOLUTION:-10000}.significances.txt}"

# Output root, shared by steps 1.1 - 1.4.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"

# This repo, holding the .py files.
workingDir="${WORKING_DIR:-$(pwd)}"

resolution="${RESOLUTION:-10000}"
fdrThreshold="${FDR_THRESHOLD:-0.01}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
curr_date=$(date +"%y%m%d")
resKb=$((resolution / 1000))
inputDir="${resultsRoot}/filtered_loops_${resKb}kb_fdr${fdrThreshold}"
outputDir="${resultsRoot}/fithic_filtered_loops_bioreplicates_${resKb}kb_fdr${fdrThreshold}/loops"
pythonFile="${workingDir}/1.4_intersect_replicates.py"

mkdir -p ${outputDir}

# Check if Python script exists
if [ ! -f "${pythonFile}" ]; then
    echo "Error: Python script not found: ${pythonFile}"
    exit 1
fi

if [ ! -d "${inputDir}" ]; then
    echo "Error: 1.3 output directory not found: ${inputDir}"
    echo "  Run 1.3_concat_loops.sh first, or set RESULTS_DIR= to where it wrote."
    exit 1
fi

echo "=================================================="
echo "Intersecting loops with per-replicate fithic calls"
echo "=================================================="
echo "Input directory:  ${inputDir}"
echo "Output directory: ${outputDir}"
echo ""

for subset in ${subsets[@]}; do
    inputFile="${inputDir}/${subset}/${subset}.coords.fdr${fdrThreshold}.txt"

    if [ ! -f "${inputFile}" ]; then
        echo "Warning: No concatenated loop file for ${subset}, skipping..."
        echo "  Expected: ${inputFile}"
        continue
    fi

    # Build the list of fithic files for this subset.
    subsetFithicFiles=()
    if [ ${#fithicFiles[@]} -gt 0 ]; then
        subsetFithicFiles=("${fithicFiles[@]}")
    else
        for replicateName in ${replicateNamesList[@]}; do
            f="${fithicDir}/${replicateName}/fithic/${resolution}/${replicateName}.${fithicTemplate}"
            if [ -f "${f}" ]; then
                subsetFithicFiles+=("${f}")
            elif [ -f "${f}.gz" ]; then
                # Read the .gz directly; the python side decompresses on the fly.
                subsetFithicFiles+=("${f}.gz")
            else
                echo "  Warning: no fithic file for ${replicateName} (${f}[.gz])"
            fi
        done
    fi

    if [ ${#subsetFithicFiles[@]} -eq 0 ]; then
        echo "Warning: No fithic files found for ${subset}, skipping..."
        continue
    fi

    outDir="${outputDir}/${subset}"
    mkdir -p ${outDir}

    echo "Processing: ${subset}"
    echo "  Loop file: ${inputFile}"
    echo "  Fithic files: ${#subsetFithicFiles[@]}"
    echo "  Output directory: ${outDir}"

    python3 ${pythonFile} \
        --combined_loops_file ${inputFile} \
        --replicate_files ${subsetFithicFiles[@]} \
        --output_dir ${outDir} \
        --fdr_threshold ${fdrThreshold} \
        --verbose

    echo ""
done

echo "=================================================="
echo "Intersected loop tables in: ${outputDir}"
echo "=================================================="
