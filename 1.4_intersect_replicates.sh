#!/bin/bash

#SBATCH --job-name=intersect_per_replicate
#SBATCH --output=1.4_intersect_replicates_%j.out
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64g

# Script to intersect the concatenated loop set from 1.3 with a list of fithic
# significances files -- one per biological replicate. Like 1.3, this runs the
# work directly rather than generating SLURM scripts.
#
# The question it answers: of the high-confidence loops called on the combined
# sample, which are also called by fithic in each individual replicate?
#
# Input  (from 1.3): <resultsRoot>/filtered_loops_<resKb>kb_fdr<fdr>/<subset>/
#                        <subset>.coords.fdr<fdr>.txt
#        (fithic):   one significances.txt[.gz] per replicate
#
# Output: <resultsRoot>/fithic_filtered_loops_bioreplicates_<resKb>kb_fdr<fdr>/
#             loops/<subset>/<replicate>.<...>.fdr<fdr>.significances.txt
#
# Which fithic files are used, in order of precedence:
#   1. FITHIC_FILES  - a space-separated list you supply, used for every subset
#   2. the fithicFiles array below, if you fill it in
#   3. otherwise they are derived from fithicDir + replicateNamesList using the
#      standard fithic filename template
#
# .gz inputs are read directly; no decompression step is needed.

source ~/.bashrc

# Parameters
curr_date=$(date +"%y%m%d")
subsets=("pTh17-1" "npTh17" "Treg" "Th1" "Th2" "Th0")
resolution=10000
fdrThreshold=0.01
resKb=$((resolution / 1000))
# Biological replicates to intersect against each subset's loop set. Used to
# derive fithic paths when no explicit list is given.
replicateNamesList=(pTh17-1 npTh17 Treg Th1 Th2 Th0)

# Explicit list of fithic files. Leave empty to derive them (see header).
# Overridden by FITHIC_FILES= on the command line.
fithicFiles=()
if [ -n "${FITHIC_FILES:-}" ]; then
    read -r -a fithicFiles <<< "${FITHIC_FILES}"
fi

# Directories
# workingDir is this repo. Derived from the script's own location so a clone
# works anywhere; override with WORKING_DIR= if you keep the code elsewhere.
workingDir="${WORKING_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# baseDir is the DATA project. It is site-specific -- override with BASE_DIR=
# rather than editing this file.
baseDir="${BASE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay}"
# resultsRoot is where 1.2 and 1.3 wrote their output. Override with RESULTS_DIR=.
resultsRoot="${RESULTS_DIR:-${baseDir}/yard/251027_absLoopQuant/results}"
inputDir="${resultsRoot}/filtered_loops_${resKb}kb_fdr${fdrThreshold}"
outputDir="${resultsRoot}/fithic_filtered_loops_bioreplicates_${resKb}kb_fdr${fdrThreshold}/loops"
# Where per-replicate fithic output lives, when deriving paths rather than
# passing them in. Expected layout:
#   ${fithicDir}/<replicate>/fithic/<resolution>/<replicate>.<template>
fithicDir="${FITHIC_DIR:-${baseDir}/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"
fithicTemplate="L20000.U3000000.p2.b200.spline_pass2.res${resolution}.significances.txt"
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
