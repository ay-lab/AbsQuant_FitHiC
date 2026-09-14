#!/bin/bash

#SBATCH --job-name=concat_loops_all_chr
#SBATCH --output=1.3_concat_loops_%j.out
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g

# Script to concatenate the per-chromosome output of 1.2 into one loop table
# per subset. Unlike 1.1 and 1.2 this does the work directly rather than
# generating SLURM scripts -- it is a cat/awk pass over ~20 small files.
#
# Input  (from 1.2): <resultsRoot>/filtered_loops_per_chr_fdr<fdr>/<subset>/
#                        filtered_loops.fdr<fdr>.<chrom>.txt
#                    columns: chr  left  right  size
#
# Output (for 1.4): <resultsRoot>/filtered_loops_<resKb>kb_fdr<fdr>/<subset>/
#                        <subset>.coords.fdr<fdr>.txt
#                    columns: chr1  fragmentMid1  chr2  fragmentMid2  size
#
# The output column names match fithic's, which is what lets 1.4 join the two
# tables on coordinates.
#
# Set FORCE=1 to overwrite an existing output file instead of skipping it.

source ~/.bashrc

# Parameters
curr_date=$(date +"%y%m%d")
subsets=("pTh17-1" "npTh17" "Treg" "Th1" "Th2" "Th0")
resolution=10000
fdrThreshold=0.01
resKb=$((resolution / 1000))
# Chromosomes to concatenate, in output order. This is a mouse (mm10) autosome
# set; for human use chr1..chr22, and add chrX / chrY if 1.2 was run on them.
chroms=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19)

# Directories
# workingDir is this repo. Derived from the script's own location so a clone
# works anywhere; override with WORKING_DIR= if you keep the code elsewhere.
workingDir="${WORKING_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# baseDir is the DATA project. It is site-specific -- override with BASE_DIR=
# rather than editing this file.
baseDir="${BASE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay}"
# resultsRoot is where 1.2 wrote its output. Override with RESULTS_DIR= to keep
# steps 1.2-1.4 pointed at one place.
resultsRoot="${RESULTS_DIR:-${baseDir}/yard/251027_absLoopQuant/results}"
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
