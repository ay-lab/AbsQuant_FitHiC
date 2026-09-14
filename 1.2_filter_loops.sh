#!/bin/bash

#SBATCH --job-name=filter_loops_per_chromosome
#SBATCH --output=1.2_filter_loops_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Writes one SLURM job script per chromosome per condition for loop filtering.
# This script generates jobs; it does not filter anything itself.
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (run it from inside the repo -- workingDir defaults to $PWD)
#
#   sbatch 1.2_filter_loops.sh
#   bash   1.2_filter_loops.sh
#
# bash is fine for THIS script: it only writes job scripts and takes seconds.
# It is NOT fine for the jobs it generates -- each chromosome runs for hours and
# must be submitted with sbatch. Running them with bash puts the whole filter on
# a login node, where it will be throttled or killed.
#
# Override any INPUT VARIABLE below on the command line:
#
#   SUBSETS="condA condB" FDR_THRESHOLD=0.001 \
#     bash 1.2_filter_loops.sh
#
#   sbatch --export=ALL,SUBSETS="condA condB",FDR_THRESHOLD=0.001 \
#     1.2_filter_loops.sh
# ---------------------------------------------------------------------------

source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================

# Conditions to process; one job script per chromosome is written per entry.
read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}"

# Chromosomes to write jobs for. Mouse autosomes by default; use chr1..chr22
# for human, and add chrX / chrY if you want them.
read -r -a chroms <<< "${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19}"

# Per-replicate matrix directory. Expected layout:
#   ${perReplicateDir}/<condition>-<rep>/cool/<condition>-<rep>.mcool
perReplicateDir="${PER_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"

# Combined-replicate matrix directory. Expected layout:
#   ${combinedReplicateDir}/<condition>/cool/<condition>.mcool
#   ${combinedReplicateDir}/<condition>/fithic/<resolution>/<condition>.<fithicTemplate>
combinedReplicateDir="${COMBINED_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix}"

# Filename of the fithic call table, after the leading "<condition>.".
fithicTemplate="${FITHIC_TEMPLATE:-L20000.U3000000.p2.b200.spline_pass2.res${RESOLUTION:-10000}.significances.txt.gz}"

# Output root, shared by steps 1.1 - 1.4.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"
# P(s) curves written by 1.1. Same default as 1.1 uses, so the two line up.
pScurvesDir="${P_S_CURVES_DIR:-${resultsRoot}/P_s_curves}"

# This repo, holding the .py files. Defaults to the directory you run from.
workingDir="${WORKING_DIR:-$(pwd)}"

resolution="${RESOLUTION:-10000}"
nproc="${NPROC:-30}"
fdrThreshold="${FDR_THRESHOLD:-0.01}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
curr_date=$(date +"%y%m%d")
pythonFile="${workingDir}/1.2_filter_loops.py"
scriptsDir="${workingDir}/qshs/${curr_date}_filter_loops_per_chr_fdr${fdrThreshold}"
resultsDir="${resultsRoot}/filtered_loops_per_chr_fdr${fdrThreshold}"

mkdir -p ${resultsDir}
mkdir -p ${scriptsDir}

# Check if Python script exists
if [ ! -f "${pythonFile}" ]; then
    echo "Error: Python script not found: ${pythonFile}"
    exit 1
fi

echo "=================================================="
echo "Generating loop filtering scripts per chromosome"
echo "=================================================="
echo "Scripts directory: ${scriptsDir}"
echo "Results directory: ${resultsDir}"
echo ""

for subset in ${subsets[@]}; do
    # Find per-replicate directories
    replicateDirs=($(find ${perReplicateDir} -maxdepth 1 -type d -name "${subset}*" | sort))
    
    
    if [ ${#replicateDirs[@]} -eq 0 ]; then
        echo "Warning: No replicates found for ${subset}, skipping..."
        continue
    fi
    
    # Build list of per-replicate .mcool files and names
    perReplicateFiles=()
    replicateNamesList=()
    for repDir in "${replicateDirs[@]}"; do
        repName=$(basename ${repDir})
        mcoolFile="${repDir}/cool/${repName}.mcool"
        if [ -f "${mcoolFile}" ]; then
            perReplicateFiles+=("${mcoolFile}")
            replicateNamesList+=("${repName}")
        fi
    done
    
    if [ ${#perReplicateFiles[@]} -eq 0 ]; then
        echo "Warning: No .mcool files found for ${subset}, skipping..."
        continue
    fi
    
    # Combined replicate file
    combinedMcoolFile="${combinedReplicateDir}/${subset}/cool/${subset}.mcool"
    
    if [ ! -f "${combinedMcoolFile}" ]; then
        echo "Warning: Combined .mcool file not found for ${subset}: ${combinedMcoolFile}"
        echo "  Skipping ${subset}..."
        continue
    fi
    
    # Fithic combined loop file
    # Assuming standard fithic output structure
    fithicFile="${combinedReplicateDir}/${subset}/fithic/${resolution}/${subset}.${fithicTemplate}"
    
    if [ ! -f "${fithicFile}" ]; then
        echo "Warning: Fithic file not found for ${subset}: ${fithicFile}"
        echo "  Skipping ${subset}..."
        continue
    fi
    
    outputDir="${resultsDir}/${subset}"
    mkdir -p ${outputDir}
    scriptsDirChr="${scriptsDir}/${subset}"
    mkdir -p ${scriptsDirChr}
    pScurvesSubDir="${pScurvesDir}/${subset}"
    
    echo "Processing: ${subset}"
    echo "  Per-replicate files: ${#perReplicateFiles[@]}"
    echo "  Combined file: ${combinedMcoolFile}"
    echo "  Fithic file: ${fithicFile}"
    echo "  Output directory: ${outputDir}"
    
    # Create scripts for each chromosome
    for chrom in ${chroms[@]}; do
        cat <<EOF > ${scriptsDirChr}/filter_loops_${subset}_${chrom}.sh
#!/bin/bash
#SBATCH --job-name=filter_loops_${subset}_${chrom}
#SBATCH --output=${scriptsDirChr}/filter_loops_${subset}_${chrom}_%j.out
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=${nproc}
#SBATCH --nodes=1
#SBATCH --mem=300g

source ~/.bashrc
mamba activate absloopquantTB
cd ${workingDir}

python3 ${pythonFile} \\
    --per_replicate_files ${perReplicateFiles[@]} \\
    --combined_replicate_file ${combinedMcoolFile} \\
    --per_replicate_names ${replicateNamesList[@]} \\
    --combined_replicate_name ${subset} \\
    --P_s_curves_dir ${pScurvesSubDir} \\
    --fithic_combined_loop_file ${fithicFile} \\
    --chromosome ${chrom} \\
    --output_dir ${outputDir} \\
    --resolution ${resolution} \\
    --nproc ${nproc} \\
    --fdr_threshold ${fdrThreshold} \\
    --verbose

EOF
    done
    
    chmod +x ${scriptsDirChr}/filter_loops_${subset}_*.sh
    echo "  ✓ Created scripts for ${#chroms[@]} chromosomes"
    echo ""
done

echo "=================================================="
echo "Generated scripts in: ${scriptsDir}"
echo "Submit jobs using:"
echo "  sbatch ${scriptsDir}/*/filter_loops_*.sh"
echo "=================================================="

