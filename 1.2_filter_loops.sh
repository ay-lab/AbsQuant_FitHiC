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
# HOW TO RUN   (workingDir defaults to $PWD)
#
#   sbatch 1.2_filter_loops.sh
#   bash   1.2_filter_loops.sh
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

read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}" # Conditions to process.

# Chromosomes to write jobs for. Mouse autosomes by default; use chr1..chr22 for human.
# for human, and add chrX / chrY if you want them.
read -r -a chroms <<< "${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19}"

# Per-replicate matrix directory. Expected <dir>/<rep>/cool/<rep>.{mcool,cool}
perReplicateDir="${PER_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"

# Combined-replicate matrix directory. Expected <dir>/<cond>/cool/<cond>.{mcool,cool}
combinedReplicateDir="${COMBINED_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix}"

# Cooler file extensions to look for, in order of preference. A multi-resolution
# .mcool is tried first; a single-resolution .cool is used if no .mcool exists.
# The python side opens whichever it is given correctly.
read -r -a coolExts <<< "${COOL_EXTS:-mcool cool}"

# Filename of the fithic call table.
fithicTemplate="${FITHIC_TEMPLATE:-L20000.U3000000.p2.b200.spline_pass2.res${RESOLUTION:-10000}.significances.txt.gz}"

# Output root, shared by steps 1.1 - 1.4.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"
# P(s) curves written by 1.1.
pScurvesDir="${P_S_CURVES_DIR:-${resultsRoot}/P_s_curves}"

# This repo, holding the .py files.
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
    
    # Build list of per-replicate cooler files and names
    perReplicateFiles=()
    replicateNamesList=()
    for repDir in "${replicateDirs[@]}"; do
        repName=$(basename ${repDir})
        coolFile=""
        for ext in ${coolExts[@]}; do
            cand="${repDir}/cool/${repName}.${ext}"
            [ -f "${cand}" ] && { coolFile="${cand}"; break; }
        done
        if [ -n "${coolFile}" ]; then
            perReplicateFiles+=("${coolFile}")
            replicateNamesList+=("${repName}")
        fi
    done
    
    if [ ${#perReplicateFiles[@]} -eq 0 ]; then
        echo "Warning: No cooler files found for ${subset}, skipping..."
        continue
    fi
    
    # Combined replicate file
    combinedMcoolFile=""
    for ext in ${coolExts[@]}; do
        cand="${combinedReplicateDir}/${subset}/cool/${subset}.${ext}"
        [ -f "${cand}" ] && { combinedMcoolFile="${cand}"; break; }
    done
    
    if [ -z "${combinedMcoolFile}" ]; then
        echo "Warning: No combined cooler (${coolExts[*]}) for ${subset} under ${combinedReplicateDir}/${subset}/cool/"
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

