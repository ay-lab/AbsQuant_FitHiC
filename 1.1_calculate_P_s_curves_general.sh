#!/bin/bash

#SBATCH --job-name=calculate_P_s_curves
#SBATCH --output=1.1_calculate_P_s_curves_general_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Writes one SLURM job script per condition for the P(s) curve calculation.
# This script generates jobs; it does not compute anything itself.
#
# ---------------------------------------------------------------------------
# HOW TO RUN   (run it from inside the repo -- workingDir defaults to $PWD)
#
#   sbatch 1.1_calculate_P_s_curves_general.sh
#   bash   1.1_calculate_P_s_curves_general.sh
#
# bash is fine for THIS script: it only writes job scripts and takes seconds.
# It is NOT fine for the jobs it generates -- those run for hours per sample and
# must be submitted with sbatch. Running them with bash puts the whole
# calculation on a login node, where it will be throttled or killed.
#
# Override any INPUT VARIABLE below on the command line:
#
#   SUBSETS="condA condB" PER_REPLICATE_DIR=/path/to/matrix \
#     bash 1.1_calculate_P_s_curves_general.sh
#
#   sbatch --export=ALL,SUBSETS="condA condB",PER_REPLICATE_DIR=/path/to/matrix \
#     1.1_calculate_P_s_curves_general.sh
# ---------------------------------------------------------------------------

source ~/.bashrc

# ===========================================================================
# INPUT VARIABLES
# ===========================================================================

# Conditions to process; one job script is written per entry.
read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}"

# Per-replicate matrix directory. Expected layout:
#   ${perReplicateDir}/<condition>-<rep>/cool/<condition>-<rep>.mcool
perReplicateDir="${PER_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"

# Combined-replicate matrix directory. Expected layout:
#   ${combinedReplicateDir}/<condition>/cool/<condition>.mcool
combinedReplicateDir="${COMBINED_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix}"

# Output root, shared by steps 1.1 - 1.4. Set RESULTS_DIR once and the whole
# chain stays in one place.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"
# Where the P(s) curves land. Step 1.2 reads this same default.
pScurvesDir="${P_S_CURVES_DIR:-${resultsRoot}/P_s_curves}"

# This repo, holding the .py files. Defaults to the directory you run from.
workingDir="${WORKING_DIR:-$(pwd)}"

resolution="${RESOLUTION:-10000}"
nproc="${NPROC:-30}"

# ===========================================================================
# Derived - no need to edit below here
# ===========================================================================
curr_date=$(date +"%y%m%d")
pythonFile="${workingDir}/1.1_calculate_P_s_curves_general.py"
scriptsDir="${workingDir}/qshs/${curr_date}_calculate_P_s_curves"
resultsDir="${pScurvesDir}"

mkdir -p ${resultsDir}
mkdir -p ${scriptsDir}

# Check if Python script exists
if [ ! -f "${pythonFile}" ]; then
    echo "Error: Python script not found: ${pythonFile}"
    exit 1
fi

echo "=================================================="
echo "Generating P(s) curve calculation scripts"
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
    
    # Build list of per-replicate .mcool files
    perReplicateFiles=()
    for repDir in "${replicateDirs[@]}"; do
        repName=$(basename ${repDir})
        mcoolFile="${repDir}/cool/${repName}.mcool"
        if [ -f "${mcoolFile}" ]; then
            perReplicateFiles+=("${mcoolFile}")
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
    
    outputDir="${resultsDir}/${subset}"
    mkdir -p ${outputDir}
    
    echo "Processing: ${subset}"
    echo "  Per-replicate files: ${#perReplicateFiles[@]}"
    echo "  Combined file: ${combinedMcoolFile}"
    echo "  Output directory: ${outputDir}"
    
    # Create SLURM script
    cat <<EOF > ${scriptsDir}/P_s_curves_${subset}.sh
#!/bin/bash
#SBATCH --job-name=P_s_curves-${subset}
#SBATCH --output=${scriptsDir}/P_s_curves-${subset}-%j.out
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
    --output_dir ${outputDir} \\
    --resolution ${resolution} \\
    --nproc ${nproc} \\
    --verbose

EOF
    
    chmod +x ${scriptsDir}/P_s_curves_${subset}.sh
    echo "  ✓ Created script: P_s_curves_${subset}.sh"
    echo ""
done

echo "=================================================="
echo "Generated scripts in: ${scriptsDir}"
echo "Submit jobs using:"
echo "  sbatch ${scriptsDir}/*.sh"
echo "=================================================="

