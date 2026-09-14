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
# HOW TO RUN   (workingDir defaults to $PWD)
#
#   sbatch 1.1_calculate_P_s_curves_general.sh
#   bash   1.1_calculate_P_s_curves_general.sh
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
read -r -a subsets <<< "${SUBSETS:-pTh17-1 npTh17 Treg Th1 Th2 Th0}" # Conditions to process
# Per-replicate matrix directory. Expected <dir>/<rep>/cool/<rep>.{mcool,cool}
perReplicateDir="${PER_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro/results/hicpro/hic_results/matrix}"

# Combined-replicate matrix directory. Expected <dir>/<cond>/cool/<cond>.{mcool,cool}
combinedReplicateDir="${COMBINED_REPLICATE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix}"

# Cooler file extensions to look for, in order of preference. A multi-resolution
# .mcool is tried first; a single-resolution .cool is used if no .mcool exists.
# The python side opens whichever it is given correctly.
read -r -a coolExts <<< "${COOL_EXTS:-mcool cool}"

# Output root, shared by steps 1.1 - 1.4.
resultsRoot="${RESULTS_DIR:-$(pwd)/results}"
# Where the P(s) curves land.
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
    
    # Build list of per-replicate cooler files
    perReplicateFiles=()
    for repDir in "${replicateDirs[@]}"; do
        repName=$(basename ${repDir})
        coolFile=""
        for ext in ${coolExts[@]}; do
            cand="${repDir}/cool/${repName}.${ext}"
            [ -f "${cand}" ] && { coolFile="${cand}"; break; }
        done
        if [ -n "${coolFile}" ]; then
            perReplicateFiles+=("${coolFile}")
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
        echo "Warning: No combined .${coolExts[0]}/.${coolExts[1]} file for ${subset} under ${combinedReplicateDir}/${subset}/cool/"
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

