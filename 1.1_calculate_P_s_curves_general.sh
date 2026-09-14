#!/bin/bash

#SBATCH --job-name=calculate_P_s_curves
#SBATCH --output=1.1_calculate_P_s_curves_general_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Script to generate SLURM scripts for calculating P(s) curves
# Creates example scripts in qshs directory

source ~/.bashrc

# Parameters
curr_date=$(date +"%y%m%d")
subsets=("pTh17-1" "npTh17" "Treg" "Th1" "Th2" "Th0")
resolution=10000
nproc=30

# Directories
# workingDir is this repo. Derived from the script's own location so a clone
# works anywhere; override with WORKING_DIR= if you keep the code elsewhere.
workingDir="${WORKING_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# baseDir is the DATA project. It is site-specific -- override with BASE_DIR=
# rather than editing this file.
baseDir="${BASE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay}"
perReplicateDir="${baseDir}/yard/251014_HiCPro/results/hicpro/hic_results/matrix"
combinedReplicateDir="${baseDir}/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix"
pythonFile="${workingDir}/1.1_calculate_P_s_curves_general.py"
scriptsDir="${workingDir}/qshs/${curr_date}_calculate_P_s_curves"
resultsDir="${workingDir}/results/P_s_curves"

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

