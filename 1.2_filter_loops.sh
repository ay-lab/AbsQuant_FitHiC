#!/bin/bash

#SBATCH --job-name=filter_loops_per_chromosome
#SBATCH --output=1.2_filter_loops_%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150g

# Script to generate SLURM scripts for filtering loops per chromosome
# Creates example scripts in qshs directory

source ~/.bashrc

# Parameters
curr_date=$(date +"%y%m%d")
subsets=("pTh17-1" "npTh17" "Treg" "Th1" "Th2" "Th0")
resolution=10000
nproc=30
fdrThreshold=0.01

# Directories
# workingDir is this repo. Derived from the script's own location so a clone
# works anywhere; override with WORKING_DIR= if you keep the code elsewhere.
workingDir="${WORKING_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# baseDir is the DATA project. It is site-specific -- override with BASE_DIR=
# rather than editing this file.
baseDir="${BASE_DIR:-/mnt/BioAdHoc/Groups/vd-ay/bbabatunde/projects/25-06-Kuchroo-Ay}"
perReplicateDir="${baseDir}/yard/251014_HiCPro/results/hicpro/hic_results/matrix"
combinedReplicateDir="${baseDir}/yard/251014_HiCPro_Combined/results/hicpro/hic_results/matrix"
pythonFile="${workingDir}/1.2_filter_loops.py"
scriptsDir="${workingDir}/qshs/${curr_date}_filter_loops_per_chr_fdr${fdrThreshold}"
resultsDir="${baseDir}/yard/251027_absLoopQuant/results/filtered_loops_per_chr_fdr${fdrThreshold}"
pScurvesDir="${baseDir}/yard/251027_absLoopQuant/results/P_s_curves"

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
    fithicFile="${combinedReplicateDir}/${subset}/fithic/${resolution}/${subset}.L20000.U3000000.p2.b200.spline_pass2.res${resolution}.significances.txt.gz"
    
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
    for chrom in {1..19}; do
        cat <<EOF > ${scriptsDirChr}/filter_loops_${subset}_chr${chrom}.sh
#!/bin/bash
#SBATCH --job-name=filter_loops_${subset}_chr${chrom}
#SBATCH --output=${scriptsDirChr}/filter_loops_${subset}_chr${chrom}_%j.out
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
    --chromosome chr${chrom} \\
    --output_dir ${outputDir} \\
    --resolution ${resolution} \\
    --nproc ${nproc} \\
    --fdr_threshold ${fdrThreshold} \\
    --verbose

EOF
    done
    
    chmod +x ${scriptsDirChr}/filter_loops_${subset}_chr*.sh
    echo "  ✓ Created scripts for chromosomes 1-19"
    echo ""
done

echo "=================================================="
echo "Generated scripts in: ${scriptsDir}"
echo "Submit jobs using:"
echo "  sbatch ${scriptsDir}/*/filter_loops_*.sh"
echo "=================================================="

