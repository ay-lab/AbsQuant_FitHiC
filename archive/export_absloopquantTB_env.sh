#!/bin/bash

# Script to export absloopquantTB mamba environment to a .yml file

source ~/.bashrc

# Environment name
envName="absloopquantTB"
scriptDir="/home/bbabatunde/packages/25-09-absloopquant/AbsLoopQuant_TB"
outputYml="${scriptDir}/absloopquantTB_env.yml"

echo "=================================================="
echo "Exporting mamba environment: ${envName}"
echo "=================================================="

# Check if environment exists
if ! mamba env list | grep -q "[[:space:]]*${envName}[[:space:]]"; then
    echo "Error: Environment '${envName}' not found"
    echo "Available environments:"
    mamba env list
    exit 1
fi

# Export environment to yml file
echo "Exporting environment to: ${outputYml}"
mamba env export -n ${envName} > ${outputYml}

if [ $? -eq 0 ]; then
    echo "✓ Successfully exported environment to: ${outputYml}"
    echo ""
    echo "To recreate this environment, run:"
    echo "  bash ${scriptDir}/create_absloopquantTB_env.sh"
else
    echo "✗ Error exporting environment"
    exit 1
fi

