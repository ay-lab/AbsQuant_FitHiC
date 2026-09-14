#!/bin/bash

# Script to create absloopquantTB mamba environment from .yml file

source ~/.bashrc

# Environment name
envName="absloopquantTB"
scriptDir="/home/bbabatunde/packages/25-09-absloopquant/AbsLoopQuant_TB"
ymlFile="${scriptDir}/absloopquantTB_env.yml"

echo "=================================================="
echo "Creating mamba environment: ${envName}"
echo "=================================================="

# Check if yml file exists
if [ ! -f "${ymlFile}" ]; then
    echo "Error: Environment file not found: ${ymlFile}"
    echo "Please export the environment first using:"
    echo "  bash ${scriptDir}/export_absloopquantTB_env.sh"
    exit 1
fi

# Check if environment already exists
if mamba env list | grep -q "[[:space:]]*${envName}[[:space:]]"; then
    echo "Warning: Environment '${envName}' already exists"
    read -p "Do you want to remove it and recreate? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        mamba env remove -n ${envName} -y
    else
        echo "Aborting. Environment already exists."
        exit 1
    fi
fi

# Create environment from yml file
echo "Creating environment from: ${ymlFile}"
mamba env create -n ${envName} -f ${ymlFile}

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Successfully created environment: ${envName}"
    echo ""
    echo "To activate the environment, run:"
    echo "  mamba activate ${envName}"
else
    echo "✗ Error creating environment"
    exit 1
fi

