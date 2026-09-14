#!/bin/bash

# Script to export absloopquantTB mamba environment to a .yml file

source ~/.bashrc

# Environment name
envName="absloopquantTB"
# Resolve the repo root from this script's location (archive/ -> repo root).
scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Writes the LOCK file, not absloopquantTB_env.yml. The latter is the
# hand-maintained portable spec and must not be clobbered by an export.
outputYml="${scriptDir}/absloopquantTB_env.lock.yml"

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
    echo "This is a build-pinned, platform-specific record of an already-built"
    echo "environment. It is gitignored and is NOT how users install the toolkit;"
    echo "they use absloopquantTB_env.yml. To rebuild from this lock file:"
    echo "  mamba env create -n absloopquantTB -f ${outputYml}"
else
    echo "✗ Error exporting environment"
    exit 1
fi

