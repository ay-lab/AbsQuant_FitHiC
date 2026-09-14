#!/bin/bash

# Create the absloopquantTB mamba environment from absloopquantTB_env.yml.
#
#   bash 0.0_create_absloopquantTB_env.sh
#   ENV_NAME=absloopquantTB2 bash 0.0_create_absloopquantTB_env.sh   # side-by-side
#   ENV_PREFIX=/path/to/envs/absloopquantTB bash 0.0_create_absloopquantTB_env.sh
#
# ENV_PREFIX builds a prefix env at a path you choose instead of a named env in
# the conda root. Prefer it when home is small, quota'd, or flaky.
#
# Requires mamba (or conda) on PATH. Everything resolves from conda-forge and
# bioconda; no cluster, module system or site configuration is needed.

set -uo pipefail

[ -f ~/.bashrc ] && source ~/.bashrc

# Resolve this repo's location rather than hardcoding it, so the script works
# from a clone in any directory.
scriptDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

envName="${ENV_NAME:-absloopquantTB}"
ymlFile="${scriptDir}/absloopquantTB_env.yml"

for arg in "$@"; do
    case "$arg" in
        -h|--help) sed -n '3,14p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

echo "=================================================="
echo "Creating mamba environment: ${envName}"
echo "  spec: $(basename "${ymlFile}")"
echo "=================================================="

if [ ! -f "${ymlFile}" ]; then
    echo "Error: Environment file not found: ${ymlFile}" >&2
    exit 1
fi

if ! command -v mamba >/dev/null 2>&1; then
    echo "Error: mamba not on PATH." >&2
    exit 1
fi

# Build as a prefix env if ENV_PREFIX is set, otherwise a named env.
if [ -n "${ENV_PREFIX:-}" ]; then
    target=(-p "${ENV_PREFIX}")
    exists=$([ -d "${ENV_PREFIX}" ] && echo yes || echo no)
    label="${ENV_PREFIX}"
else
    target=(-n "${envName}")
    exists=$(mamba env list | awk -v n="${envName}" '$1==n {print "yes"; exit}')
    exists="${exists:-no}"
    label="${envName}"
fi

if [ "${exists}" = "yes" ]; then
    echo "Warning: Environment '${label}' already exists"
    if [ "${FORCE_RECREATE:-0}" = "1" ]; then
        REPLY=y
    else
        read -p "Remove it and recreate? (y/n): " -n 1 -r
        echo
    fi
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        mamba env remove "${target[@]}" -y
    else
        echo "Aborting. Environment already exists."
        exit 1
    fi
fi

echo "Creating environment from: ${ymlFile}"
mamba env create "${target[@]}" -f "${ymlFile}" || {
    echo "✗ Error creating environment" >&2
    exit 1
}

# Verify the imports the scripts in this repo actually make, rather than
# trusting that the solve succeeded.
if [ -n "${ENV_PREFIX:-}" ]; then py="${ENV_PREFIX}/bin/python"; else py="$(mamba run "${target[@]}" which python 2>/dev/null)"; fi
echo ""
echo "Verifying imports..."
rc=0
for m in numpy pandas scipy cooler cooltools cv2 matplotlib; do
    if "${py}" -c "import ${m}" >/dev/null 2>&1; then
        printf '  %-12s ok\n' "${m}"
    else
        printf '  %-12s FAIL\n' "${m}"; rc=1
    fi
done
"${py}" -c "import sys; sys.path.insert(0,'${scriptDir}'); import looptools" >/dev/null 2>&1 \
    && printf '  %-12s ok\n' "looptools" || { printf '  %-12s FAIL\n' "looptools"; rc=1; }

echo ""
if [ "${rc}" -eq 0 ]; then
    echo "✓ Successfully created environment: ${label}"
    echo ""
    echo "To activate the environment, run:"
    if [ -n "${ENV_PREFIX:-}" ]; then
        echo "  mamba activate ${ENV_PREFIX}"
    else
        echo "  mamba activate ${envName}"
    fi
else
    echo "✗ Environment created but some imports failed - see above." >&2
    exit 1
fi
