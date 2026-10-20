#!/usr/bin/env bash
# Launch the autoscore Voila app in the sleephd conda environment.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SLEEPFM_ROOT="${SLEEPFM_ROOT:-$HOME/WorkGit/projects/ext/Hypnodensities_TET/external/sleepfm-clinical}"
export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"

# conda activate needs the shell hook; source it then activate
eval "$(conda shell.bash hook)"
conda activate sleephd

echo "Starting autoscore at http://localhost:8866 ..."
echo "SLEEPFM_ROOT=$SLEEPFM_ROOT"
voila "$SCRIPT_DIR/score_edf_voila.ipynb" \
    --port=8866 --strip_sources=True --show_tracebacks=True
