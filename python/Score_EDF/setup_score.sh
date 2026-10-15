#!/usr/bin/env bash
# Create (or update) the autoscore conda environment.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_NAME="autoscore"
ENV_FILE="$SCRIPT_DIR/environment_autoscore.yml"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Updating existing '$ENV_NAME' environment..."
  conda env update --name "$ENV_NAME" --file "$ENV_FILE" --prune
else
  echo "Creating '$ENV_NAME' environment..."
  conda env create --file "$ENV_FILE"
fi

# SleepFM-Clinical checkout
SLEEPFM_ROOT="${SLEEPFM_ROOT:-}"
if [ -z "$SLEEPFM_ROOT" ]; then
  echo ""
  echo "NOTE: Set SLEEPFM_ROOT to your sleepfm-clinical checkout to enable SleepFM."
  echo "  export SLEEPFM_ROOT=/path/to/sleepfm-clinical"
  echo ""
else
  echo "SLEEPFM_ROOT=$SLEEPFM_ROOT"
  if [ ! -d "$SLEEPFM_ROOT/sleepfm/checkpoints" ]; then
    echo "WARNING: checkpoints not found at $SLEEPFM_ROOT/sleepfm/checkpoints"
  fi
fi

echo ""
echo "Ready. Launch with:  bash autoscore.sh"
