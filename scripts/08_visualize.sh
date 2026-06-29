#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="${1:-config/pipeline_config.sh}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/common.sh"
load_config "${CONFIG_FILE}"
create_output_tree

if [[ "${ENV_BACKEND:-none}" == "conda" && -n "${VISUALIZATION_CONDA_ENV:-}" ]]; then
  export CONDA_ENV="${VISUALIZATION_CONDA_ENV}"
fi

STEP_DIR="090-visualizations"
if is_done "${STEP_DIR}" "visualize"; then
  log "Visualizations already completed; skipping"
  exit 0
fi

run_python_step visualize
mark_done "${STEP_DIR}" "visualize"
