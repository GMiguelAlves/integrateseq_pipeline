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

activate_runtime
require_cmd "${RSCRIPT_BIN}"
"${RSCRIPT_BIN}" "${SCRIPT_DIR}/r/visualize_integrative.R" \
  --project-dir "${INTEGRATION_OUTPUT_DIR}" \
  --outdir "${INTEGRATION_OUTPUT_DIR}/${STEP_DIR}"
mark_done "${STEP_DIR}" "visualize"
