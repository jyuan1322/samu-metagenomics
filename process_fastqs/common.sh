#!/bin/bash
# =============================================================================
# common.sh — shared helper functions, sourced by the pipeline scripts.
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

# Activate a micromamba environment by name. Requires MAMBA_EXE and
# MAMBA_ROOT_PREFIX (see config.sh). Used inside SLURM job scripts; in an
# interactive local run the environment is assumed already active.
activate_env() {
    local env_name="$1"
    local setup
    setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2>/dev/null || true)"
    if [ -n "$setup" ]; then
        eval "$setup"
    else
        echo "[common.sh] Failed to initialize micromamba" >&2
        exit 1
    fi
    micromamba activate "$env_name" || { echo "[common.sh] Failed to activate $env_name" >&2; exit 1; }
}
