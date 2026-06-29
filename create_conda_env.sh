#!/usr/bin/env bash

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${PIPELINE_DIR}/environment.yaml"
ENV_NAME="${ENV_NAME:-ebv-typing}"
FRONTEND="${CONDA_FRONTEND:-}"
SKIP_VERIFY="${SKIP_VERIFY:-0}"

usage() {
  cat <<'EOF'
Create or update the conda environment for the EBV typing Snakemake pipeline.

Usage:
  bash create_conda_env.sh
  ENV_NAME=my-ebv-env bash create_conda_env.sh
  CONDA_FRONTEND=mamba bash create_conda_env.sh
  SKIP_VERIFY=1 bash create_conda_env.sh

After creation:
  conda activate ebv-typing
  snakemake --cores 8 --configfile config/config.yaml

Notes:
  - The environment is intentionally all-in-one, so you can run without
    --use-conda and avoid nested rule-environment solver issues.
  - If you prefer Snakemake rule environments, activate this environment first
    and run:
      snakemake --use-conda --cores 8 --configfile config/config.yaml
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

detect_frontend() {
  if [[ -n "${FRONTEND}" ]]; then
    command -v "${FRONTEND}" >/dev/null 2>&1 || {
      echo "ERROR: CONDA_FRONTEND=${FRONTEND} was requested but is not on PATH." >&2
      exit 1
    }
    printf '%s\n' "${FRONTEND}"
    return
  fi

  if command -v mamba >/dev/null 2>&1; then
    printf 'mamba\n'
  elif command -v conda >/dev/null 2>&1; then
    printf 'conda\n'
  elif command -v micromamba >/dev/null 2>&1; then
    printf 'micromamba\n'
  else
    echo "ERROR: Could not find mamba, conda, or micromamba on PATH." >&2
    exit 1
  fi
}

FRONTEND="$(detect_frontend)"

env_exists() {
  "${FRONTEND}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"
}

run_in_env() {
  if [[ "${FRONTEND}" == "micromamba" ]]; then
    micromamba run -n "${ENV_NAME}" "$@"
  else
    "${FRONTEND}" run -n "${ENV_NAME}" "$@"
  fi
}

echo "[setup] Pipeline directory: ${PIPELINE_DIR}"
echo "[setup] Environment file:   ${ENV_FILE}"
echo "[setup] Environment name:   ${ENV_NAME}"
echo "[setup] Conda frontend:     ${FRONTEND}"

if [[ ! -s "${ENV_FILE}" ]]; then
  echo "ERROR: Missing environment file: ${ENV_FILE}" >&2
  exit 1
fi

if [[ "${FRONTEND}" != "micromamba" ]]; then
  "${FRONTEND}" config --set channel_priority strict >/dev/null 2>&1 || true
fi

if env_exists; then
  echo "[setup] Existing environment found; updating with --prune."
  "${FRONTEND}" env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
  echo "[setup] Creating environment."
  "${FRONTEND}" env create -n "${ENV_NAME}" -f "${ENV_FILE}"
fi

if [[ "${SKIP_VERIFY}" != "1" ]]; then
  echo "[verify] Checking command-line tools."
  run_in_env bash -lc '
    set -euo pipefail
    for tool in snakemake fastp hisat2 hisat2-build samtools bcftools bedtools featureCounts tabix; do
      command -v "$tool" >/dev/null
      printf "  OK  %s -> %s\n" "$tool" "$(command -v "$tool")"
    done
  '

  echo "[verify] Checking Python imports."
  run_in_env python - <<'PY'
import Bio
import mappy
import matplotlib
import numpy
import pandas
import pysam
import seaborn

print("  OK  Python imports: Bio, mappy, matplotlib, numpy, pandas, pysam, seaborn")
PY

  echo "[verify] Snakemake version:"
  run_in_env snakemake --version
fi

cat <<EOF

[done] Environment is ready.

Activate it:
  conda activate ${ENV_NAME}

Recommended run command from ${PIPELINE_DIR}:
  cd ${PIPELINE_DIR}
  snakemake --cores 8 --configfile config/config.yaml

If you want Snakemake to create per-rule environments anyway:
  snakemake --use-conda --cores 8 --configfile config/config.yaml
EOF

