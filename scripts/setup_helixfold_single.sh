#!/usr/bin/env bash
set -euo pipefail

# Clones/updates HelixFold-single and downloads the trained model.
# Dependency installation is done manually in a clean Python 3.11 venv.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

HELIXFOLD_REPO_URL="https://github.com/EricwanAR/HelixFold-single.git"
HELIXFOLD_MODEL_URL="https://baidu-nlp.bj.bcebos.com/PaddleHelix/HelixFold-Single/helixfold-single.pdparams"

REPO_DIR="${REPO_ROOT}/backend/external/HelixFold-single"
MODEL_PATH="${REPO_ROOT}/backend/models/helixfold-single.pdparams"

usage() {
  cat <<'EOF'
Usage: setup_helixfold_single.sh [options]

Options:
  --repo-dir <path>      Target HelixFold-single repository directory.
  --model-path <path>    Target model .pdparams file path.
  -h, --help             Show this help message.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo-dir)
      REPO_DIR="$2"
      shift 2
      ;;
    --model-path)
      MODEL_PATH="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

download_file() {
  local url="$1"
  local out="$2"

  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 --retry-delay 1 -o "$out" "$url"
    return
  fi

  if command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
    return
  fi

  echo "Neither curl nor wget is available to download: $url" >&2
  exit 1
}

need_cmd git
need_cmd python3

REPO_DIR="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$REPO_DIR")"
MODEL_PATH="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$MODEL_PATH")"

mkdir -p "$(dirname "$REPO_DIR")"
mkdir -p "$(dirname "$MODEL_PATH")"

if [[ -d "$REPO_DIR" ]]; then
  if [[ -d "$REPO_DIR/.git" ]]; then
    echo "[1/2] HelixFold-single repo already present at $REPO_DIR"
  else
    echo "[1/2] HelixFold-single directory already exists at $REPO_DIR"
  fi
else
  echo "[1/2] Cloning HelixFold-single into $REPO_DIR"
  git clone "$HELIXFOLD_REPO_URL" "$REPO_DIR"
fi

if [[ -f "$MODEL_PATH" ]]; then
  echo "[2/2] HelixFold model already present at $MODEL_PATH"
else
  echo "[2/2] Downloading model to $MODEL_PATH"
  download_file "$HELIXFOLD_MODEL_URL" "$MODEL_PATH"
fi

cat <<EOF

Done.

Create a clean Python 3.11 venv manually, then install HelixFold requirements there.
Before installing, remove version pins from HelixFold's requirements.txt and add paddlepaddle.

Suggested next steps:
  cd "$REPO_DIR"
  python3.11 -m venv .venv
  source .venv/bin/activate
  python -m pip install --upgrade pip
  # edit requirements.txt: drop version pins, add paddlepaddle
  python -m pip install -r requirements.txt

Add these values to backend/.env:

HELIXFOLD_SINGLE_ENABLED=true
HELIXFOLD_SINGLE_REPO_DIR=$REPO_DIR
HELIXFOLD_SINGLE_MODEL_PATH=$MODEL_PATH
HELIXFOLD_SINGLE_PYTHON_BIN=$REPO_DIR/.venv/bin/python
HELIXFOLD_SINGLE_SCRIPT_RELPATH=helixfold_single_inference.py
HELIXFOLD_SINGLE_TIMEOUT_SECONDS=7200
EOF
