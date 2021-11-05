#!/usr/bin/env bash
cd "$(dirname "${BASH_SOURCE[0]}")"
if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
pip uninstall -y alpaqa
./wheel-local.sh
pip install --find-links=../../wheelhouse alpaqa
