#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
pushd /tmp

[ -d sphinx ] || \
    git clone https://github.com/tttapa/sphinx --branch 4.x --depth 1
pushd sphinx
python setup.py install
popd

popd
