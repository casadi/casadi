#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ] && [ $1 != "--force-install" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

set -ex
pushd /tmp

[ -d breathe ] || git clone https://github.com/michaeljones/breathe
pushd breathe
git fetch origin pull/711/head:pr711 ||:
git checkout pr711
python setup.py install
popd

popd
