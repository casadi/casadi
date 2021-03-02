#!/usr/bin/env bash

root_dir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"/..

export CUTEST="${root_dir}/py-venv/opt/CUTEst/cutest"
export SIFDECODE="${root_dir}/py-venv/opt/CUTEst/sifdecode"
export ARCHDEFS="${root_dir}/py-venv/opt/CUTEst/archdefs"
export MASTSIF="${root_dir}/py-venv/opt/CUTEst/sif"
export MYARCH="pc64.lnx.gfo"

prob_name="${1:-ROSENBR}"

tmp_dir=$(mktemp -d)
pushd "${tmp_dir}"
"${root_dir}/py-venv/opt/CUTEst/sifdecode/bin/sifdecoder" ${prob_name}
LD_LIBRARY_PATH="${root_dir}/py-venv/lib" \
    "${root_dir}/py-venv/opt/CUTEst/cutest/bin/runcutest" \
        -A $MYARCH \
        -p ipopt \
        -L"${root_dir}/py-venv/lib"
popd
rm -rf "${tmp_dir}"
