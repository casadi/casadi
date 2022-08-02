if ( $null -eq $env:VIRTUAL_ENV ) {
    echo "No active virtual environment, refusing to install."
    exit 1
}

$ErrorActionPreference = "Stop"

$env:CMAKE_PREFIX_PATH = $env:VIRTUAL_ENV + ";" + $env:CMAKE_PREFIX_PATH

pushd $env:Temp

if ( -not ( Test-Path casadi) ) {
    git clone --branch "3.5.5" --depth 1 --recursive `
        https://github.com/casadi/casadi
}
pushd casadi
cmake -Bbuild -S. `
    -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
    -D CMAKE_POSITION_INDEPENDENT_CODE=On `
    -D WITH_COMMON=Off `
    -D WITH_PYTHON=Off `
    -D WITH_PYTHON3=Off `
    -D WITH_OPENMP=Off `
    -D WITH_THREAD=On `
    -D WITH_DL=On `
    -D WITH_IPOPT=Off `
    -D ENABLE_STATIC=On `
    -D ENABLE_SHARED=Off
cmake --build build -j --config Release
cmake --install build --config Release
popd

popd
