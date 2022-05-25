if ( $null -eq $env:VIRTUAL_ENV ) {
    echo "No active virtual environment, refusing to install."
    exit 1
}

$ErrorActionPreference = "Stop"

$env:CMAKE_PREFIX_PATH = $env:VIRTUAL_ENV + ";" + $env:CMAKE_PREFIX_PATH

pushd $env:Temp

if (Test-Path eigen) {
    rm -r -fo eigen
}
git clone --single-branch --depth=1 --branch 3.4.0 `
    "https://gitlab.com/libeigen/eigen"
pushd eigen
cmake -Bbuild -S. `
    -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
    -D CMAKE_Fortran_COMPILER="NOTFOUND" `
    -D BUILD_TESTING=Off
cmake --build build -j --config RelWithDebInfo
cmake --install build --config RelWithDebInfo
popd

popd
