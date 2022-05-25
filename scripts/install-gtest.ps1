param([String]$config="RelWithDebInfo")

if ( $null -eq $env:VIRTUAL_ENV ) {
    echo "No active virtual environment, refusing to install."
    exit 1
}

$ErrorActionPreference = "Stop"

$env:CMAKE_PREFIX_PATH = $env:VIRTUAL_ENV + ";" + $env:CMAKE_PREFIX_PATH

pushd $env:Temp

if (Test-Path googletest) {
    rm -r -fo googletest
}
git clone --single-branch --depth=1 --branch main `
    "https://github.com/google/googletest.git"
pushd googletest
cmake -Bbuild -S. `
    -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
    -D gtest_force_shared_crt=On
cmake --build build -j --config $config
cmake --install build --config $config
popd

popd
