param([String]$config="RelWithDebInfo")

if ([string]::IsNullOrEmpty($env:VIRTUAL_ENV)) {
    echo "No active virtual environment, refusing to install."
    exit 1
}

$env:CMAKE_PREFIX_PATH = $env:VIRTUAL_ENV + ";" + $env:CMAKE_PREFIX_PATH

$ErrorActionPreference = "Stop"
try {
    $srcdir=$env:VIRTUAL_ENV + "\src"
    mkdir -fo $srcdir
    pushd $srcdir

    if ( -not (Test-Path googletest) ) {
        git clone --single-branch --depth=1 --branch main `
            "https://github.com/google/googletest.git"
    }
    pushd googletest
    if (Test-Path build/CMakeCache.txt) {
        rm build/CMakeCache.txt
    }
    cmake -Bbuild -S. `
        -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
        -D gtest_force_shared_crt=On
    cmake --build build -j --config "$config"
    cmake --install build --config "$config"
    popd # googletest

    popd # $srcdir
} catch {
    throw
}
