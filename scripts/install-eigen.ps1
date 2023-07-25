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
        git clone --single-branch --depth=1 --branch master `
            "https://gitlab.com/libeigen/eigen"
    }
    pushd eigen
    if (Test-Path build/CMakeCache.txt) {
        rm build/CMakeCache.txt
    }
    cmake -Bbuild -S. `
        -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
        -D CMAKE_Fortran_COMPILER="NOTFOUND" `
        -D BUILD_TESTING=Off
    cmake --build build -j --config "$config"
    cmake --install build --config "$config"
    popd # eigen

    popd # $srcdir
} catch {
    throw
}
