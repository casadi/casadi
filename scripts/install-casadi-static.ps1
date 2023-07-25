param([String]$config="Release")

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

    if ( -not (Test-Path casadi) ) {
        git clone --branch "3.6.3" --depth 1 --recursive `
            https://github.com/casadi/casadi
    }
    pushd casadi
    if (Test-Path build/CMakeCache.txt) {
        rm build/CMakeCache.txt
    }
    cmake -Bbuild -S. `
        -D CMAKE_INSTALL_PREFIX="$env:VIRTUAL_ENV" `
        -D CMAKE_Fortran_COMPILER="NOTFOUND" `
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
    cmake --build build -j --config $config
    cmake --install build --config $config
    popd # casadi

    popd # $srcdir
} catch {
    throw
}
