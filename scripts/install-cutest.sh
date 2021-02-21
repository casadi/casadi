#!/usr/bin/env bash
cd "$( dirname "${BASH_SOURCE[0]}" )"/..

if [ -z "${VIRTUAL_ENV+x}" ]; then
    echo "No active virtual environment, refusing to install."
    exit 1
fi

export CUTEST_ROOT="${VIRTUAL_ENV}/opt/CUTEst"
mkdir -p "$CUTEST_ROOT"
cd "$CUTEST_ROOT"

rm -rf archdefs
rm -rf sifdecode
rm -rf cutest
[ -d archdefs ] || git clone https://github.com/ralna/ARCHDefs --depth 1 ./archdefs
[ -d sifdecode ] || git clone https://github.com/ralna/SIFDecode --depth 1 ./sifdecode
[ -d cutest ] || git clone https://github.com/ralna/CUTEst --depth 1 ./cutest
[ -d sif ] || git clone https://bitbucket.org/optrove/sif --depth 1 ./sif

export ARCH="$CUTEST_ROOT/archdefs"
export ARCHDEFS="$CUTEST_ROOT/archdefs"
export CCOMPILER="ccompiler.pc64.lnx.gcc"
export COMPILER="compiler.pc64.lnx.gfo"
export CMP="gfo"
export CUTEST="$CUTEST_ROOT/cutest"
export GALAHAD_REQPREC="D"
export MACHINE="Intel-like PC with a 64-bit processor"
export MCH="pc64"
export OPSYS="Linux"
export OS="lnx"
export SIFDECODE="$CUTEST_ROOT/sifdecode"
export VERSION="${MCH}.${OS}.${CMP}"

source $ARCH/system.$OS
source $ARCH/$COMPILER
source $ARCH/$CCOMPILER

yesno_default_no() {
    echo "$1 (y/N)? "
    return 0
}
yesno_default_yes() {
    echo "$1 (Y/n)? "
    return 1
}
warning() {
    tput bold
    echo -e "\n WARNING: $1\n"
    tput sgr0
}
error() {
    tput bold
    echo -e "\n ERROR: $1\n"
    tput sgr0
}
success() {
    tput bold
    echo -e "\n SUCCESS: $1\n"
    tput sgr0
}
message() {
    tput bold
    echo -e "\n MESSAGE: $1\n"
    tput sgr0
}

COMPUSED=`$LS $ARCH/compiler.${MCH}.${OS}.gfo 2>/dev/null`
CCOMPUSED=`$LS $ARCH/ccompiler.${MCH}.${OS}.gcc 2>/dev/null`

dirs -c

echo ' Installing SIFDecode ...'
cd "$SIFDECODE"
source ./bin/install_sifdecode_main
if [[ ! $? ]]; then
    error 'An error occurred while installing SIFDecode.'
    exit $?
fi

echo ' Installing CUTEst ...'
cd "$CUTEST"
source ./bin/install_cutest_main
if [[ ! $? ]]; then
    error 'An error occurred while installing CUTEst.'
    exit $?
fi
