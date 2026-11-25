#!/bin/bash
# Build and run the C# tests against a casadi build:
#   test/dotnet/run.sh <build>/swig/dotnet [test name substrings...]
# <dir> is <build>/swig/dotnet or an installed casadi dir (CasADi.dll, casadi_dotnet and casadi side by side).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
DIR=$(cd "${1:?usage: run.sh <build>/swig/dotnet or <install>/casadi [tests...]}" && pwd); shift
# A build dir (<build>/swig/dotnet: assembly/CasADi.dll, libs in <build>/lib) or an install dir (CasADi.dll next to the libs)
ASM="$DIR/assembly/CasADi.dll"; [ -f "$DIR/CasADi.dll" ] && ASM="$DIR/CasADi.dll"
# Git Bash on Windows: dotnet and casadi want Windows paths, native libraries are found on PATH
win() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s' "$1"; fi; }
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
export CASADI_TEST_DATA=${CASADI_TEST_DATA:-$(cd "$HERE/../data" 2>/dev/null && win "$(pwd)")}
OUT=${CASADI_DOTNET_TEST_OUT:-$(mktemp -d)}
dotnet build "$HERE/Tests.csproj" -c Release -nologo -v q -p:CasadiAssembly="$(win "$ASM")" -o "$(win "$OUT")" \
  -p:BaseIntermediateOutputPath="$(win "$OUT")/obj/"
LIBDIR=$(cd "$DIR/../../lib" 2>/dev/null && pwd || echo "$DIR")
export LD_LIBRARY_PATH="$OUT:$LIBDIR:$DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$LD_LIBRARY_PATH"
export PATH="$OUT:$LIBDIR:$DIR:$PATH"
dotnet "$(win "$OUT")/Tests.dll" "$@"
