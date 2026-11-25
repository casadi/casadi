#!/bin/bash
# Run every C# example against a casadi build:
#   docs/examples/csharp/run_all.sh <build>/swig/dotnet [example ...]
# An example that needs an absent plugin prints "SKIP: <reason>" and exits 0.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
DIR=$(cd "${1:?usage: run_all.sh <build>/swig/dotnet or <install>/casadi [example ...]}" && pwd); shift
# A build dir (<build>/swig/dotnet: assembly/CasADi.dll, libs in <build>/lib) or an install dir (CasADi.dll next to the libs)
ASM="$DIR/assembly/CasADi.dll"; [ -f "$DIR/CasADi.dll" ] && ASM="$DIR/CasADi.dll"
# Git Bash on Windows: dotnet and casadi want Windows paths, native libraries are found on PATH
win() { if command -v cygpath >/dev/null 2>&1; then cygpath -m "$1"; else printf '%s' "$1"; fi; }
export DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1
OUT=${CASADI_DOTNET_EXAMPLES_OUT:-$(mktemp -d)}
dotnet build "$HERE/Examples.csproj" -c Release -nologo -v q -p:CasadiAssembly="$(win "$ASM")" -o "$(win "$OUT")" \
  -p:BaseIntermediateOutputPath="$(win "$OUT")/obj/"
LIBDIR=$(cd "$DIR/../../lib" 2>/dev/null && pwd || echo "$DIR")
export LD_LIBRARY_PATH="$OUT:$LIBDIR:$DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$LD_LIBRARY_PATH"
export PATH="$OUT:$LIBDIR:$DIR:$PATH"
cd "$HERE"
[ $# -eq 0 ] && set -- $(ls *.cs | grep -v ExampleRunner.cs | sed 's/\.cs$//' | sort)
fail=0
for e in "$@"; do
  echo "== $e =="
  dotnet "$(win "$OUT")/Examples.dll" "$e" || { echo "EXAMPLE FAILED: $e"; fail=1; }
done
exit $fail
