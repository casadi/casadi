#!/bin/bash
# Run every Java example (single-file source launch) against a built casadi.jar.
#   docs/examples/java/run_all.sh <dir with casadi.jar> [example.java ...]
# An example that needs an absent plugin prints "SKIP: <reason>" and exits 0.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
# Windows (Git Bash): Java wants native paths and ';' between class path entries
nat() { if command -v cygpath >/dev/null; then cygpath -m "$1"; else echo "$1"; fi; }
SEP=:; command -v cygpath >/dev/null && SEP=';'
JARDIR=$(nat "$(cd "${1:?usage: run_all.sh <dir with casadi.jar> [example.java ...]}" && pwd)"); shift
cd "$HERE"
[ $# -eq 0 ] && set -- $(find . -name '*.java' | sort)
fail=0
for f in "$@"; do
  echo "== $f =="
  (cd "$(dirname "$f")" && java -Xcheck:jni -cp "$JARDIR/casadi.jar${CLASSPATH:+$SEP$CLASSPATH}" "$(basename "$f")") || { echo "EXAMPLE FAILED: $f"; fail=1; }
done
exit $fail
