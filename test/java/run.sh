#!/bin/bash
# Compile and run the Java tests (JUnit 5) against a built casadi.jar; CASADI_JAVA_RUN_SLOW=1 adds the slow ones.
#   test/java/run.sh <dir with casadi.jar> [junit console args, e.g. --select-class org.casadi.test.SparsityTest]
# The native library is found by NativeLoader (next to casadi.jar, CASADI_JAVA_LIBRARY_PATH, or a natives jar on CLASSPATH).
set -e
# Windows (Git Bash): Java wants native paths and ';' between class path entries
nat() { if command -v cygpath >/dev/null; then cygpath -m "$1"; else echo "$1"; fi; }
SEP=:; command -v cygpath >/dev/null && SEP=';'
HERE=$(nat "$(cd "$(dirname "$0")" && pwd)")
JARDIR=$(nat "$(cd "${1:?usage: run.sh <dir with casadi.jar> [junit args]}" && pwd)"); shift
V=1.11.4
SHA1=44650aaff63a4dbb63cb0bf0f3725a4eb9406da1
CACHE=$(nat "${CASADI_JAVA_TEST_CACHE:-$HOME/.cache/casadi-java-test}")
JUNIT=$CACHE/junit-platform-console-standalone-$V.jar
if [ ! -f "$JUNIT" ]; then
  mkdir -p "$CACHE"
  curl -fsSL -o "$JUNIT.part" https://repo1.maven.org/maven2/org/junit/platform/junit-platform-console-standalone/$V/junit-platform-console-standalone-$V.jar
  echo "$SHA1  $JUNIT.part" | sha1sum -c - >/dev/null
  mv "$JUNIT.part" "$JUNIT"
fi
CP="$JARDIR/casadi.jar${CLASSPATH:+$SEP$CLASSPATH}"
OUT=$(nat "${CASADI_JAVA_TEST_OUT:-$(mktemp -d)}")
mkdir -p "$OUT/classes"
# CASADI_JAVA_TEST_FILES="FooTest.java ..." compiles only those (plus CasadiTestCase.java)
if [ -n "$CASADI_JAVA_TEST_FILES" ]; then
  SRC="$HERE/src/org/casadi/test/CasadiTestCase.java"
  for f in $CASADI_JAVA_TEST_FILES; do SRC="$SRC $HERE/src/org/casadi/test/$f"; done
else
  SRC=$(find "$HERE/src" -name '*.java')
fi
javac --release 11 -encoding UTF-8 -Xlint:all -cp "$CP$SEP$JUNIT" -d "$OUT/classes" $SRC
[ $# -eq 0 ] && set -- --scan-class-path "$OUT/classes"
cd "$HERE"
exec java -Xcheck:jni ${JAVA_OPTS} -jar "$JUNIT" execute --disable-banner --disable-ansi-colors --details=${JUNIT_DETAILS:-tree} \
  --class-path "$CP$SEP$OUT/classes" "$@"
