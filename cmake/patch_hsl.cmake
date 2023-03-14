file(READ ${SRC}/libcoinhsl.sym INPUT_TEXT)


string(REGEX REPLACE "(m[ac][0-9][0-9][a-z])s_" "\\1_" INPUT_TEXT "${INPUT_TEXT}")

file(WRITE ${SRC}/libcoinhsl.sym "${INPUT_TEXT}")

file(READ ${SRC}/Makefile.in INPUT_TEXT)
string(REGEX REPLACE "(m__append_[0-9][0-9] =) -no-undefined" "\\1 -avoid-version -no-undefined" INPUT_TEXT "${INPUT_TEXT}")

file(WRITE ${SRC}/Makefile.in "${INPUT_TEXT}")
