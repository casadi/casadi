file(READ ${FILE} INPUT_TEXT)

string(REPLACE "libopenblas" "lib${CASADI_TP_SHARED_PREFIX}openblas" INPUT_TEXT "${INPUT_TEXT}")
string(REPLACE "-fpic -shared" "-fpic -L. -shared" INPUT_TEXT "${INPUT_TEXT}")
string(REPLACE "-Wl,-all_load" "-Wl,-all_load -Wl,-L." INPUT_TEXT "${INPUT_TEXT}")

file(WRITE ${FILE} "${INPUT_TEXT}")
file(APPEND ${FILE} "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES LIBRARY_OUTPUT_NAME_DEBUG \"\${OpenBLAS_LIBNAME}\")\n")
file(APPEND ${FILE} "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES PREFIX \"\${CMAKE_SHARED_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}\")\n")
file(APPEND ${FILE} "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES IMPORT_PREFIX \"\${CMAKE_IMPORT_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}\")\n")

