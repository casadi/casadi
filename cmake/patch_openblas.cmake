file(READ CMakeLists.txt INPUT_TEXT)

string(REPLACE "libopenblas" "lib${CASADI_TP_SHARED_PREFIX}openblas" INPUT_TEXT "${INPUT_TEXT}")
string(REPLACE "-fpic -shared" "-fpic -L. -shared" INPUT_TEXT "${INPUT_TEXT}")
string(REPLACE "-Wl,-all_load" "-Wl,-all_load -Wl,-L." INPUT_TEXT "${INPUT_TEXT}")
string(REPLACE "add_subdirectory(utest)" "" INPUT_TEXT "${INPUT_TEXT}")

file(WRITE CMakeLists.txt "${INPUT_TEXT}")
file(APPEND CMakeLists.txt "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES LIBRARY_OUTPUT_NAME_DEBUG \"\${OpenBLAS_LIBNAME}\")\n")
file(APPEND CMakeLists.txt "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES PREFIX \"\${CMAKE_SHARED_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}\")\n")
file(APPEND CMakeLists.txt "set_target_properties( \${OpenBLAS_LIBS} PROPERTIES IMPORT_PREFIX \"\${CMAKE_IMPORT_LIBRARY_PREFIX}${CASADI_TP_SHARED_PREFIX}\")\n")
file(APPEND CMakeLists.txt "configure_file(\${PROJECT_SOURCE_DIR}/cmake/openblas.pc.in \${PROJECT_BINARY_DIR}/blas\${SUFFIX64}.pc @ONLY)\n")
file(APPEND CMakeLists.txt "install (FILES \${PROJECT_BINARY_DIR}/blas\${SUFFIX64}.pc DESTINATION \${CMAKE_INSTALL_LIBDIR}/pkgconfig/)\n")
file(APPEND CMakeLists.txt "configure_file(\${PROJECT_SOURCE_DIR}/cmake/openblas.pc.in \${PROJECT_BINARY_DIR}/lapack\${SUFFIX64}.pc @ONLY)\n")
file(APPEND CMakeLists.txt "install (FILES \${PROJECT_BINARY_DIR}/lapack\${SUFFIX64}.pc DESTINATION \${CMAKE_INSTALL_LIBDIR}/pkgconfig/)\n")

file(READ cmake/openblas.pc.in INPUT_TEXT)
string(REPLACE "lopenblas" "l${CASADI_TP_SHARED_PREFIX}openblas" INPUT_TEXT "${INPUT_TEXT}")
file(WRITE cmake/openblas.pc.in "${INPUT_TEXT}")
