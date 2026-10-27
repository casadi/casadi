cmake_minimum_required(VERSION 3.10.2)

# Three-stage build, mirroring swig/julia:
#   SWIG_EXPORT  -- run SWIG only; emit casadiJAVA_wrap.cxx + java/org/casadi/*.java
#                   into target/source (no JDK needed).
#   SWIG_IMPORT  -- compile the pre-generated wrapper; SWIG not required.
#   neither      -- local all-in-one: run SWIG inline, then build.
# The native library needs only the JNI headers; casadi.jar is built when a JDK is found.

set(_iface    ${PROJECT_SOURCE_DIR}/swig/casadi.i)
set(_wrap_dir ${CMAKE_CURRENT_SOURCE_DIR}/target/source)

if(SWIG_IMPORT)
  message(STATUS "SWIG_IMPORT: compiling pre-generated java wrapper; SWIG not invoked")
  set(_out_dir ${_wrap_dir})
else()
  find_package(SWIG REQUIRED)
  execute_process(COMMAND ${SWIG_EXECUTABLE} -swiglib OUTPUT_VARIABLE _swig_lib OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT EXISTS "${_swig_lib}/java/javaguest.swg")
    message(FATAL_ERROR
      "SWIG at ${SWIG_EXECUTABLE} lacks the Java guest-object mode (java/javaguest.swg). "
      "Use the casadi-patched SWIG (ghcr.io/casadi/ci-swig).")
  endif()
  if(SWIG_EXPORT)
    set(_out_dir ${_wrap_dir})
  else()
    set(_out_dir ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()

set(_swig_cpp  ${_out_dir}/casadiJAVA_wrap.cxx)
set(_java_dir  ${_out_dir}/java)
set(_swig_stamp ${_out_dir}/java.stamp)

if(NOT SWIG_IMPORT)
  configure_file(../swig_config.h.cmake "${CMAKE_CURRENT_BINARY_DIR}/swig_config.h" ESCAPE_QUOTES)
  set(_swig_defs)
  if(WITH_DEPRECATED_FEATURES)
    list(APPEND _swig_defs -DWITH_DEPRECATED_FEATURES)
  endif()
  add_custom_command(
    OUTPUT  ${_swig_cpp} ${_swig_stamp}
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${_java_dir}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${_java_dir}/org/casadi
    COMMAND ${SWIG_EXECUTABLE} -java -c++ -module ca -package org.casadi ${_swig_defs}
            -I${PROJECT_SOURCE_DIR} -I${PROJECT_SOURCE_DIR}/swig -I${CMAKE_CURRENT_BINARY_DIR}
            -outdir ${_java_dir}/org/casadi -o ${_swig_cpp} ${_iface}
    COMMAND ${CMAKE_COMMAND} -E touch ${_swig_stamp}
    DEPENDS ${_iface}
    COMMENT "swig -java -> casadiJAVA_wrap.cxx + org/casadi/*.java"
    VERBATIM
  )
  add_custom_target(java_source DEPENDS ${_swig_cpp} ${_swig_stamp})
  if(SWIG_EXPORT)
    return()
  endif()
endif()

# ---- native library: JNI headers only, libjvm is not linked ----
# Cross builds pass -DJAVA_INCLUDE_PATH=<jdk>/include -DJAVA_INCLUDE_PATH2=<jdk>/include/<os>
if(NOT JAVA_INCLUDE_PATH OR NOT JAVA_INCLUDE_PATH2)
  set(JAVA_AWT_LIBRARY NotNeeded)
  set(JAVA_JVM_LIBRARY NotNeeded)
  set(JAVA_AWT_INCLUDE_PATH NotNeeded)
  find_package(JNI REQUIRED)
endif()
set(_jni_include ${JAVA_INCLUDE_PATH} ${JAVA_INCLUDE_PATH2})
message(STATUS "JNI headers: ${_jni_include}")

add_library(casadi_java SHARED ${_swig_cpp})
if(NOT SWIG_IMPORT)
  add_dependencies(casadi_java java_source)
endif()
if(WIN32)
  set_target_properties(casadi_java PROPERTIES PREFIX "" CXX_STANDARD 11)
else()
  set_target_properties(casadi_java PROPERTIES PREFIX "lib" CXX_STANDARD 11)
endif()
target_include_directories(casadi_java PRIVATE ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${_jni_include})
target_compile_definitions(casadi_java PRIVATE WITH_DEPRECATED_FEATURES)
if(MSVC)
  target_compile_options(casadi_java PRIVATE "/bigobj")
elseif(MINGW)
  target_compile_options(casadi_java PRIVATE "-Wa,-mbig-obj")
endif()
target_link_libraries(casadi_java casadi)
install(TARGETS casadi_java DESTINATION ${LIB_PREFIX})

# ---- casadi.jar: generated proxies + hand-written classes ----
find_package(Java 11 COMPONENTS Development)
if(Java_JAVAC_EXECUTABLE AND Java_JAR_EXECUTABLE)
  set(_jar ${CMAKE_CURRENT_BINARY_DIR}/casadi.jar)
  file(GLOB _java_extra ${CMAKE_CURRENT_SOURCE_DIR}/*.java)
  set(_java_deps ${_java_extra})
  if(NOT SWIG_IMPORT)
    list(APPEND _java_deps ${_swig_stamp})
  endif()
  add_custom_command(
    OUTPUT ${_jar}
    COMMAND ${CMAKE_COMMAND} -DJAVAC=${Java_JAVAC_EXECUTABLE} -DJAR=${Java_JAR_EXECUTABLE}
            -DGEN=${_java_dir} -DEXTRA=${CMAKE_CURRENT_SOURCE_DIR}
            -DWORK=${CMAKE_CURRENT_BINARY_DIR}/classes -DOUT=${_jar}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/jar.cmake
    DEPENDS ${_java_deps} ${CMAKE_CURRENT_SOURCE_DIR}/jar.cmake
    COMMENT "javac + jar -> casadi.jar"
    VERBATIM
  )
  add_custom_target(casadi_jar ALL DEPENDS ${_jar})
  if(NOT SWIG_IMPORT)
    add_dependencies(casadi_jar java_source)
  endif()
  install(FILES ${_jar} DESTINATION ${LIB_PREFIX})
else()
  message(STATUS "No JDK found: casadi.jar is not built (only libcasadi_java)")
endif()
