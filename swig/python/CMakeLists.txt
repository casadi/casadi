cmake_minimum_required(VERSION 3.10.2)

# Workaround for OS X
if(APPLE)
  if(WITH_PYTHON3)
  set(PYTHON_CONFIG_NAME python3-config)
  else()
  set(PYTHON_CONFIG_NAME python-config)
  endif()
  find_program(PYTHON_CONFIG_EXECUTABLE
               NAMES ${PYTHON_CONFIG_NAME} DOC "python-config executable")
  if(PYTHON_CONFIG_EXECUTABLE)
    execute_process(COMMAND ${PYTHON_CONFIG_EXECUTABLE} --prefix
                    OUTPUT_VARIABLE PYTHON_PREFIX_STRING
                    RESULT_VARIABLE PYTHON_PREFIX_FAILED
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT PYTHON_PREFIX_FAILED)
      file(GLOB LOCAL_PYTHON_NAME RELATIVE ${PYTHON_PREFIX_STRING}/lib "${PYTHON_PREFIX_STRING}/lib/python*.*")
      find_library(PYTHON_LIBRARY NAMES ${LOCAL_PYTHON_NAME}
        PATHS ${PYTHON_PREFIX_STRING}/lib NO_DEFAULT_PATH)
      file(GLOB LOCAL_PYTHON_NAME "${PYTHON_PREFIX_STRING}/include/python*")
      find_path(PYTHON_INCLUDE_DIR
        NAMES Python.h
        PATHS ${LOCAL_PYTHON_NAME} NO_DEFAULT_PATH)
    endif()
    message(STATUS ${PYTHON_INCLUDE_DIR})
  endif()
endif()

function(CasadiGetPythonVersionMajor VERSION_MAJOR_RESULT INCLUDE_PATH)
  set(MATCH_STRING "#define[ \t]+PY_MAJOR_VERSION[ \t]+([0-9]+)")
  file(STRINGS "${INCLUDE_PATH}/patchlevel.h" PY_MAJOR_VERSION REGEX ${MATCH_STRING})
  string(REGEX REPLACE "${MATCH_STRING}" "\\1" PYTHON_VERSION_MAJOR "${PY_MAJOR_VERSION}")
  set(${VERSION_MAJOR_RESULT} ${PYTHON_VERSION_MAJOR} PARENT_SCOPE)
endfunction()

if(SWIG_EXPORT)
  if(WITH_PYTHON3)
    set(PYTHON_VERSION_MAJOR 3)
  else()
    set(PYTHON_VERSION_MAJOR 2)
  endif()
else()
  if(PYTHON_INCLUDE_DIR)
    set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_DIR})
  else()
    # Find packages
    if(WITH_PYTHON3)
    set(MINPYVERSION "3")
    else()
    set(MINPYVERSION "")
    endif()

    # We don't really need the libs, but the headers.
    find_package(PythonInterp ${MINPYVERSION})
    find_package(PythonLibs ${MINPYVERSION})
    if(NOT PYTHON_INCLUDE_PATH)
      message(SEND_ERROR "Python include headers not found, but required." )
    endif()
  endif()
  
  CasadiGetPythonVersionMajor(PYTHON_VERSION_MAJOR "${PYTHON_INCLUDE_PATH}")
  include_directories(${PYTHON_INCLUDE_PATH})
  
  
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)

if(WITH_PYTHON_INTERRUPTS)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()
if(WITH_COPYSIGN_UNDEF)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DCASADI_WITH_COPYSIGN_UNDEF")
endif()

set(PYTHONFLAG "")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DPy_USING_UNICODE")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-noproxydel")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-fastunpack")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-modernargs")
if(WITH_PYTHON_GIL_RELEASE)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-threads")
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DCASADI_WITH_PYTHON_GIL_RELEASE")
endif()
if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
set(PYTHONFLAG "3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-py3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON3")
endif()

message("PYTHONFLAGS: ${CMAKE_SWIG_FLAGS} ${PYTHONFLAG}")

if(SWIG_EXPORT)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra)
  set(CMAKE_SWIG_WRAP_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/source/dummy)
  
  configure_file(setup.py.in "${PROJECT_SOURCE_DIR}/setup.py" @ONLY)
  file(COPY requirements.txt DESTINATION ${PROJECT_SOURCE_DIR})
  file(COPY MANIFEST.in DESTINATION ${PROJECT_SOURCE_DIR})
endif()

if(SWIG_IMPORT)
  set(PYTHON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/source/casadiPYTHON_wrap.cxx)
  set_source_files_properties(${PYTHON_FILE}  PROPERTIES  CPLUSPLUS ON)
else()
  # Generate SWIG wrapper
  set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_module_initialize(casadi python)
  swig_add_source_to_module(casadi FALSE PYTHON_FILE ../casadi.i)
  get_directory_property(swig_extra_clean_files ADDITIONAL_MAKE_CLEAN_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${PYTHON_FILE}")
endif()

if(SWIG_EXPORT)
  add_custom_target(python_source DEPENDS ${PYTHON_FILE})
endif()

# Make target
add_custom_target(python DEPENDS _casadi)

add_library(_casadi MODULE ${PYTHON_FILE})

set_target_properties(_casadi PROPERTIES PREFIX "")
if(WIN32 AND NOT CYGWIN)
  set_target_properties(_casadi PROPERTIES SUFFIX ".pyd")
  set(CASADI_PYTHON_LIBRARY_SUFFIX ".pyd")
else()
  set(CASADI_PYTHON_LIBRARY_SUFFIX ${CMAKE_SHARED_MODULE_SUFFIX})
endif()

if(APPLE AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
  set_target_properties(_casadi PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR} -Wno-maybe-uninitialized")
  if(WIN32)
    # Fix for error: '_hypot' has not been declared in 'std'
    # cfr https://github.com/tudat-team/tudatpy/commit/ea9ef90a462825a65741f89ac1776413868f52cb
    set(COMPILE_FLAGS "${COMPILE_FLAGS} -D_hypot=hypot")
    # Fix for error: casadiPYTHON_wrap.cxx:(.text+0xd25cc): undefined reference to `__imp_Py_InitModule4'
    # cfr https://github.com/tpaviot/pythonocc-core/commit/9aef6c8e86c5b77c33742acf36e45db0b2f9ca23#diff-1e7de1ae2d059d21e1dd75d5812d5a34b0222cef273b7c3a2af62eb747f9d20a
    if(${CMAKE_SIZEOF_VOID_P} MATCHES "8")
      set(COMPILE_FLAGS "${COMPILE_FLAGS} -DMS_WIN64")
    endif()
  endif()
else()
  set(COMPILE_FLAGS)
endif()

if(WITH_PYTHON_GIL_RELEASE)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} -DCASADI_WITH_PYTHON_GIL_RELEASE")
endif()

message("COMPILE_FLAGS: ${COMPILE_FLAGS}")

set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")

if(WIN32) # dlls don't allow missing symbols
  if(NOT PYTHON_LIBRARY)
    message(SEND_ERROR "Python libraries not found, but required." )
  endif()
  target_link_libraries(_casadi ${PYTHON_LIBRARY})
endif()
target_link_libraries(_casadi casadi)

# Custom installation command for Python
add_custom_target(install_python
COMMAND ${CMAKE_COMMAND}
  -D COMPONENT=install_python
  -D CMAKE_INSTALL_PREFIX="${PYTHON_PREFIX}"
  -P cmake_install.cmake
)
add_dependencies(install_python _casadi)

# Install C++ wrapper library
install(TARGETS _casadi
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

if (SWIG_IMPORT)
# Install Python proxy classes
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)
else()
# Install Python proxy classes
install(FILES ${PROJECT_BINARY_DIR}/swig/python/casadi.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)
endif()

# Install Python tools
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
)

# Install Python package initialization
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION "${PYTHON_PREFIX}/casadi"
  COMPONENT install_python
)

# Install pip metadata files to ensure that casadi installed via CMake is listed by pip list
# See https://packaging.python.org/specifications/recording-installed-packages/
# and https://packaging.python.org/en/latest/specifications/core-metadata/#core-metadata
option(CASADI_PYTHON_PIP_METADATA_INSTALL "Use CMake to install Python pip metadata. Set to off if some other tool already installs it." OFF)
mark_as_advanced(CASADI_PYTHON_PIP_METADATA_INSTALL)
set(CASADI_PYTHON_PIP_METADATA_INSTALLER "cmake" CACHE STRING "Specify the string to identify the pip Installer. Default: cmake, change this if you are using another tool.")
mark_as_advanced(CASADI_PYTHON_PIP_METADATA_INSTALLER)
if (CASADI_PYTHON_PIP_METADATA_INSTALL)
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/METADATA "")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Metadata-Version: 2.1\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Name: casadi\n")
  file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/METADATA "Version: ${CASADI_VERSION}\n")
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/INSTALLER "${CASADI_PYTHON_PIP_METADATA_INSTALLER}\n")
  install(
    FILES "${CMAKE_CURRENT_BINARY_DIR}/METADATA" "${CMAKE_CURRENT_BINARY_DIR}/INSTALLER"
    DESTINATION ${PYTHON_PREFIX}/casadi-${CASADI_VERSION}.dist-info
    COMPONENT install_python)
endif()

# Example of how to extend CasADi with additional features
if (WITH_EXTENDING_CASADI)
  set_source_files_properties(../extending_casadi/extending_casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_add_module(extending_casadi python ../extending_casadi/extending_casadi.i)
  swig_link_libraries(extending_casadi extending_casadi)
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  endif()
  set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_extending_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

  add_custom_target(extending_casadi_python DEPENDS _extending_casadi extending_casadi)

  install(TARGETS _extending_casadi
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${PROJECT_BINARY_DIR}/swig/python/extending_casadi.py
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../extending_casadi/__init__.py
    DESTINATION "${PYTHON_PREFIX}/extending_casadi"
    COMPONENT install_python
  )

endif()
