cmake_minimum_required(VERSION 3.12)

# Find Python
find_package(Python REQUIRED COMPONENTS Interpreter Development)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Python_INCLUDE_DIRS})

# a python library is built in the build directory inside swig/python
make_directory(${PROJECT_BINARY_DIR}/python/casadi)

if(WITH_PYTHON_INTERRUPTS)
  set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON_INTERRUPTS")
endif()

set(PYTHONFLAG "")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DPy_USING_UNICODE")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-noproxydel")
if (${Python_VERSION_MAJOR} EQUAL 3)
set(PYTHONFLAG "3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-py3")
set(CMAKE_SWIG_FLAGS ${CMAKE_SWIG_FLAGS} "-DWITH_PYTHON3")
endif()

if(SWIG_EXPORT)
  set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra)
  set(CMAKE_SWIG_WRAP_OUTDIR ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/source/dummy)
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

if(APPLE AND ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
  set_target_properties(_casadi PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR}")
endif()

if(WIN32) # dlls don't allow missing symbols
  if(NOT Python_LIBRARIES)
    message(SEND_ERROR "Python libraries not found, but required." )
  endif()
  target_link_libraries(_casadi ${Python_LIBRARIES})
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
  DESTINATION "${Python_SITELIB}/casadi"
  COMPONENT install_python
)

if (SWIG_IMPORT)
# Install Python proxy classes
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/target${PYTHONFLAG}/extra/casadi.py
  DESTINATION "${Python_SITELIB}/casadi"
  COMPONENT install_python
)
else()
# Install Python proxy classes
install(FILES ${PROJECT_BINARY_DIR}/swig/python/casadi.py
  DESTINATION "${Python_SITELIB}/casadi"
  COMPONENT install_python
)
endif()

# Install Python tools
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools
  DESTINATION "${Python_SITELIB}/casadi"
  COMPONENT install_python
  USE_SOURCE_PERMISSIONS
  PATTERN .pyc EXCLUDE
  PATTERN .svn EXCLUDE
)

# Install Python package initialization
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION "${Python_SITELIB}/casadi"
  COMPONENT install_python
)

# Example of how to extend CasADi with additional features
if (WITH_EXTENDING_CASADI)
  set_source_files_properties(../extending_casadi/extending_casadi.i  PROPERTIES  CPLUSPLUS ON)
  swig_add_module(extending_casadi python ../extending_casadi/extending_casadi.i)
  swig_link_libraries(extending_casadi extending_casadi)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set_target_properties(_extending_casadi PROPERTIES COMPILE_FLAGS "-Wno-dynamic-class-memaccess -Wno-self-assign -Wno-maybe-uninitialized ${MAYBE_WERROR} -fno-builtin-strdup -fno-builtin-strndup")
  endif()
  set(CASADI_PYTHON_LIBRARY ${SWIG_MODULE_extending_casadi_REAL_NAME}${CASADI_PYTHON_LIBRARY_SUFFIX})

  add_custom_target(extending_casadi_python DEPENDS _extending_casadi extending_casadi)

  install(TARGETS _extending_casadi
    DESTINATION "${Python_SITELIB}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${PROJECT_BINARY_DIR}/swig/python/extending_casadi.py
    DESTINATION "${Python_SITELIB}/extending_casadi"
    COMPONENT install_python
  )

  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/../extending_casadi/__init__.py
    DESTINATION "${Python_SITELIB}/extending_casadi"
    COMPONENT install_python
  )

endif()
