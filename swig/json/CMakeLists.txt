cmake_minimum_required(VERSION 2.8.6)

# JSON/Haskell front-end
add_custom_target(json DEPENDS casadi_json)

# fooXML_wrap.xml are a big files which Joris wants to delete immedately after using.
# These files are dependencies of json outputs, but they don't exist according to cmake.
# So we silence these warnings. THIS IS A HACK AND WE NEED A BETTER SOLUTION
if(POLICY CMP0046)
  cmake_policy(SET CMP0046 OLD)
endif()

# Generate SWIG wrapper
set_source_files_properties(../casadi.i  PROPERTIES  CPLUSPLUS ON)
swig_add_module_xml(casadi_xml ../casadi.i)
add_custom_target(casadi_json
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/generate.py ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml casadi
  COMMAND ${CMAKE_COMMAND} -E remove ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml
)
add_dependencies(casadi_json casadi_xml
                 ${CMAKE_CURRENT_SOURCE_DIR}/generate.py
                 ${PROJECT_BINARY_DIR}/swig/casadiXML_wrap.xml)

#install(FILES ${PROJECT_BINARY_DIR}/swig/casadi.json
#  DESTINATION share/casadi
#  COMPONENT componentjson
#)

#set(CPACK_MONOLITHIC_INSTALL FALSE)
set(CPACK_TGZ_COMPONENT_INSTALL TRUE)

set(CPACK_COMPONENT_COMPONENTJSON_DISPLAY_NAME "Json meta-data")
set(CPACK_COMPONENT_COMPONENTJSON_DESCRIPTION
  "Installs meta-data about CasADi in json form.")

set(CPACK_COMPONENTS_GROUPING "ALL_COMPONENTS_IN_ONE")
set(CPACK_COMPONENTS_ALL Unspecified componentjson)
# This will not actually work since we are in a EXCLUDE_FROM_ALL
# http://cmake.3232098.n2.nabble.com/Yet-another-install-vs-exclude-from-all-thread-td4967491.html

# Pass meta information
configure_file(
  ${CMAKE_SOURCE_DIR}/swig/linklist.txt.cmake
  ${PROJECT_BINARY_DIR}/swig/linklist.txt
  ESCAPE_QUOTES
)

include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_CONTACT "${PACKAGE_MAINTAINER} <contact@casadi.org>")
set(CPACK_PACKAGE_NAME "libcasadi")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PACKAGE_DESCRIPTION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../../LICENSE.txt")
set(CPACK_DEBIAN_PACKAGE_LICENSE "LGPL")

set(THREE_PART_VERSION_REGEX "[0-9]+\\.[0-9]+\\.[0-9]+")

# Breaks up a string in the form n1.n2.n3 into three parts and stores
# them in major, minor, and patch.  version should be a value, not a
# variable, while major, minor and patch should be variables.
macro(THREE_PART_VERSION_TO_VARS version major minor patch)
  if(${version} MATCHES ${THREE_PART_VERSION_REGEX})
    string(REGEX REPLACE "^([0-9]+)\\.[0-9]+\\.[0-9]+\\+?" "\\1" ${major} "${version}")
    string(REGEX REPLACE "^[0-9]+\\.([0-9])+\\.[0-9]+\\+?" "\\1" ${minor} "${version}")
    string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+)\\+?" "\\1" ${patch} "${version}")
  else()
    message("macro(THREE_PART_VERSION_TO_VARS ${version} ${major} ${minor} ${patch}")
    message(FATAL_ERROR "Problem parsing version string, I can't parse it properly.")
  endif()
endmacro()

THREE_PART_VERSION_TO_VARS(${PACKAGE_VERSION} CPACK_PACKAGE_VERSION_MAJOR CPACK_PACKAGE_VERSION_MINOR CPACK_PACKAGE_VERSION_PATCH)

if(ENABLE_STATIC)
  set(CPACK_PACKAGE_FILE_NAME "libcasadi-static-${PACKAGE_VERSION_FULL}-${CMAKE_SYSTEM_PROCESSOR}")
else()
  set(CPACK_PACKAGE_FILE_NAME "libcasadi-shared-${PACKAGE_VERSION_FULL}-${CMAKE_SYSTEM_PROCESSOR}")
endif()

# if version end on a plus, it is a late build
if(${PACKAGE_VERSION} MATCHES "\\+$")
set(CPACK_PACKAGE_VERSION_PATCH "99")
endif()

include(CPack)
