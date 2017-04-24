# Turn a runtime file into a file with strings
macro(CASADI_STRINGIFY STRFILE)
  # Start with an empty file
  file(WRITE ${STRFILE} "")
  foreach(FILE ${ARGN})
    # Add declaration of string
    get_filename_component(FILENAME ${FILE} NAME_WE)
    file(APPEND ${STRFILE} "const char* ${FILENAME}_str =")
    # Append file as strings
    file(STRINGS ${FILE} FILE_CONTENTS)
    foreach(LINE ${FILE_CONTENTS})
      string(REPLACE "\\" "\\\\" LINE ${LINE}) # Replace \ with \\
      string(REPLACE "\"" "\\\"" LINE ${LINE}) # Replace " with \"
      file(APPEND ${STRFILE} "\n  \"${LINE}\\n\"")
    endforeach()
    # End declaration
    file(APPEND ${STRFILE} ";\n\n")
  endforeach()
endmacro()

# this function handles libraries and plugins
# there is so much similarity in plugins/libraries, we just
# use a flag and wrap this with casadi_plugin or casadi_library (above)
macro(casadi_library name)
  if(ENABLE_STATIC)
    add_library(${name} STATIC ${ARGN})
  endif()
  if(ENABLE_SHARED)
    add_definitions(-DSHARED_LIBRARY)
    add_library(${name} SHARED ${ARGN})
  endif()
  generate_export_header(${name})

  # Link with the core
  if(NOT "${name}" STREQUAL "casadi")
    target_link_libraries(${name} casadi)
  endif()

  if(WITH_SELFCONTAINED)
    if (WITH_PYTHON)
      install(TARGETS ${name}
        LIBRARY DESTINATION casadi
        ARCHIVE DESTINATION casadi
        RUNTIME DESTINATION casadi
      )
    endif()
    if (WITH_MATLAB OR WITH_OCTAVE)
      install(TARGETS ${name}
        LIBRARY DESTINATION .
        ARCHIVE DESTINATION .
        RUNTIME DESTINATION .
      )
    endif()
  else()
    install(TARGETS ${name}
      LIBRARY DESTINATION lib
      ARCHIVE DESTINATION lib
      RUNTIME DESTINATION lib
    )
    SET_TARGET_PROPERTIES(
      ${name}
    PROPERTIES
      SOVERSION ${PACKAGE_SO_VERSION}
    )
  endif()

  if(WITH_LINT)
    add_style_check_target(${name} "${ARGN}")
    set_property(GLOBAL APPEND PROPERTY LINT_TARGETS  "lint_${name}")
  endif()

  if(WITH_SPELL)
    add_spell_check_target(${name} "${ARGN}")
    set_property(GLOBAL APPEND PROPERTY SPELL_TARGETS  "spell_${name}")
    set(SPELL_TARGETS ${SPELL_TARGETS} "spell_${name}" PARENT_SCOPE)
  endif()

  set_property(GLOBAL APPEND PROPERTY CASADI_MODULES ${name})
endmacro()

macro(casadi_plugin Type name)
  string(TOLOWER ${Type} type)
  casadi_library(casadi_${type}_${name} ${ARGN})
  set_property(GLOBAL APPEND PROPERTY CASADI_PLUGINS "${Type}::${name}")
endmacro()

macro(casadi_plugin_link_libraries Type name)
  string(TOLOWER ${Type} type)
  target_link_libraries(casadi_${type}_${name} ${ARGN})
endmacro()

if(ENABLE_SHARED)
  if(NOT ENABLE_EXPORT_ALL)
    add_compiler_export_flags()
  endif()
endif()

add_definitions(-DSHARED_LIBRARY_PREFIX="\\"${CMAKE_SHARED_LIBRARY_PREFIX}\\"")
add_definitions(-DSHARED_LIBRARY_SUFFIX="\\"${CMAKE_SHARED_LIBRARY_SUFFIX}\\"")

add_subdirectory(core)
add_subdirectory(interfaces)
add_subdirectory(solvers)

if(WITH_LINT)
  get_property(LINT_TARGETS GLOBAL PROPERTY LINT_TARGETS)
  add_custom_target(lint DEPENDS ${LINT_TARGETS} VERBATIM)
endif()
if(WITH_SPELL)
  get_property(SPELL_TARGETS GLOBAL PROPERTY SPELL_TARGETS)
  add_custom_target(spell DEPENDS ${SPELL_TARGETS} VERBATIM)
endif()

get_property(CASADI_MODULES GLOBAL PROPERTY CASADI_MODULES)
get_property(CASADI_PLUGINS GLOBAL PROPERTY CASADI_PLUGINS)

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
set(CASADI_PLUGINS ${CASADI_PLUGINS} PARENT_SCOPE)

add_custom_target(libs ALL DEPENDS ${CASADI_MODULES})

install(FILES casadi.hpp DESTINATION include/casadi)
