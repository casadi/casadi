macro(casadi_plugin_link_libraries Type name)
  string(TOLOWER ${Type} type)
  target_link_libraries(casadi_${type}_${name} ${ARGN})
endmacro()

macro(casadi_plugin Type name)
  string(TOLOWER ${Type} type)
  casadi_library_or_plugin(casadi_${type}_${name} TRUE ${ARGN})
  set_property(GLOBAL APPEND PROPERTY CASADI_PLUGINS "${Type}::${name}")
endmacro()

macro(casadi_library name)
  casadi_library_or_plugin(${name} FALSE ${ARGN})
endmacro()

if(ENABLE_SHARED)
  if(NOT ENABLE_EXPORT_ALL)
    add_compiler_export_flags()
  endif()
endif()

# this function handles libraries and plugins
# there is so much similarity in plugins/libraries, we just
# use a flag and wrap this with casadi_plugin or casadi_library (above)
macro(casadi_library_or_plugin name is_plugin)
  if(ENABLE_STATIC)
    add_library(${name} STATIC ${ARGN})
  endif()
  if(ENABLE_SHARED)
    add_definitions(-DSHARED_LIBRARY)
    add_library(${name} SHARED ${ARGN})
  endif()
  if(NOT "${name}" STREQUAL "casadi")
    generate_export_header(${name})

    # Make sure the exported header is installed
    file(RELATIVE_PATH RELPATH "${CMAKE_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${name}_export.h")
    get_filename_component(RELPATH ${RELPATH} PATH)

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

  add_style_check_target(${name} "${ARGN}")
  set_property(GLOBAL APPEND PROPERTY LINT_TARGETS  "lint_${name}")
  add_spell_check_target(${name} "${ARGN}")
  set_property(GLOBAL APPEND PROPERTY SPELL_TARGETS  "spell_${name}")
  set(SPELL_TARGETS ${SPELL_TARGETS} "spell_${name}" PARENT_SCOPE)

  set_property(GLOBAL APPEND PROPERTY CASADI_MODULES ${name})
endmacro()

add_definitions(-DSHARED_LIBRARY_PREFIX="\\"${CMAKE_SHARED_LIBRARY_PREFIX}\\"")
add_definitions(-DSHARED_LIBRARY_SUFFIX="\\"${CMAKE_SHARED_LIBRARY_SUFFIX}\\"")

add_subdirectory(core)        # needed by all except external_packages
add_subdirectory(interfaces)  # needs external_packages
add_subdirectory(solvers)     # needs interfaces

get_property(LINT_TARGETS GLOBAL PROPERTY LINT_TARGETS)
get_property(SPELL_TARGETS GLOBAL PROPERTY SPELL_TARGETS)

# add a dummy custom target which does nothing but
# depend on the other lint targets
add_custom_target(lint
  DEPENDS ${LINT_TARGETS}
  VERBATIM)

add_custom_target(spell
  DEPENDS ${SPELL_TARGETS}
  VERBATIM)

get_property(CASADI_MODULES GLOBAL PROPERTY CASADI_MODULES)
get_property(CASADI_PLUGINS GLOBAL PROPERTY CASADI_PLUGINS)

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
set(CASADI_PLUGINS ${CASADI_PLUGINS} PARENT_SCOPE)

add_custom_target(libs ALL DEPENDS ${CASADI_MODULES})

install(FILES casadi.hpp DESTINATION include/casadi)
