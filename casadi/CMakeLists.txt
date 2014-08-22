
#define_property(TARGET PROPERTY PKG_CONFIG_REQUIRES
#  BRIEF_DOCS "add requires info for target"
#  FULL_DOCS "add requires info for target"
#)

set(LINT_TARGETS)
set(SPELL_TARGETS)

macro(casadi_plugin Type name)
  string(TOLOWER ${Type} type)
  casadi_library(casadi_${type}_${name} ${ARGN})
  set_property(GLOBAL APPEND PROPERTY CASADI_PLUGINS "${Type}::${name}")
endmacro()

macro(casadi_plugin_link_libraries Type name)
  string(TOLOWER ${Type} type)
  target_link_libraries(casadi_${type}_${name} ${ARGN})
endmacro()

macro(casadi_library name)
  if(ENABLE_STATIC)
    add_library(${name} STATIC ${ARGN})
  endif()
  if(ENABLE_SHARED)
    add_definitions(-DSHARED_LIBRARY)
    if(NOT ENABLE_EXPORT_ALL)
      add_compiler_export_flags()
    endif()
    add_library(${name} SHARED ${ARGN})
  endif()
  generate_export_header(${name})
  
  # Make sure the exported header is installed
  file(RELATIVE_PATH RELPATH "${CMAKE_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/${name}_export.h")
  get_filename_component(RELPATH ${RELPATH} PATH)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${name}_export.h
  DESTINATION include/${RELPATH})
  
  if(NOT "${name}" STREQUAL "casadi_core")
    target_link_libraries(${name} casadi_core)
  endif()
 
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
  add_style_check_target(${name} "${ARGN}")
  set(LINT_TARGETS ${LINT_TARGETS} "lint_${name}" PARENT_SCOPE)
  add_spell_check_target(${name} "${ARGN}")
  set(SPELL_TARGETS ${SPELL_TARGETS} "spell_${name}" PARENT_SCOPE)

  set_property(GLOBAL APPEND PROPERTY CASADI_MODULES ${name})
endmacro()

add_subdirectory(core)  # needed by all except external_packages
add_subdirectory(interfaces)             # needs external_packages
add_subdirectory(solvers)                # needs interfaces

# add a dummy custom target which does nothing but
# depend on the other lint targets
add_custom_target(lint
  DEPENDS ${LINT_TARGETS}
  VERBATIM)
  
add_custom_target(spell
  DEPENDS ${SPELL_TARGETS}
  VERBATIM)

get_property(CASADI_MODULES GLOBAL PROPERTY CASADI_MODULES)
get_property(CASADI_SWIG_MODULES GLOBAL PROPERTY CASADI_SWIG_MODULES)
get_property(CASADI_PLUGINS GLOBAL PROPERTY CASADI_PLUGINS)

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
set(CASADI_PLUGINS ${CASADI_PLUGINS} PARENT_SCOPE)
set(CASADI_SWIG_MODULES ${CASADI_SWIG_MODULES} PARENT_SCOPE)

foreach(name ${CASADI_MODULES})
  set(PC_NAME "${name}")
  #set(PC_CXX_FLAGS "{CXX_FLAGS}")
  
  function(join VALUES GLUE OUTPUT)
    string (REPLACE ";" "${GLUE}" _TMP_STR "${VALUES}")
    set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
  endfunction()
  
  function(obtain_interface_link_libraries TARGET OUTPUT)
    set(LIBLIST)
    get_target_property(TEMP_INTERFACE_LINK_LIBRARIES ${TARGET} "INTERFACE_LINK_LIBRARIES")
    if(TEMP_INTERFACE_LINK_LIBRARIES)
    foreach(LIB ${TEMP_INTERFACE_LINK_LIBRARIES})
      list(APPEND LIBLIST ${LIB})
      #if(TARGET ${LIB})
      #  obtain_interface_link_libraries(${LIB} OUTPUT2)
      #  if(OUTPUT2)
      #    list(APPEND LIBLIST ${OUTPUT2})
      #  endif()
      #endif()
    endforeach()
    endif()
    set(${OUTPUT} ${LIBLIST} PARENT_SCOPE)
  endfunction()
  
  obtain_interface_link_libraries(${name} OUTPUT)
  set(PC_REQUIRES)
  set(PC_LIBRARIES)
  
  foreach(LIB ${OUTPUT})
    if(TARGET ${LIB})
      list(FIND CASADI_MODULES ${LIB} FOUND)
      if (FOUND EQUAL -1 )
        list(APPEND PC_LIBRARIES ${LIB})
      else()
        list(APPEND PC_REQUIRES "${LIB}")
      endif()
    else()
      list(APPEND PC_LIBRARIES ${LIB})
    endif()
  endforeach()
  join("${PC_REQUIRES}" "," PC_REQUIRES)

  set(PC_LIBRARIES_LIST "-l${name}")
  foreach(LIB ${PC_LIBRARIES})
    string(REGEX MATCH "[\\/](lib)?([a-zA-Z0-9_]+)\\.[a-zA-Z0-9_]+$" LIBMATCH ${LIB})
    if(LIBMATCH)
      set(PC_LIBRARIES_LIST "${PC_LIBRARIES_LIST} -l${CMAKE_MATCH_2}")
    else()
      set(PC_LIBRARIES_LIST "${PC_LIBRARIES_LIST} -l${LIB}")
    endif()
  endforeach()
  get_target_property(PC_LINK_LIBRARIES ${name} "LINK_LIBRARIES")
  get_target_property(PC_PKG_CONFIG_REQUIRES ${name} "PKG_CONFIG_REQUIRES")
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/pkgconfig.pc.in
    ${CMAKE_BINARY_DIR}/lib/pkgconfig/${name}.pc
    ESCAPE_QUOTES
    @ONLY)
  install(FILES ${CMAKE_BINARY_DIR}/lib/pkgconfig/${name}.pc
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig
  )
endforeach()

install(FILES casadi.hpp DESTINATION include/casadi)
