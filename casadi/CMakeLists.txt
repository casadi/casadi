
set(LINT_TARGETS)

macro(casadi_library name)
if(ENABLE_STATIC)
add_library(${name} STATIC ${ARGN})
endif()
if(ENABLE_SHARED)
add_compiler_export_flags()
add_library(${name} SHARED ${ARGN})
generate_export_header(${name})
endif()
if(NOT "${name}" STREQUAL "casadi_symbolic")
target_link_libraries(${name} casadi_symbolic)
endif()
set(CASADI_MODULES ${CASADI_MODULES} ${name} PARENT_SCOPE)
install(TARGETS ${name}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
add_style_check_target(${name} "${ARGN}")
set(LINT_TARGETS ${LINT_TARGETS} "lint_${name}" PARENT_SCOPE)
endmacro()

add_subdirectory(symbolic)  # needed by all except external_packages
add_subdirectory(interfaces)             # needs external_packages
add_subdirectory(nonlinear_programming)  # needs interfaces
add_subdirectory(convex_programming)     # needs interfaces
add_subdirectory(integration)
add_subdirectory(optimal_control)        # needs integration, external_packages
add_subdirectory(control)

# add a dummy custom target which does nothing but
# depend on the other lint targets
add_custom_target(lint
  DEPENDS ${LINT_TARGETS}
  VERBATIM)


set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
