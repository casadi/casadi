
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
endmacro()

add_subdirectory(symbolic)  # needed by all except external_packages
add_subdirectory(interfaces)             # needs external_packages
add_subdirectory(nonlinear_programming)  # needs interfaces
add_subdirectory(convex_programming)     # needs interfaces
add_subdirectory(integration)
add_subdirectory(optimal_control)        # needs integration, external_packages
add_subdirectory(control)

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
