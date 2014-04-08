add_subdirectory(symbolic)  # needed by all except external_packages
add_subdirectory(interfaces)             # needs external_packages
add_subdirectory(nonlinear_programming)  # needs interfaces
add_subdirectory(convex_programming)     # needs interfaces
add_subdirectory(integration)
add_subdirectory(optimal_control)        # needs integration, external_packages
add_subdirectory(control)

set(CASADI_MODULES ${CASADI_MODULES} PARENT_SCOPE)
