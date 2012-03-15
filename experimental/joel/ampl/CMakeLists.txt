IF(IPOPT_FOUND)
add_executable(ipopt_casadi ipopt_casadi.cpp)
target_link_libraries(ipopt_casadi ipopt_interface nonlinear_programming casadi ${IPOPT_LIBRARIES} ${CASADI_DEPENDENCIES})
ENDIF(IPOPT_FOUND)
