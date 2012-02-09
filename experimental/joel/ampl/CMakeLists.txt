IF(IPOPT_FOUND)

add_executable(asl_reader asl_reader.cpp asl_internal.cpp asl_reader.hpp )
# target_link_libraries(nlc casadi)
include_directories(/home/janderss/src/ampl/netlib.org/ampl/solvers)
target_link_libraries(asl_reader /home/janderss/src/ampl/netlib.org/ampl/solvers/funcadd0.o /home/janderss/src/ampl/netlib.org/ampl/solvers/amplsolver.a ipopt_interface casadi ${IPOPT_LIBRARIES} m)


ENDIF(IPOPT_FOUND)
