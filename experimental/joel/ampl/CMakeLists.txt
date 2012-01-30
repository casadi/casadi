add_executable(ampl_reader ampl_reader.cpp ampl_reader_cops.cpp ampl_reader.hpp )
# target_link_libraries(nlc casadi)
include_directories(/home/janderss/src/ampl/netlib.org/ampl/solvers)
target_link_libraries(ampl_reader /home/janderss/src/ampl/netlib.org/ampl/solvers/funcadd0.o /home/janderss/src/ampl/netlib.org/ampl/solvers/amplsolver.a casadi m)
