add_executable(nlc nlc.c cops.c nlc.h )
# target_link_libraries(nlc casadi)
include_directories(/home/janderss/src/ampl/netlib.org/ampl/solvers)
target_link_libraries(nlc /home/janderss/src/ampl/netlib.org/ampl/solvers/funcadd0.o /home/janderss/src/ampl/netlib.org/ampl/solvers/amplsolver.a casadi m)
