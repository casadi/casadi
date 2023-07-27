cmake_minimum_required(VERSION 3.10.2)

set(NLPSOL_ALPAQA_SRCS
  alpaqa_interface.cpp
  alpaqa_interface.hpp
  alpaqa_problem.cpp
  alpaqa_problem.hpp
)

casadi_plugin(Nlpsol alpaqa ${NLPSOL_ALPAQA_SRCS})
casadi_plugin_link_libraries(Nlpsol alpaqa alpaqa)
