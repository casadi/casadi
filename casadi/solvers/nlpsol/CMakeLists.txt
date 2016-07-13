cmake_minimum_required(VERSION 2.8.6)

# SQPMethod -  A basic SQP method
casadi_plugin(Nlpsol sqpmethod
  sqpmethod.hpp sqpmethod.cpp sqpmethod_meta.cpp)

# SCPgen -  An implementation of Lifted Newton SQP
casadi_plugin(Nlpsol scpgen
  scpgen.hpp scpgen.cpp scpgen_meta.cpp)
