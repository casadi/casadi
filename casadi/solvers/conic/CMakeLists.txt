cmake_minimum_required(VERSION 2.8.6)

casadi_plugin(Conic nlpsol
  qp_to_nlp.hpp qp_to_nlp.cpp qp_to_nlp_meta.cpp)
