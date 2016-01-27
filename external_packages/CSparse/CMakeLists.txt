cmake_minimum_required(VERSION 2.8.6)

set(CSPARSE_SRC csparse_mod.h csparse_mod.c)
casadi_external_library(casadi_csparse ${CSPARSE_SRC})
