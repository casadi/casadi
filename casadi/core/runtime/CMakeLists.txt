cmake_minimum_required(VERSION 2.8.6)

# C runtime sources
set(RUNTIME_SRC
  casadi_copy.hpp
  casadi_swap.hpp
  casadi_project.hpp
  casadi_densify.hpp
  casadi_sparsify.hpp
  casadi_scal.hpp
  casadi_axpy.hpp
  casadi_dot.hpp
  casadi_fill.hpp
  casadi_mtimes.hpp
  casadi_mv.hpp
  casadi_trans.hpp
  casadi_bilin.hpp
  casadi_rank1.hpp
  casadi_iamax.hpp
  casadi_interpn.hpp
  casadi_interpn_grad.hpp
  casadi_interpn_interpolate.hpp
  casadi_interpn_weights.hpp
  casadi_norm_1.hpp
  casadi_norm_2.hpp
  casadi_norm_inf.hpp
  casadi_de_boor.hpp
  casadi_nd_boor_eval.hpp
  casadi_flip.hpp
  casadi_low.hpp
)
set(CASADI_RUNTIME_SRC "${RUNTIME_SRC}" PARENT_SCOPE)

# Stringify C runtime
CASADI_STRINGIFY("${CMAKE_CURRENT_BINARY_DIR}/casadi_runtime_str.h"
  ${RUNTIME_SRC}
)

install(FILES casadi_runtime.hpp shared.hpp ${RUNTIME_SRC}
  DESTINATION include/casadi/core/runtime
)
