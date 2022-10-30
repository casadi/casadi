cmake_minimum_required(VERSION 2.8.12)

# C runtime sources
set(RUNTIME_SRC
  casadi_axpy.hpp
  casadi_bilin.hpp
  casadi_copy.hpp
  casadi_cvx.hpp
  casadi_de_boor.hpp
  casadi_densify.hpp
  casadi_dot.hpp
  casadi_feasiblesqpmethod.hpp
  casadi_clear.hpp
  casadi_clip_min.hpp
  casadi_clip_max.hpp
  casadi_fill.hpp
  casadi_flip.hpp
  casadi_file_slurp.hpp
  casadi_getu.hpp
  casadi_iamax.hpp
  casadi_interpn.hpp
  casadi_interpn_grad.hpp
  casadi_interpn_interpolate.hpp
  casadi_interpn_weights.hpp
  casadi_kron.hpp
  casadi_low.hpp
  casadi_max_viol.hpp
  casadi_mmin.hpp
  casadi_mmax.hpp
  casadi_mtimes.hpp
  casadi_vfmin.hpp
  casadi_vfmax.hpp
  casadi_mv.hpp
  casadi_trilsolve.hpp
  casadi_triusolve.hpp
  casadi_mv_dense.hpp
  casadi_nd_boor_eval.hpp
  casadi_nd_boor_dual_eval.hpp
  casadi_norm_1.hpp
  casadi_norm_2.hpp
  casadi_norm_inf.hpp
  casadi_masked_norm_inf.hpp
  casadi_norm_inf_mul.hpp
  casadi_polyval.hpp
  casadi_project.hpp
  casadi_tri_project.hpp
  casadi_rank1.hpp
  casadi_scal.hpp
  casadi_sparsify.hpp
  casadi_sum_viol.hpp
  casadi_swap.hpp
  casadi_trans.hpp
  casadi_finite_diff.hpp
  casadi_ldl.hpp
  casadi_qr.hpp
  casadi_qp.hpp
  casadi_kkt.hpp
  casadi_ipqp.hpp
  casadi_nlp.hpp
  casadi_sqpmethod.hpp
  casadi_bfgs.hpp
  casadi_regularize.hpp
  casadi_newton.hpp
  casadi_bound_consistency.hpp
  casadi_lsqr.hpp
  casadi_cache.hpp
  casadi_convexify.hpp
  casadi_logsumexp.hpp
)
set(CASADI_RUNTIME_SRC "${RUNTIME_SRC}" PARENT_SCOPE)

set(EXTRA_CASADI_CXX_FLAGS "${EXTRA_CASADI_CXX_FLAGS} -Wno-conversion")

install(FILES casadi_runtime.hpp shared.hpp ${RUNTIME_SRC}
  DESTINATION ${INCLUDE_PREFIX}/casadi/core/runtime
)
