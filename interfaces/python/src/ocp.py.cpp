#include <alpaqa/inner/directions/panoc-ocp/ocp-vars.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;

template <alpaqa::Config Conf>
void register_ocp(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    m.def("eval_ocp",
          [](const alpaqa::TypeErasedControlProblem<config_t> &problem, crvec u, crvec y, crvec μ) {
              alpaqa::OCPEvaluator<config_t> eval{problem};
              auto N                    = eval.vars.N;
              auto nu                   = eval.vars.nu();
              auto nx                   = eval.vars.nx();
              auto nc                   = eval.vars.nc();
              auto nc_N                 = eval.vars.nc_N();
              alpaqa::Box<config_t> U   = alpaqa::Box<config_t>::NaN(nu);
              alpaqa::Box<config_t> D   = alpaqa::Box<config_t>::NaN(nc);
              alpaqa::Box<config_t> D_N = alpaqa::Box<config_t>::NaN(nc_N);
              problem.get_U(U);     // input box constraints
              problem.get_D(D);     // general constraints
              problem.get_D_N(D_N); // general terminal constraints
              assert((nc == 0 && nc_N == 0) || μ.size() == N + 1);
              assert(y.size() == nc * N + nc_N);
              vec storage = eval.vars.create();
              alpaqa::detail::assign_interleave_xu(eval.vars, u, storage);
              problem.get_x_init(eval.vars.xk(storage, 0));
              vec qr      = eval.vars.create_qr();
              auto V      = eval.forward(storage, D, D_N, μ, y);
              vec grad(N * nu);
              vec λ(nx);
              vec w(nx);
              vec v(std::max(nc, nc_N));
              auto mut_qrk = [&](index_t k) -> rvec { return eval.vars.qrk(qr, k); };
              auto mut_q_N = [&]() -> rvec { return eval.vars.qk(qr, N); };
              eval.backward(storage, grad, λ, w, v, mut_qrk, mut_q_N, D, D_N, μ, y);
              return std::tuple{V, std::move(grad)};
          });
}

template void register_ocp<alpaqa::EigenConfigd>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigf>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_ocp<alpaqa::EigenConfigq>(py::module_ &);
#endif
