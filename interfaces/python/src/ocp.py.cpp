#include <alpaqa/inner/directions/panoc-ocp/lqr.hpp>
#include <alpaqa/inner/directions/panoc-ocp/ocp-vars.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/check-dim.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>
#include <alpaqa/util/index-set.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

template <alpaqa::Config Conf>
void register_ocp(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    struct OCPEvaluator {
        using Problem = alpaqa::TypeErasedControlProblem<config_t>;
        using OCPVars = alpaqa::OCPVariables<config_t>;
        alpaqa::util::copyable_unique_ptr<Problem> problem;
        alpaqa::OCPEvaluator<config_t> eval;
        alpaqa::Box<config_t> U{alpaqa::Box<config_t>::NaN(eval.vars.nu())};
        alpaqa::Box<config_t> D{alpaqa::Box<config_t>::NaN(eval.vars.nc())};
        alpaqa::Box<config_t> D_N{alpaqa::Box<config_t>::NaN(eval.vars.nc_N())};
        OCPEvaluator(const Problem &p) : problem{std::make_unique<Problem>(p)}, eval{*problem} {
            eval.problem->get_U(U);
            eval.problem->get_D(D);
            eval.problem->get_D_N(D_N);
        }

        auto prepare_y_μ(std::optional<vec> &&y, std::optional<vec> &&μ) const {
            auto N    = eval.vars.N;
            auto nu   = eval.vars.nu();
            auto nc   = eval.vars.nc();
            auto nc_N = eval.vars.nc_N();
            if (y)
                alpaqa::util::check_dim<config_t>("y", *y, nc * N + nc_N);
            else if (nc * N + nc_N == 0)
                y = vec{};
            else
                throw std::invalid_argument("Missing argument y");
            if (μ)
                alpaqa::util::check_dim<config_t>("μ", *μ, nc * N + nc_N);
            else if (nc * N + nc_N == 0)
                μ = vec{};
            else
                throw std::invalid_argument("Missing argument μ");
            return std::tuple{std::move(*y), std::move(*μ)};
        }

        auto prepare_storage(crvec u) const {
            auto N      = eval.vars.N;
            auto nu     = eval.vars.nu();
            auto nc     = eval.vars.nc();
            auto nc_N   = eval.vars.nc_N();
            vec storage = eval.vars.create();
            alpaqa::detail::assign_interleave_xu(eval.vars, u, storage);
            eval.problem->get_x_init(eval.vars.xk(storage, 0));
            return storage;
        }

        auto forward_backward(crvec u, std::optional<vec> y_, std::optional<vec> μ_) const {
            auto N        = eval.vars.N;
            auto nu       = eval.vars.nu();
            auto &&[y, μ] = prepare_y_μ(std::move(y_), std::move(μ_));
            vec storage   = prepare_storage(u);
            vec qr        = eval.vars.create_qr();
            vec grad(N * nu);
            auto mut_qrk = [&](index_t k) -> rvec { return eval.vars.qrk(qr, k); };
            auto mut_q_N = [&]() -> rvec { return eval.vars.qk(qr, N); };
            auto V       = eval.forward(storage, D, D_N, μ, y);
            eval.backward(storage, grad, mut_qrk, mut_q_N, D, D_N, μ, y);
            return std::tuple{V, std::move(grad)};
        }

        auto Qk(index_t k, crvec u, std::optional<vec> y_, std::optional<vec> μ_) const {
            auto nx       = eval.vars.nx();
            mat out       = mat::Zero(nx, nx);
            auto &&[y, μ] = prepare_y_μ(std::move(y_), std::move(μ_));
            vec storage   = prepare_storage(u);
            eval.forward_simulate(storage);
            eval.Qk(storage, y, μ, D, D_N, k, out);
            return out;
        }

        auto Rk(index_t k, crvec u, crindexvec mask) {
            mat out     = mat::Zero(mask.size(), mask.size());
            vec storage = prepare_storage(u);
            eval.forward_simulate(storage);
            eval.Rk(storage, k, mask, out);
            return out;
        }

        auto Sk(index_t k, crvec u, crindexvec mask) {
            auto nx     = eval.vars.nx();
            mat out     = mat::Zero(mask.size(), nx);
            vec storage = prepare_storage(u);
            eval.forward_simulate(storage);
            eval.Sk(storage, k, mask, out);
            return out;
        }

        auto active_indices(crvec u, crvec grad_ψ, real_t γ, rvec q) {
            auto N                  = eval.vars.N;
            auto nu                 = eval.vars.nu();
            auto is_constr_inactive = [&](index_t t, index_t i) {
                real_t ui = u(nu * t + i);
                // Gradient descent step.
                real_t gs = ui - γ * grad_ψ(t * nu + i);
                // Check whether the box constraints are active for this index.
                bool active_lb = gs <= U.lowerbound(i);
                bool active_ub = gs >= U.upperbound(i);
                if (active_ub) {
                    q(nu * t + i) = U.upperbound(i) - ui;
                    return false;
                } else if (active_lb) {
                    q(nu * t + i) = U.lowerbound(i) - ui;
                    return false;
                } else { // Store inactive indices
                    return true;
                }
            };
            alpaqa::detail::IndexSet<config_t> J{N, nu};
            J.update(is_constr_inactive);
            return J;
        }

        auto lqr_factor_solve(crvec u, real_t γ, std::optional<vec> y_, std::optional<vec> μ_) {
            auto &vars    = eval.vars;
            auto N        = vars.N;
            auto nu       = vars.nu();
            auto nx       = vars.nx();
            auto &&[y, μ] = prepare_y_μ(std::move(y_), std::move(μ_));
            // Compute x, h, c
            vec storage = prepare_storage(u);
            eval.forward_simulate(storage);
            // Compute gradients
            vec grad(N * nu);
            vec qr = vars.create_qr();
            eval.backward(storage, grad, vars.qr_mut(qr), vars.qN_mut(qr), D, D_N, μ, y);
            // Find active indices
            vec q(N * nu);
            auto J = active_indices(u, grad, γ, q);
            // Compute dynamics Jacobians
            mat jacs = vars.create_AB();
            for (index_t t = 0; t < N; ++t)
                problem->eval_jac_f(t, vars.xk(storage, t), vars.uk(storage, t), vars.ABk(jacs, t));
            // LQR factor
            using LQRFactor = alpaqa::StatefulLQRFactor<config_t>;
            LQRFactor lqr{{.N = N, .nx = nx, .nu = nu}};
            bool use_cholesky = false;
            auto uk_eq        = [&](index_t k) -> crvec { return q.segment(k * nu, nu); };
            auto Jk           = [&](index_t k) -> crindexvec { return J.indices(k); };
            auto Kk           = [&](index_t k) -> crindexvec { return J.compl_indices(k); };
            lqr.factor_masked(vars.AB(jacs), eval.Q(storage, y, μ, D, D_N), eval.R(storage),
                              eval.S(storage), eval.R_prod(storage), eval.S_prod(storage),
                              vars.q(qr), vars.r(qr), uk_eq, Jk, Kk, use_cholesky);
            // LQR solve
            vec work_2x(2 * nx);
            lqr.solve_masked(vars.AB(jacs), Jk, q, work_2x);
            return q;
        }
    };

    py::class_<OCPEvaluator>(m, "OCPEvaluator")
        .def(py::init<const typename OCPEvaluator::Problem &>())
        .def("forward_backward", &OCPEvaluator::forward_backward, "u"_a, "y"_a = py::none(),
             "μ"_a = py::none(),
             ":return: * Cost\n"
             "         * Gradient\n\n")
        .def("Qk", &OCPEvaluator::Qk, "k"_a, "u"_a, "y"_a = py::none(), "μ"_a = py::none())
        .def("Rk", &OCPEvaluator::Rk, "k"_a, "u"_a, "mask"_a)
        .def("Sk", &OCPEvaluator::Sk, "k"_a, "u"_a, "mask"_a)
        .def("lqr_factor_solve", &OCPEvaluator::lqr_factor_solve, "u"_a, "γ"_a, "y"_a = py::none(),
             "μ"_a = py::none());
}

template void register_ocp<alpaqa::EigenConfigd>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigf>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_ocp<alpaqa::EigenConfigq>(py::module_ &);
#endif
