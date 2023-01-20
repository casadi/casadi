#include <alpaqa/inner/directions/panoc-ocp/ocp-vars.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/check-dim.hpp>
#include <alpaqa/util/copyable_unique_ptr.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;

template <alpaqa::Config Conf>
void register_ocp(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    struct OCPEvaluator {
        using Problem = alpaqa::TypeErasedControlProblem<config_t>;
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
                alpaqa::util::check_dim<config_t>("μ", *μ, (nc == 0 && nc_N == 0) ? 0 : (N + 1));
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
    };

    py::class_<OCPEvaluator>(m, "OCPEvaluator")
        .def("forward_backward", &OCPEvaluator::forward_backward, "u"_a, "y"_a = py::none(),
             "μ"_a = py::none(),
             ":return: * Cost\n"
             "         * Gradient\n\n")
        .def("Qk", &OCPEvaluator::Qk, "k"_a, "u"_a, "y"_a = py::none(), "μ"_a = py::none())
        .def("Rk", &OCPEvaluator::Rk, "k"_a, "u"_a, "mask"_a)
        .def("Sk", &OCPEvaluator::Sk, "k"_a, "u"_a, "mask"_a);
}

template void register_ocp<alpaqa::EigenConfigd>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigf>(py::module_ &);
template void register_ocp<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_ocp<alpaqa::EigenConfigq>(py::module_ &);
#endif
