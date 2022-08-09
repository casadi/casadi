#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/problem.hpp>
#include <alpaqa/problem/wrapped-problem-with-counters.hpp>

#include "trampolines.hpp"

template <class FuncProb, auto py_f, auto f, class Ret, class... Args>
void functional_setter_ret(FuncProb &p, std::optional<py::object> o) {
    if (o) {
        p.*py_f = *std::move(o);
        p.*f    = [&pf {p.*py_f}](Args... x) -> Ret { return py::cast<Ret>(pf(x...)); };
    } else {
        p.*py_f = py::none();
        p.*f    = [](Args...) -> Ret {
            throw std::runtime_error("FunctionalProblem function is None");
        };
    }
};

template <class FuncProb, auto py_f, auto f, class Out, class Ret, class... Args>
void functional_setter_out(FuncProb &p, std::optional<py::object> o) {
    if (o) {
        p.*py_f = *std::move(o);
        p.*f    = [&pf {p.*py_f}](Args... x, Out r) -> void { r = py::cast<Ret>(pf(x...)); };
    } else {
        p.*py_f = py::none();
        p.*f    = [](Args..., Out) -> void {
            throw std::runtime_error("FunctionalProblem function is None");
        };
    }
};

template <alpaqa::Config Conf>
void register_problems(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using Box = alpaqa::Box<config_t>;
    py::class_<Box>(m, "Box", "C++ documentation: :cpp:class:`alpaqa::Box`")
        .def(py::init([](length_t n) {
                 return Box {vec::Constant(n, alpaqa::inf<config_t>),
                             vec::Constant(n, -alpaqa::inf<config_t>)};
             }),
             "n"_a,
             "Create an :math:`n`-dimensional box at with bounds at "
             ":math:`\\pm\\infty` (no constraints).")
        .def(py::init([](vec ub, vec lb) {
                 if (ub.size() != lb.size())
                     throw std::invalid_argument("Upper and lower bound dimensions do not match");
                 return Box {std::move(ub), std::move(lb)};
             }),
             "ub"_a, "lb"_a, "Create a box with the given bounds.")
        .def_readwrite("upperbound", &Box::upperbound)
        .def_readwrite("lowerbound", &Box::lowerbound);

    using ProblemBase = alpaqa::ProblemBase<config_t>;
    py::class_<ProblemBase, std::shared_ptr<ProblemBase>>(m, "ProblemBase")
        .def("__copy__", [](const ProblemBase &p) { return p.clone(); })
        .def("__deepcopy__", [](const ProblemBase &p, py::dict) { return p.clone(); })
        .def_readwrite("n", &ProblemBase::n, "Number of unknowns, dimension of :math:`x`")
        .def_readwrite("m", &ProblemBase::m,
                       "Number of general constraints, dimension of :math:`g(x)`")
        .def("eval_f", &ProblemBase::eval_f, "x"_a)
        .def("eval_grad_f", &ProblemBase::eval_grad_f, "x"_a, "grad_fx"_a)
        .def(
            "eval_grad_f",
            [](const ProblemBase &p, crvec x) {
                vec g(p.n);
                p.eval_grad_f(x, g);
                return g;
            },
            "x"_a)
        .def("eval_g", &ProblemBase::eval_g, "x"_a, "gx"_a)
        .def(
            "eval_g",
            [](const ProblemBase &p, crvec x) {
                vec g(p.m);
                p.eval_g(x, g);
                return g;
            },
            "x"_a)
        .def("eval_grad_g_prod", &ProblemBase::eval_grad_g_prod, "x"_a, "y"_a, "grad_gxy"_a)
        .def(
            "eval_grad_g_prod",
            [](const ProblemBase &p, crvec x, crvec y) {
                vec g(p.n);
                p.eval_grad_g_prod(x, y, g);
                return g;
            },
            "x"_a, "y"_a)
        .def("eval_ψ_ŷ", &ProblemBase::eval_ψ_ŷ, "x"_a, "y"_a, "Σ"_a, "ŷ"_a)
        .def(
            "eval_ψ_ŷ",
            [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
                vec ŷ(p.m);
                auto ψ = p.eval_ψ_ŷ(x, y, Σ, ŷ);
                return std::make_tuple(ψ, ŷ);
            },
            "x"_a, "y"_a, "Σ"_a)
        .def("eval_grad_ψ_from_ŷ", &ProblemBase::eval_grad_ψ_from_ŷ, "x"_a, "ŷ"_a, "grad_ψ"_a,
             "work_n"_a)
        .def(
            "eval_grad_ψ_from_ŷ",
            [](const ProblemBase &p, crvec x, crvec ŷ) {
                vec grad_ψ(p.n), work(p.n);
                p.eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work);
                return grad_ψ;
            },
            "x"_a, "ŷ"_a)
        .def("eval_grad_ψ", &ProblemBase::eval_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a,
             "work_m"_a)
        .def(
            "eval_grad_ψ",
            [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.n), work_n(p.n), work_m(p.m);
                p.eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return grad_ψ;
            },
            "x"_a, "y"_a, "Σ"_a)
        .def("eval_ψ_grad_ψ", &ProblemBase::eval_ψ_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a,
             "work_n"_a, "work_m"_a)
        .def(
            "eval_ψ_grad_ψ",
            [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.n), work_n(p.n), work_m(p.m);
                auto ψ = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return std::make_tuple(ψ, grad_ψ);
            },
            "x"_a, "y"_a, "Σ"_a);

    using Problem = alpaqa::Problem<config_t>;
    py::class_<Problem, ProblemBase, ProblemTrampoline<Problem>, std::shared_ptr<Problem>>(
        m, "Problem", "C++ documentation: :cpp:class:`alpaqa::Problem`")
        // .def(py::init())
        .def(py::init<length_t, length_t, length_t>(), "n"_a, "m"_a, "p"_a = 0,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints\n"
             ":param p: Number of parameters")
        .def_readwrite("C", &Problem::C, "Box constraints on :math:`x`")
        .def_readwrite("D", &Problem::D, "Box constraints on :math:`g(x)`")
        .def_property(
            "param", py::overload_cast<>(&Problem::get_param),
            [](Problem &p, crvec param) {
                if (param.size() != p.get_param().size())
                    throw std::invalid_argument("Invalid parameter dimension: got " +
                                                std::to_string(param.size()) + ", should be " +
                                                std::to_string(p.get_param().size()) + ".");
                p.set_param(param);
            },
            "Parameter vector :math:`p` of the problem");

    struct FunctionalProblem : alpaqa::FunctionalProblem<config_t> {
        FunctionalProblem(length_t n, length_t m, length_t p)
            : alpaqa::FunctionalProblem<config_t> {n, m, p} {
            this->f = [](crvec) -> real_t {
                throw std::runtime_error("FunctionalProblem.f is uninitialized");
            };
            this->grad_f = [](crvec, rvec) -> void {
                throw std::runtime_error("FunctionalProblem.grad_f is uninitialized");
            };
            this->g = [](crvec, rvec) -> void {
                throw std::runtime_error("FunctionalProblem.g is uninitialized");
            };
            this->grad_g_prod = [](crvec, crvec, rvec) -> void {
                throw std::runtime_error("FunctionalProblem.grad_g_prod is uninitialized");
            };
            this->grad_gi = [](crvec, index_t, rvec) -> void {
                throw std::runtime_error("FunctionalProblem.grad_gi is uninitialized");
            };
            this->hess_L_prod = [](crvec, crvec, crvec, rvec) -> void {
                throw std::runtime_error("FunctionalProblem.hess_L_prod is uninitialized");
            };
            this->hess_L = [](crvec, crvec, rmat) -> void {
                throw std::runtime_error("FunctionalProblem.hess_L is uninitialized");
            };
        }
        py::object py_f, py_grad_f, py_g, py_grad_g_prod, py_grad_gi, py_hess_L_prod, py_hess_L;
    };

    py::class_<FunctionalProblem, Problem, ProblemTrampoline<FunctionalProblem>,
               std::shared_ptr<FunctionalProblem>>(
        m, "FunctionalProblem", "C++ documentation: :cpp:class:`alpaqa::FunctionalProblem`")
        .def(py::init<length_t, length_t, length_t>(), "n"_a, "m"_a, "p"_a = 0,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints\n"
             ":param p: Number of parameters")
        .def_property(
            "f", [](const FunctionalProblem &p) { return p.py_f; },
            functional_setter_ret<FunctionalProblem, &FunctionalProblem::py_f,
                                  &FunctionalProblem::f, real_t, crvec>)
        .def_property(
            "grad_f", [](const FunctionalProblem &p) { return p.py_grad_f; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_grad_f,
                                  &FunctionalProblem::grad_f, rvec, crvec, crvec>)
        .def_property(
            "g", [](const FunctionalProblem &p) { return p.py_g; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_g,
                                  &FunctionalProblem::g, rvec, crvec, crvec>)
        .def_property(
            "grad_g_prod", [](const FunctionalProblem &p) { return p.py_grad_g_prod; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_grad_g_prod,
                                  &FunctionalProblem::grad_g_prod, rvec, crvec, crvec, crvec>)
        .def_property(
            "grad_gi", [](const FunctionalProblem &p) { return p.py_grad_gi; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_grad_gi,
                                  &FunctionalProblem::grad_gi, rvec, crvec, crvec, index_t>)
        .def_property(
            "hess_L_prod", [](const FunctionalProblem &p) { return p.py_hess_L_prod; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_hess_L_prod,
                                  &FunctionalProblem::hess_L_prod, rvec, crvec, crvec, crvec,
                                  crvec>)
        .def_property(
            "hess_L", [](const FunctionalProblem &p) { return p.py_hess_L; },
            functional_setter_out<FunctionalProblem, &FunctionalProblem::py_hess_L,
                                  &FunctionalProblem::hess_L, rmat, crmat, crvec, crvec>);

    using CountedProblem =
        alpaqa::WrappedProblemWithCounters<config_t, std::shared_ptr<const ProblemBase>>;
    py::class_<CountedProblem, ProblemBase, ProblemTrampoline<CountedProblem>,
               std::shared_ptr<CountedProblem>>(
        m, "CountedProblem", "C++ documentation: :cpp:class:`alpaqa::WrappedProblemWithCounters`")
        .def_readwrite("evaluations", &CountedProblem::evaluations);

    m.def(
        "with_counters",
        [](std::shared_ptr<ProblemBase> prob) { return CountedProblem {std::move(prob)}; },
        "prob"_a, "Return a counted version of the given problem.");
}

template void register_problems<alpaqa::EigenConfigd>(py::module_ &);
template void register_problems<alpaqa::EigenConfigf>(py::module_ &);
template void register_problems<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_problems<alpaqa::EigenConfigq>(py::module_ &);
#endif