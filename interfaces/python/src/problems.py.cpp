#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/problem.hpp>
#include <alpaqa/problem/wrapped-problem-with-counters.hpp>

#include "trampolines.hpp"

template <alpaqa::Config Conf>
void register_problems(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using Box = alpaqa::Box<config_t>;
    py::class_<Box>(m, "Box", "C++ documentation: :cpp:class:`alpaqa::Box`")
        .def(py::init([](length_t n) {
                 return Box{vec::Constant(n, alpaqa::inf<config_t>),
                            vec::Constant(n, -alpaqa::inf<config_t>)};
             }),
             "n"_a,
             "Create an :math:`n`-dimensional box at with bounds at "
             ":math:`\\pm\\infty` (no constraints).")
        .def(py::init([](vec ub, vec lb) {
                 if (ub.size() != lb.size())
                     throw std::invalid_argument("Upper and lower bound dimensions do not match");
                 return Box{std::move(ub), std::move(lb)};
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

    using CountedProblem =
        alpaqa::WrappedProblemWithCounters<config_t, std::shared_ptr<const ProblemBase>>;
    py::class_<CountedProblem, ProblemBase, ProblemTrampoline<CountedProblem>,
               std::shared_ptr<CountedProblem>>(
        m, "CountedProblem", "C++ documentation: :cpp:class:`alpaqa::WrappedProblemWithCounters`")
        .def_readwrite("evaluations", &CountedProblem::evaluations);

    m.def(
        "with_counters",
        [](std::shared_ptr<ProblemBase> prob) { return CountedProblem{std::move(prob)}; }, "prob"_a,
        "Return a counted version of the given problem.");
}

template void register_problems<alpaqa::EigenConfigd>(py::module_ &);
template void register_problems<alpaqa::EigenConfigf>(py::module_ &);
template void register_problems<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_problems<alpaqa::EigenConfigq>(py::module_ &);
#endif