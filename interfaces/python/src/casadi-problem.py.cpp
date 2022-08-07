#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/problem.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/interop/casadi/CasADiLoader.hpp>
#endif

#include "trampolines.hpp"

template <alpaqa::Config Conf>
void register_casadi_problem(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using Problem     = alpaqa::Problem<config_t>;
    using ProblemBase = alpaqa::ProblemBase<config_t>;
#if ALPAQA_HAVE_CASADI
    using CasADiProblem      = alpaqa::CasADiProblem<config_t>;
    auto load_CasADi_problem = [](const char *so_name, unsigned n, unsigned m, unsigned p,
                                  bool second_order) {
        return std::make_shared<CasADiProblem>(so_name, n, m, p, second_order);
    };
#else
    class CasADiProblem : public Problem {};
    auto load_CasADi_problem = [](const char *, unsigned, unsigned, unsigned,
                                  bool) -> std::unique_ptr<CasADiProblem> {
        throw std::runtime_error("This version of alpaqa was compiled without CasADi support");
    };
#endif

    py::class_<CasADiProblem, Problem, ProblemTrampoline<CasADiProblem>,
               std::shared_ptr<CasADiProblem>>(
        m, "CasADiProblem",
        "C++ documentation: :cpp:class:`alpaqa::CasADiProblem`\n\n"
        "See :py:class:`alpaqa._alpaqa.Problem` for the full documentation.");

    m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a, "n"_a = 0, "m"_a = 0, "p"_a = 0,
          "second_order"_a = false, "Load a compiled CasADi problem.\n\n");
}

template void register_casadi_problem<alpaqa::EigenConfigd>(py::module_ &);
