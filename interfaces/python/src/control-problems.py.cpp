#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <variant>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/util/check-dim.hpp>
#if ALPAQA_HAVE_CASADI_OCP
#include <alpaqa/casadi/CasADiControlProblem.hpp>
#endif

template <alpaqa::Config Conf>
void register_control_problems(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);
    using alpaqa::util::check_dim;

    using ControlProblem = alpaqa::TypeErasedControlProblem<config_t>;
    py::class_<ControlProblem> te_problem(
        m, "ControlProblem", "C++ documentation: :cpp:class:`alpaqa::TypeErasedControlProblem`");
    te_problem //
        .def(py::init<const ControlProblem &>())
        .def("__copy__", [](const ControlProblem &self) { return ControlProblem{self}; })
        .def(
            "__deepcopy__",
            [](const ControlProblem &self, py::dict) { return ControlProblem{self}; }, "memo"_a)
        // TODO
        ;

    // ProblemWithCounters
    static constexpr auto te_pwc = []<class P>(P &&p) {
        using PwC = alpaqa::ControlProblemWithCounters<std::remove_cvref_t<P>>;
        auto te_p = ControlProblem::template make<PwC>(std::forward<P>(p));
        auto eval = te_p.template as<PwC>().evaluations;
        return std::make_tuple(std::move(te_p), std::move(eval));
    };

    if constexpr (std::is_same_v<typename Conf::real_t, double>) {
#if ALPAQA_HAVE_CASADI_OCP
        using CasADiControlProblem       = alpaqa::CasADiControlProblem<config_t>;
        auto load_CasADi_control_problem = [](const char *so_name, unsigned N) {
            return std::make_unique<CasADiControlProblem>(so_name, N);
        };
#else
        struct CasADiControlProblem {};
        auto load_CasADi_control_problem = [](const char *so_name,
                                              unsigned N) -> std::unique_ptr<CasADiControlProblem> {
            throw std::runtime_error("This version of alpaqa was compiled without CasADi support");
        };
#endif

        py::class_<CasADiControlProblem>(
            m, "CasADiControlProblem",
            "C++ documentation: :cpp:class:`alpaqa::CasADiControlProblem`\n\n"
            "See :py:class:`alpaqa._alpaqa.float64.TEControlProblem` for the full documentation.")
            .def("__copy__",
                 [](const CasADiControlProblem &self) { return CasADiControlProblem{self}; })
            .def(
                "__deepcopy__",
                [](const CasADiControlProblem &self, py::dict) {
                    return CasADiControlProblem{self};
                },
                "memo"_a)
#if ALPAQA_HAVE_CASADI_OCP
            .def_readonly("N", &CasADiControlProblem::N)
            .def_readonly("nx", &CasADiControlProblem::nx)
            .def_readonly("nu", &CasADiControlProblem::nu)
            .def_readonly("nh", &CasADiControlProblem::nh)
            .def_readonly("nh_N", &CasADiControlProblem::nh_N)
            .def_readonly("nc", &CasADiControlProblem::nc)
            .def_readonly("nc_N", &CasADiControlProblem::nc_N)
            .def_readwrite("U", &CasADiControlProblem::U)
            .def_readwrite("D", &CasADiControlProblem::D)
            .def_readwrite("D_N", &CasADiControlProblem::D_N)
            .def_property(
                "x_init", [](CasADiControlProblem &p) -> rvec { return p.x_init; },
                [](CasADiControlProblem &p, crvec x_init) {
                    if (x_init.size() != p.x_init.size())
                        throw std::invalid_argument("Invalid x_init dimension: got " +
                                                    std::to_string(x_init.size()) + ", should be " +
                                                    std::to_string(p.x_init.size()) + ".");
                    p.x_init = x_init;
                },
                "Initial state vector :math:`x^0` of the problem")
            .def_property(
                "param", [](CasADiControlProblem &p) -> rvec { return p.param; },
                [](CasADiControlProblem &p, crvec param) {
                    if (param.size() != p.param.size())
                        throw std::invalid_argument("Invalid parameter dimension: got " +
                                                    std::to_string(param.size()) + ", should be " +
                                                    std::to_string(p.param.size()) + ".");
                    p.param = param;
                },
                "Parameter vector :math:`p` of the problem")
#endif
            ;
#if ALPAQA_HAVE_CASADI_OCP
        te_problem.def(py::init<const CasADiControlProblem &>());
        py::implicitly_convertible<CasADiControlProblem, ControlProblem>();
#endif
        m.def("load_casadi_control_problem", load_CasADi_control_problem, "so_name"_a, "N"_a,
              "Load a compiled CasADi optimal control problem.\n\n");

#if ALPAQA_HAVE_CASADI_OCP
        m.def(
            "control_problem_with_counters",
            [](const CasADiControlProblem &p) { return te_pwc(p); }, "problem"_a,
            "Wrap the problem to count all function evaluations.\n\n"
            ":param problem: The original problem to wrap. Copied.\n"
            ":return: * Wrapped problem.\n"
            "         * Counters for wrapped problem.\n\n");
#endif
    }
}

template void register_control_problems<alpaqa::EigenConfigd>(py::module_ &);
template void register_control_problems<alpaqa::EigenConfigf>(py::module_ &);
template void register_control_problems<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_control_problems<alpaqa::EigenConfigq>(py::module_ &);
#endif