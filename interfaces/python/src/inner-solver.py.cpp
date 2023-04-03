#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace py::literals;

#include "copy.hpp"
#include "inner-solve.hpp"
#include "member.hpp"
#include "type-erased-inner-solver.hpp"

template <alpaqa::Config Conf>
void register_inner_solver(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using TEProblem   = alpaqa::TypeErasedProblem<config_t>;
    using InnerSolver = alpaqa::TypeErasedInnerSolver<config_t, TEProblem>;
    py::class_<InnerSolver> inner_solver(m, "InnerSolver");
    default_copy_methods(inner_solver);
    inner_solver //
        .def_property_readonly("name", &InnerSolver::get_name)
        .def("stop", &InnerSolver::stop)
        .def("__str__",
             [](const InnerSolver &self) { return "InnerSolver<" + self.get_name() + ">"; });
    inner_solver_class<InnerSolver>.initialize(std::move(inner_solver));

    using TEOCProblem    = alpaqa::TypeErasedControlProblem<config_t>;
    using InnerOCPSolver = alpaqa::TypeErasedInnerSolver<config_t, TEOCProblem>;
    py::class_<InnerOCPSolver> inner_ocp_solver(m, "InnerOCPSolver");
    default_copy_methods(inner_ocp_solver);
    inner_ocp_solver.def_property_readonly("name", &InnerOCPSolver::get_name);
    inner_ocp_solver.def("stop", &InnerOCPSolver::stop);
    inner_ocp_solver.def("__str__", [](const InnerOCPSolver &self) {
        return "InnerOCPSolver<" + self.get_name() + ">";
    });
    inner_solver_class<InnerOCPSolver>.initialize(std::move(inner_ocp_solver));
}

template void register_inner_solver<alpaqa::EigenConfigd>(py::module_ &);
template void register_inner_solver<alpaqa::EigenConfigf>(py::module_ &);
template void register_inner_solver<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_inner_solver<alpaqa::EigenConfigq>(py::module_ &);
#endif
