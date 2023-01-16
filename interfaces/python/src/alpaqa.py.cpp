#include <alpaqa/config/config.hpp>
#include <alpaqa/util/not-implemented.hpp>

#include <pybind11/pybind11.h>
namespace py = pybind11;
using namespace py::literals;

void register_counters(py::module_ &m);
void register_enums(py::module_ &m);

template <alpaqa::Config Conf>
void register_problems(py::module_ &m);

template <alpaqa::Config Conf>
void register_control_problems(py::module_ &m);

template <alpaqa::Config Conf>
void register_lbfgs(py::module_ &m);

template <alpaqa::Config Conf>
void register_panoc_directions(py::module_ &m);

template <alpaqa::Config Conf>
void register_panoc(py::module_ &m);

template <alpaqa::Config Conf>
void register_panoc_ocp(py::module_ &m);

template <alpaqa::Config Conf>
void register_alm(py::module_ &m);

template <alpaqa::Config Conf>
void register_classes_for(py::module_ &m) {
    register_problems<Conf>(m);
    register_control_problems<Conf>(m);
    register_lbfgs<Conf>(m);
    register_panoc_directions<Conf>(m);
    register_panoc<Conf>(m);
    register_panoc_ocp<Conf>(m);
    register_alm<Conf>(m);
}

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc()               = "Python interface to alpaqa's C++ implementation.";
    m.attr("__version__") = VERSION_INFO;
    m.attr("build_time")  = __DATE__ " - " __TIME__;
#if ALPAQA_HAVE_CASADI
    m.attr("with_casadi") = true;
#else
    m.attr("with_casadi") = false;
#endif
#if ALPAQA_HAVE_CASADI_OCP
    m.attr("with_casadi_ocp") = true;
#else
    m.attr("with_casadi_ocp") = false;
#endif

    py::register_exception<alpaqa::not_implemented_error>(m, "not_implemented_error",
                                                          PyExc_NotImplementedError);

    register_counters(m);
    register_enums(m);

    auto m_single = m.def_submodule("float32", "Single precision");
    register_classes_for<alpaqa::EigenConfigf>(m_single);
    auto m_double = m.def_submodule("float64", "Double precision");
    register_classes_for<alpaqa::EigenConfigd>(m_double);
    auto m_long_double = m.def_submodule("longdouble", "Long double precision");
    register_classes_for<alpaqa::EigenConfigl>(m_long_double);
#ifdef ALPAQA_WITH_QUAD_PRECISION
    auto m_quad = m.def_submodule("float128", "Quadruple precision");
    register_classes_for<alpaqa::EigenConfigq>(m_quad);
    // Note: this is usually disabled because NumPy doesn't support it.
#endif
}
