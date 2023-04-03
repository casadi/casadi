#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/accelerators/anderson.hpp>
#include <alpaqa/util/check-dim.hpp>

#include <dict/stats-to-dict.hpp>
#include <params/params.hpp>

template <alpaqa::Config Conf>
void register_anderson(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    using Anderson = alpaqa::AndersonAccel<config_t>;
    py::class_<Anderson> anderson(m, "AndersonAccel",
                                  "C++ documentation :cpp:class:`alpaqa::AndersonAccel`");
    using AndersonParams = typename Anderson::Params;
    register_dataclass<AndersonParams>(
        anderson, "Params", "C++ documentation :cpp:class:`alpaqa::AndersonAccelParams`");

    anderson //
        .def(py::init([](params_or_dict<AndersonParams> params) {
                 return Anderson{var_kwargs_to_struct(params)};
             }),
             "params"_a)
        .def(py::init([](params_or_dict<AndersonParams> params, length_t n) {
                 return Anderson{var_kwargs_to_struct(params), n};
             }),
             "params"_a, "n"_a)
        .def_property_readonly("params", &Anderson::get_params)
        .def("__str__", &Anderson::get_name);
}

template void register_anderson<alpaqa::EigenConfigd>(py::module_ &);
template void register_anderson<alpaqa::EigenConfigf>(py::module_ &);
template void register_anderson<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_anderson<alpaqa::EigenConfigq>(py::module_ &);
#endif
