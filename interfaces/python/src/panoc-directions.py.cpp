#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>

#include "params/params.hpp"
#include "stats-to-dict.hpp"
#include "type-erased-panoc-direction.hpp"

template <alpaqa::Config Conf>
void register_panoc_directions(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    py::class_<TypeErasedPANOCDirection> te_direction(m, "PANOCDirection");
    te_direction //
        .def_property_readonly("params", &TypeErasedPANOCDirection::template get_params<>)
        .def("__str__", &TypeErasedPANOCDirection::template get_name<>);

    // ----------------------------------------------------------------------------------------- //
    using LBFGSDir       = alpaqa::LBFGSDirection<config_t>;
    using LBFGSParams    = alpaqa::LBFGSParams<config_t>;
    using LBFGSDirParams = alpaqa::LBFGSDirectionParams<config_t>;

    py::class_<LBFGSDir> lbfgs(m, "LBFGSDirection",
                               "C++ documentation: :cpp:class:`alpaqa::LBFGSDirection`");
    py::class_<LBFGSDirParams> lbfgs_params(
        lbfgs, "DirectionParams",
        "C++ documentation: :cpp:class:`alpaqa::LBFGSDirection::DirectionParams`");
    lbfgs_params //
        .def(py::init(&kwargs_to_struct<LBFGSDirParams>))
        .def("to_dict", &struct_to_dict<LBFGSDirParams>)
        .def_readwrite("rescale_when_γ_changes", &LBFGSDirParams::rescale_when_γ_changes);
    lbfgs //
        .def(py::init([](params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<LBFGSDirParams> direction_params) {
                 return LBFGSDir{var_kwargs_to_struct(lbfgs_params),
                                 var_kwargs_to_struct(direction_params)};
             }),
             "lbfgs_params"_a = py::dict{}, "direction_params"_a = py::dict{})
        .def_property_readonly(
            "params",
            py::cpp_function(&LBFGSDir::get_params, py::return_value_policy::reference_internal))
        .def("__str__", &LBFGSDir::get_name);

    te_direction.def(
        py::init(&alpaqa::erase_direction_with_params_dict<LBFGSDir, const LBFGSDir &>));
    py::implicitly_convertible<LBFGSDir, TypeErasedPANOCDirection>();

    // ----------------------------------------------------------------------------------------- //
    using StructuredLBFGSDir  = alpaqa::StructuredLBFGSDirection<config_t>;
    using StrucLBFGSDirParams = alpaqa::StructuredLBFGSDirectionParams<config_t>;

    py::class_<StructuredLBFGSDir> struc_lbfgs(
        m, "StructuredLBFGSDirection",
        "C++ documentation: :cpp:class:`alpaqa::StructuredLBFGSDirection`");
    py::class_<StrucLBFGSDirParams> struc_lbfgs_params(
        struc_lbfgs, "DirectionParams",
        "C++ documentation: :cpp:class:`alpaqa::StructuredLBFGSDirection::DirectionParams`");
    struc_lbfgs_params //
        .def(py::init(&kwargs_to_struct<StrucLBFGSDirParams>))
        .def("to_dict", &struct_to_dict<StrucLBFGSDirParams>)
        .def_readwrite("hessian_vec", &StrucLBFGSDirParams::hessian_vec)
        .def_readwrite("hessian_vec_finite_differences",
                       &StrucLBFGSDirParams::hessian_vec_finite_differences)
        .def_readwrite("full_augmented_hessian", &StrucLBFGSDirParams::full_augmented_hessian);
    struc_lbfgs //
        .def(py::init([](params_or_dict<LBFGSParams> lbfgs_params,
                         params_or_dict<StrucLBFGSDirParams> direction_params) {
                 return StructuredLBFGSDir{var_kwargs_to_struct(lbfgs_params),
                                           var_kwargs_to_struct(direction_params)};
             }),
             "lbfgs_params"_a = py::dict{}, "direction_params"_a = py::dict{})
        .def_property_readonly("params",
                               py::cpp_function(&StructuredLBFGSDir::get_params,
                                                py::return_value_policy::reference_internal))
        .def("__str__", &StructuredLBFGSDir::get_name);

    te_direction.def(py::init(
        &alpaqa::erase_direction_with_params_dict<StructuredLBFGSDir, const StructuredLBFGSDir &>));
    py::implicitly_convertible<StructuredLBFGSDir, TypeErasedPANOCDirection>();

    // ----------------------------------------------------------------------------------------- //
    // Catch-all, must be last
    te_direction //
        .def(py::init([](py::object o) {
            struct {
                using Problem = alpaqa::TypeErasedProblem<Conf>;
                void initialize(const Problem &problem, crvec y, crvec Σ, real_t γ_0, crvec x_0,
                                crvec x̂_0, crvec p_0, crvec grad_ψx_0) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("initialize")(problem, y, Σ, γ_0, x_0, x̂_0, p_0, grad_ψx_0);
                }
                bool update(real_t γₖ, real_t γₙₑₓₜ, crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ,
                            crvec grad_ψxₖ, crvec grad_ψxₙₑₓₜ) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(
                        o.attr("update")(γₖ, γₙₑₓₜ, xₖ, xₙₑₓₜ, pₖ, pₙₑₓₜ, grad_ψxₖ, grad_ψxₙₑₓₜ));
                }
                bool has_initial_direction() const {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(o.attr("has_initial_direction")());
                }
                bool apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ, rvec qₖ) const {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<bool>(o.attr("apply")(γₖ, xₖ, x̂ₖ, pₖ, grad_ψxₖ, qₖ));
                }
                void changed_γ(real_t γₖ, real_t old_γₖ) {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("changed_γ")(γₖ, old_γₖ);
                }
                void reset() {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    o.attr("reset")();
                }
                std::string get_name() const {
                    py::gil_scoped_acquire gil;
                    return py::cast<std::string>(py::str(o));
                }
                py::object get_params() const {
                    py::gil_scoped_acquire gil;
                    return py::getattr(o, "params");
                }

                py::object o;
            } s{std::move(o)};
            return TypeErasedPANOCDirection{std::move(s)};
        }));
}

template void register_panoc_directions<alpaqa::EigenConfigd>(py::module_ &);
template void register_panoc_directions<alpaqa::EigenConfigf>(py::module_ &);
template void register_panoc_directions<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_panoc_directions<alpaqa::EigenConfigq>(py::module_ &);
#endif
