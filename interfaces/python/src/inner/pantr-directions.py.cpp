#include <alpaqa/util/quadmath/quadmath.hpp>

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/inner/directions/pantr/newton-tr.hpp>

#include <dict/stats-to-dict.hpp>
#include <inner/type-erased-tr-direction.hpp>
#include <params/params.hpp>

template <alpaqa::Config Conf>
void register_pantr_directions(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedTRDirection = alpaqa::TypeErasedTRDirection<Conf>;
    py::class_<TypeErasedTRDirection> te_direction(m, "PANTRDirection");
    te_direction //
        .def_property_readonly("params", &TypeErasedTRDirection::template get_params<>)
        .def("__str__", &TypeErasedTRDirection::template get_name<>);

    // ----------------------------------------------------------------------------------------- //
    using SteihaugCGParams = alpaqa::SteihaugCGParams<config_t>;
    register_dataclass<SteihaugCGParams>(
        m, "SteihaugCGParams", "C++ documentation: :cpp:class:`alpaqa::SteihaugCGParams`");

    using NewtonTRDirParams = alpaqa::NewtonTRDirectionParams<config_t>;
    register_dataclass<NewtonTRDirParams>(
        m, "NewtonTRDirectionParams",
        "C++ documentation: :cpp:class:`alpaqa::NewtonTRDirectionParams`");

    using NewtonTRDir = alpaqa::NewtonTRDirection<config_t>;
    py::class_<NewtonTRDir> newton_dir(m, "NewtonTRDirection",
                                       "C++ documentation: :cpp:class:`alpaqa::NewtonTRDirection`");
    newton_dir //
        .def(py::init([](params_or_dict<SteihaugCGParams> accelerator_params,
                         params_or_dict<NewtonTRDirParams> direction_params) {
                 return NewtonTRDir{var_kwargs_to_struct(accelerator_params),
                                    var_kwargs_to_struct(direction_params)};
             }),
             "accelerator_params"_a = py::dict{}, "direction_params"_a = py::dict{})
        .def_property_readonly(
            "params",
            py::cpp_function(&NewtonTRDir::get_params, py::return_value_policy::reference_internal))
        .def("__str__", &NewtonTRDir::get_name);

    te_direction.def(
        py::init(&alpaqa::erase_tr_direction_with_params_dict<NewtonTRDir, const NewtonTRDir &>));
    py::implicitly_convertible<NewtonTRDir, TypeErasedTRDirection>();

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
                real_t apply(real_t γₖ, crvec xₖ, crvec x̂ₖ, crvec pₖ, crvec grad_ψxₖ, real_t radius,
                             rvec qₖ) const {
                    alpaqa::ScopedMallocAllower ma;
                    py::gil_scoped_acquire gil;
                    return py::cast<real_t>(o.attr("apply")(γₖ, xₖ, x̂ₖ, pₖ, grad_ψxₖ, radius, qₖ));
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
            return TypeErasedTRDirection{std::move(s)};
        }));
}

template void register_pantr_directions<alpaqa::EigenConfigd>(py::module_ &);
template void register_pantr_directions<alpaqa::EigenConfigf>(py::module_ &);
template void register_pantr_directions<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_pantr_directions<alpaqa::EigenConfigq>(py::module_ &);
#endif
