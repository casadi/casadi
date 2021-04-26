#include <panoc-alm/inner/directions/lbfgs.hpp>
#include <panoc-alm/inner/panoc.hpp>
#include <panoc-alm/util/solverstatus.hpp>

#include <pybind11/attr.h>
#include <pybind11/cast.h>
#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace py = pybind11;

namespace pa {

class PolymorphicPANOCDirectionBase
    : public std::enable_shared_from_this<PolymorphicPANOCDirectionBase> {
  public:
    virtual ~PolymorphicPANOCDirectionBase()         = default;
    virtual void initialize(const vec &x₀, const vec &x̂₀, const vec &p₀,
                            const vec &grad₀)        = 0;
    virtual bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ,
                        const vec &pₖ₊₁, const vec &grad_new, const Box &C,
                        real_t γ_new)                = 0;
    virtual bool apply(const vec &xₖ, const vec &x̂ₖ, const vec &pₖ, real_t γ,
                       vec &qₖ)                      = 0;
    virtual void changed_γ(real_t γₖ, real_t old_γₖ) = 0;
    virtual void reset()                             = 0;
    virtual std::string get_name() const             = 0;
    vec apply_ret(const vec &xₖ, const vec &x̂ₖ, const vec &pₖ, real_t γ) {
        vec qₖ(pₖ.size());
        apply(xₖ, x̂ₖ, pₖ, γ, qₖ);
        return qₖ;
    }
};

class PolymorphicPANOCDirectionTrampoline
    : public PolymorphicPANOCDirectionBase {
  public:
    void initialize(const vec &x₀, const vec &x̂₀, const vec &p₀,
                    const vec &grad₀) override {
        PYBIND11_OVERRIDE_PURE(void, PolymorphicPANOCDirectionBase, initialize,
                               x₀, x̂₀, p₀, grad₀);
    }
    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &grad_new, const Box &C, real_t γ_new) override {
        PYBIND11_OVERRIDE_PURE(bool, PolymorphicPANOCDirectionBase, update, xₖ,
                               xₖ₊₁, pₖ, pₖ₊₁, grad_new, C, γ_new);
    }
    bool apply(const vec &xₖ, const vec &x̂ₖ, const vec &pₖ, real_t γ,
               vec &qₖ) override {
        PYBIND11_OVERRIDE_PURE(bool, PolymorphicPANOCDirectionBase, apply, xₖ,
                               x̂ₖ, pₖ, γ, qₖ);
    }
    void changed_γ(real_t γₖ, real_t old_γₖ) override {
        PYBIND11_OVERRIDE_PURE(void, PolymorphicPANOCDirectionBase, changed_γ,
                               γₖ, old_γₖ);
    }
    void reset() override {
        PYBIND11_OVERRIDE_PURE(void, PolymorphicPANOCDirectionBase, reset, );
    }
    std::string get_name() const override {
        PYBIND11_OVERRIDE_PURE(std::string, PolymorphicPANOCDirectionBase,
                               get_name, );
    }
};

template <>
struct PANOCDirection<PolymorphicPANOCDirectionBase> {
    using DirectionPtr = std::shared_ptr<PolymorphicPANOCDirectionBase>;
    DirectionPtr direction;

    PANOCDirection(const DirectionPtr &direction) : direction(direction) {}
    PANOCDirection(DirectionPtr &&direction)
        : direction(std::move(direction)) {}

    void initialize(const vec &x₀, const vec &x̂₀, const vec &p₀,
                    const vec &grad₀) {
        direction->initialize(x₀, x̂₀, p₀, grad₀);
    }
    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &grad_new, const Box &C, real_t γ_new) {
        return direction->update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, grad_new, C, γ_new);
    }
    bool apply(const vec &xₖ, const vec &x̂ₖ, const vec &pₖ, real_t γ, vec &qₖ) {
        return direction->apply(xₖ, x̂ₖ, pₖ, γ, qₖ);
    }
    void changed_γ(real_t γₖ, real_t old_γₖ) {
        direction->changed_γ(γₖ, old_γₖ);
    }
    void reset() { direction->reset(); }
    std::string get_name() const { return direction->get_name(); }
};

template <class DirectionProviderT>
class PolymorphicPANOCDirection : public PolymorphicPANOCDirectionBase {

  public:
    using DirectionProvider = DirectionProviderT;

    PolymorphicPANOCDirection(DirectionProvider &&direction)
        : direction_provider(std::forward<DirectionProvider>(direction)) {}
    PolymorphicPANOCDirection(const DirectionProvider &direction)
        : direction_provider(direction) {}

    void initialize(const vec &x₀, const vec &x̂₀, const vec &p₀,
                    const vec &grad₀) override {
        direction_provider.initialize(x₀, x̂₀, p₀, grad₀);
    }
    bool update(const vec &xₖ, const vec &xₖ₊₁, const vec &pₖ, const vec &pₖ₊₁,
                const vec &grad_new, const Box &C, real_t γ_new) override {
        return direction_provider.update(xₖ, xₖ₊₁, pₖ, pₖ₊₁, grad_new, C,
                                         γ_new);
    }
    bool apply(const vec &xₖ, const vec &x̂ₖ, const vec &pₖ, real_t γ,
               vec &qₖ) override {
        return direction_provider.apply(xₖ, x̂ₖ, pₖ, γ, qₖ);
    }
    void changed_γ(real_t γₖ, real_t old_γₖ) override {
        direction_provider.changed_γ(γₖ, old_γₖ);
    }
    void reset() override { direction_provider.reset(); }
    std::string get_name() const override {
        return direction_provider.get_name();
    }

  private:
    PANOCDirection<DirectionProvider> direction_provider;
};

using PolymorphicLBFGSDirection = PolymorphicPANOCDirection<LBFGS>;

template <class DirectionProviderT>
auto PolymorphicPANOCConstructor() {
    return [](const pa::PANOCParams &pp, const DirectionProviderT &dir) {
        using Base = PolymorphicPANOCDirectionBase;
        static_assert(std::is_base_of_v<Base, DirectionProviderT>);
        auto full_python_copy = std::make_shared<py::object>(py::cast(dir));
        auto base_copy        = full_python_copy->template cast<Base *>();
        return PANOCSolver<Base>{
            pp,
            std::shared_ptr<Base>(full_python_copy, base_copy),
        };
    };
}
template <class DirectionProviderT, class... DirectionArgumentsT>
auto PolymorphicPANOCConversion() {
    return [](const pa::PANOCParams &pp, const DirectionArgumentsT &...args) {
        using Base = PolymorphicPANOCDirectionBase;
        static_assert(std::is_base_of_v<Base, DirectionProviderT>);
        static_assert(std::is_constructible_v<DirectionProviderT,
                                              DirectionArgumentsT...>);
        DirectionProviderT dir{args...};
        return PolymorphicPANOCConstructor<DirectionProviderT>()(pp, dir);
    };
}

} // namespace pa

template <class T, class A>
auto attr_setter(A T::*attr) {
    return [attr](T &t, const py::handle &h) { t.*attr = h.cast<A>(); };
}
template <class T>
using attr_setter_fun_t = std::function<void(T &, const py::handle &)>;

template <class T>
using kwargs_to_struct_table_t = std::map<std::string, attr_setter_fun_t<T>>;

template <class T>
kwargs_to_struct_table_t<T> kwargs_to_struct_table;

template <>
const kwargs_to_struct_table_t<pa::PANOCParams>
    kwargs_to_struct_table<pa::PANOCParams>{
        {"Lipschitz", attr_setter(&pa::PANOCParams::Lipschitz)},
        {"max_iter", attr_setter(&pa::PANOCParams::max_iter)},
        {"max_time", attr_setter(&pa::PANOCParams::max_time)},
        {"τ_min", attr_setter(&pa::PANOCParams::τ_min)},
        {"max_no_progress", attr_setter(&pa::PANOCParams::max_no_progress)},
        {"print_interval", attr_setter(&pa::PANOCParams::print_interval)},
        {"quadratic_upperbound_tolerance_factor",
         attr_setter(&pa::PANOCParams::quadratic_upperbound_tolerance_factor)},
        {"update_lipschitz_in_linesearch",
         attr_setter(&pa::PANOCParams::update_lipschitz_in_linesearch)},
        {"alternative_linesearch_cond",
         attr_setter(&pa::PANOCParams::alternative_linesearch_cond)},
        {"lbfgs_stepsize", attr_setter(&pa::PANOCParams::lbfgs_stepsize)},
    };

template <>
const kwargs_to_struct_table_t<decltype(pa::PANOCParams::Lipschitz)>
    kwargs_to_struct_table<decltype(pa::PANOCParams::Lipschitz)>{
        {"L0", attr_setter(&decltype(pa::PANOCParams::Lipschitz)::L₀)},
        {"δ", attr_setter(&decltype(pa::PANOCParams::Lipschitz)::δ)},
        {"ε", attr_setter(&decltype(pa::PANOCParams::Lipschitz)::ε)},
        {"Lγ_factor",
         attr_setter(&decltype(pa::PANOCParams::Lipschitz)::Lγ_factor)},
    };

template <>
const kwargs_to_struct_table_t<pa::LBFGSParams>
    kwargs_to_struct_table<pa::LBFGSParams>{
        {"memory", attr_setter(&pa::LBFGSParams::memory)},
        {"cbfgs", attr_setter(&pa::LBFGSParams::cbfgs)},
        {"rescale_when_γ_changes",
         attr_setter(&pa::LBFGSParams::rescale_when_γ_changes)},
    };

template <>
const kwargs_to_struct_table_t<decltype(pa::LBFGSParams::cbfgs)>
    kwargs_to_struct_table<decltype(pa::LBFGSParams::cbfgs)>{
        {"α", attr_setter(&decltype(pa::LBFGSParams::cbfgs)::α)},
        {"ϵ", attr_setter(&decltype(pa::LBFGSParams::cbfgs)::ϵ)},
    };

template <class T>
void kwargs_to_struct_helper(T &t, const py::kwargs &kwargs) {
    const auto &m = kwargs_to_struct_table<T>;
    for (auto &&[key, val] : kwargs) {
        auto skey = key.template cast<std::string>();
        auto it   = m.find(skey);
        if (it == m.end())
            throw py::key_error("Unknown parameter " + skey);
        it->second(t, val);
    }
}

template <class T>
T kwargs_to_struct(const py::kwargs &kwargs) {
    T t{};
    kwargs_to_struct_helper(t, kwargs);
    return t;
}

class ostream_redirect {
  private:
    py::scoped_ostream_redirect stdout_stream{
        std::cout,                                // std::ostream&
        py::module_::import("sys").attr("stdout") // Python output
    };
    py::scoped_ostream_redirect stderr_stream{
        std::cerr,                                // std::ostream&
        py::module_::import("sys").attr("stderr") // Python output
    };

  public:
    ostream_redirect() = default;
};

PYBIND11_MODULE(PANOCPY_MODULE_NAME, m) {
    using py::operator""_a;

    m.doc() = "PANOC+ALM solvers"; // TODO

    py::class_<pa::Box>(m, "Box")
        .def_readwrite("upperbound", &pa::Box::upperbound)
        .def_readwrite("lowerbound", &pa::Box::lowerbound);

    py::class_<pa::Problem>(m, "Problem")
        // .def(py::init())
        .def(py::init<unsigned, unsigned>(), "n"_a, "m"_a)
        .def_readwrite("n", &pa::Problem::n)
        .def_readwrite("m", &pa::Problem::m)
        .def_readwrite("C", &pa::Problem::C)
        .def_readwrite("D", &pa::Problem::D)
        .def_property(
            "f",
            [](const pa::Problem &p) -> std::function<pa::real_t(pa::crvec)> {
                return [n{p.n}, f{p.f}](pa::crvec x) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    return f(x);
                };
            },
            [](pa::Problem &p,
               std::function<pa::real_t(pa::crvec)> fun) -> void { p.f = fun; })
        .def_property(
            "grad_f",
            [](const pa::Problem &p) -> std::function<pa::vec(pa::crvec)> {
                return [n{p.n}, grad_f{p.grad_f}](pa::crvec x) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    pa::vec gr(n);
                    grad_f(x, gr);
                    return gr;
                };
            },
            [](pa::Problem &p, std::function<pa::vec(pa::crvec)> fun) -> void {
                p.grad_f = [n{p.n}, fun{std::move(fun)}](pa::crvec x,
                                                         pa::rvec gr) {
                    auto &&res = fun(x);
                    if (res.size() != n)
                        throw std::out_of_range(
                            "Dimension of result of grad_f not consistent "
                            "with problem dimension n");
                    gr = std::move(res);
                };
            })
        .def_property(
            "g",
            [](const pa::Problem &p) -> std::function<pa::vec(pa::crvec)> {
                return [n{p.n}, m{p.m}, g{p.g}](pa::crvec x) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    pa::vec gg(m);
                    g(x, gg);
                    return gg;
                };
            },
            [](pa::Problem &p, std::function<pa::vec(pa::crvec)> fun) -> void {
                p.g = [m{p.m}, fun{std::move(fun)}](pa::crvec x, pa::rvec gg) {
                    auto &&res = fun(x);
                    if (res.size() != m)
                        throw std::out_of_range(
                            "Dimension of result of g not consistent "
                            "with problem dimension m");
                    gg = std::move(res);
                };
            })
        .def_property(
            "grad_g_prod",
            [](const pa::Problem &p)
                -> std::function<pa::vec(pa::crvec, pa::crvec)> {
                return [n{p.n}, m{p.m},
                        grad_g_prod{p.grad_g_prod}](pa::crvec x, pa::crvec y) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    if (y.size() != m)
                        throw std::out_of_range("Dimension of y not consistent "
                                                "with problem dimension m");
                    pa::vec gy(n);
                    grad_g_prod(x, y, gy);
                    return gy;
                };
            },
            [](pa::Problem &p,
               std::function<pa::vec(pa::crvec, pa::crvec)> fun) -> void {
                p.grad_g_prod = [n{p.n}, fun{std::move(fun)}](
                                    pa::crvec x, pa::crvec y, pa::rvec gy) {
                    auto &&res = fun(x, y);
                    if (res.size() != n)
                        throw std::out_of_range(
                            "Dimension of result of grad_g_prod not consistent "
                            "with problem dimension n");
                    gy = std::move(res);
                };
            })
        .def_property(
            "grad_gi",
            [](const pa::Problem &p)
                -> std::function<pa::vec(pa::crvec, unsigned)> {
                return [n{p.n}, m{p.m}, grad_gi{p.grad_gi}](pa::crvec x,
                                                            unsigned i) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    if (i < m)
                        throw std::out_of_range("Constraint index greater or "
                                                "equal to problem dimension m");
                    pa::vec gg(n);
                    grad_gi(x, i, gg);
                    return gg;
                };
            },
            [](pa::Problem &p,
               std::function<pa::vec(pa::crvec, unsigned)> fun) -> void {
                p.grad_gi = [n{p.n}, fun{std::move(fun)}](
                                pa::crvec x, unsigned i, pa::rvec gg) {
                    auto &&res = fun(x, i);
                    if (res.size() != n)
                        throw std::out_of_range(
                            "Dimension of result of grad_gi not consistent "
                            "with problem dimension n");
                    gg = std::move(res);
                };
            })
        .def_property(
            "hess_L",
            [](const pa::Problem &p)
                -> std::function<pa::mat(pa::crvec, pa::crvec)> {
                return [n{p.n}, m{p.m}, hess_L{p.hess_L}](pa::crvec x,
                                                          pa::crvec y) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    if (y.size() != m)
                        throw std::out_of_range("Dimension of y not consistent "
                                                "with problem dimension m");
                    pa::mat H(n, n);
                    hess_L(x, y, H);
                    return H;
                };
            },
            [](pa::Problem &p,
               std::function<pa::mat(pa::crvec, pa::crvec)> fun) -> void {
                p.hess_L = [n{p.n}, fun{std::move(fun)}](
                               pa::crvec x, pa::crvec y, pa::rmat H) {
                    auto &&res = fun(x, y);
                    if (res.rows() != n)
                        throw std::out_of_range(
                            "Number of rows of result of hess_L not consistent "
                            "with problem dimension n");
                    if (res.cols() != n)
                        throw std::out_of_range("Number of columns of result "
                                                "of hess_L not consistent "
                                                "with problem dimension n");
                    H = std::move(res);
                };
            })
        .def_property(
            "hess_L_prod",
            [](const pa::Problem &p)
                -> std::function<pa::vec(pa::crvec, pa::crvec, pa::crvec)> {
                return [n{p.n}, m{p.m}, hess_L_prod{p.hess_L_prod}](
                           pa::crvec x, pa::crvec y, pa::crvec v) {
                    if (x.size() != n)
                        throw std::out_of_range("Dimension of x not consistent "
                                                "with problem dimension n");
                    if (y.size() != m)
                        throw std::out_of_range("Dimension of y not consistent "
                                                "with problem dimension m");
                    if (v.size() != n)
                        throw std::out_of_range("Dimension of v not consistent "
                                                "with problem dimension n");
                    pa::vec Hv(n);
                    hess_L_prod(x, y, v, Hv);
                    return Hv;
                };
            },
            [](pa::Problem &p,
               std::function<pa::vec(pa::crvec, pa::crvec, pa::crvec)> fun)
                -> void {
                p.hess_L_prod = [n{p.n}, fun{std::move(fun)}](
                                    pa::crvec x, pa::crvec y, pa::crvec v,
                                    pa::rvec Hv) {
                    auto &&res = fun(x, y, v);
                    if (res.rows() != n)
                        throw std::out_of_range(
                            "Dimension of result of hess_L_prod not consistent "
                            "with problem dimension n");
                    Hv = std::move(res);
                };
            });

    py::class_<pa::PolymorphicPANOCDirectionBase,
               std::shared_ptr<pa::PolymorphicPANOCDirectionBase>,
               pa::PolymorphicPANOCDirectionTrampoline>(m, "PANOCDirection")
        .def(py::init<>())
        .def("initialize", &pa::PolymorphicPANOCDirectionBase::initialize)
        .def("update", &pa::PolymorphicPANOCDirectionBase::update)
        .def("apply", &pa::PolymorphicPANOCDirectionBase::apply_ret)
        .def("changed_γ", &pa::PolymorphicPANOCDirectionBase::changed_γ)
        .def("reset", &pa::PolymorphicPANOCDirectionBase::reset)
        .def("get_name", &pa::PolymorphicPANOCDirectionBase::get_name)
        .def("__str__", &pa::PolymorphicPANOCDirectionBase::get_name);

    py::class_<pa::PolymorphicLBFGSDirection,
               std::shared_ptr<pa::PolymorphicLBFGSDirection>,
               pa::PolymorphicPANOCDirectionBase>(m, "LBFGSDirection")
        .def(py::init<pa::LBFGSParams>())
        .def("initialize", &pa::PolymorphicLBFGSDirection::initialize)
        .def("update", &pa::PolymorphicLBFGSDirection::update)
        .def("apply", &pa::PolymorphicLBFGSDirection::apply_ret)
        .def("changed_γ", &pa::PolymorphicLBFGSDirection::changed_γ)
        .def("reset", &pa::PolymorphicLBFGSDirection::reset)
        .def("get_name", &pa::PolymorphicLBFGSDirection::get_name)
        .def("__str__", &pa::PolymorphicLBFGSDirection::get_name);

    using paLBFGSParamCBFGS = decltype(pa::LBFGSParams::cbfgs);
    py::class_<paLBFGSParamCBFGS>(m, "LBFGSParamsCBFGS")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<paLBFGSParamCBFGS>))
        .def_readwrite("α", &paLBFGSParamCBFGS::α)
        .def_readwrite("ϵ", &paLBFGSParamCBFGS::ϵ);

    py::class_<pa::LBFGSParams>(m, "LBFGSParams")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::LBFGSParams>))
        .def_readwrite("memory", &pa::LBFGSParams::memory)
        .def_readwrite("cbfgs", &pa::LBFGSParams::cbfgs)
        .def_readwrite("rescale_when_γ_changes",
                       &pa::LBFGSParams::rescale_when_γ_changes);

    py::enum_<pa::LBFGSStepSize>(m, "LBFGSStepsize")
        .value("BasedOnGradientStepSize",
               pa::LBFGSStepSize::BasedOnGradientStepSize)
        .value("BasedOnCurvature", pa::LBFGSStepSize::BasedOnCurvature)
        .export_values();

    using paPANOCParamsLipschitz = decltype(pa::PANOCParams::Lipschitz);
    py::class_<paPANOCParamsLipschitz>(m, "PANOCParamsLipschitz")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<paPANOCParamsLipschitz>))
        .def_readwrite("L0", &paPANOCParamsLipschitz::L₀)
        .def_readwrite("ε", &paPANOCParamsLipschitz::ε)
        .def_readwrite("δ", &paPANOCParamsLipschitz::δ)
        .def_readwrite("Lγ_factor", &paPANOCParamsLipschitz::Lγ_factor);

    py::class_<pa::PANOCParams>(m, "PANOCParams")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<pa::PANOCParams>))
        .def_readwrite("Lipschitz", &pa::PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &pa::PANOCParams::max_iter)
        .def_readwrite("max_time", &pa::PANOCParams::max_time)
        .def_readwrite("τ_min", &pa::PANOCParams::τ_min)
        .def_readwrite("max_no_progress", &pa::PANOCParams::max_no_progress)
        .def_readwrite("print_interval", &pa::PANOCParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &pa::PANOCParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("update_lipschitz_in_linesearch",
                       &pa::PANOCParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond",
                       &pa::PANOCParams::alternative_linesearch_cond)
        .def_readwrite("lbfgs_stepsize", &pa::PANOCParams::lbfgs_stepsize);

    py::enum_<pa::SolverStatus>(m, "SolverStatus", py::arithmetic(),
                                "Solver status")
        .value("Unknown", pa::SolverStatus::Unknown, "Initial value")
        .value("Converged", pa::SolverStatus::Converged,
               "Converged and reached given tolerance")
        .value("MaxTime", pa::SolverStatus::MaxTime,
               "Maximum allowed execution time exceeded")
        .value("MaxIter", pa::SolverStatus::MaxIter,
               "Maximum number of iterations exceeded")
        .value("NotFinite", pa::SolverStatus::NotFinite,
               "Intermediate results were infinite or not-a-number")
        .value("NoProgress", pa::SolverStatus::NoProgress,
               "No progress was made in the last iteration")
        .value("Interrupted", pa::SolverStatus::Interrupted,
               "Solver was interrupted by the user")
        .export_values();

    using pyPANOCSolver = pa::PANOCSolver<pa::PolymorphicPANOCDirectionBase>;
    py::class_<pyPANOCSolver::Stats>(m, "PANOCStats")
        .def_readwrite("status", &pyPANOCSolver::Stats::status)
        .def_readwrite("ε", &pyPANOCSolver::Stats::ε)
        .def_readwrite("elapsed_time", &pyPANOCSolver::Stats::elapsed_time)
        .def_readwrite("iterations", &pyPANOCSolver::Stats::iterations)
        .def_readwrite("linesearch_failures",
                       &pyPANOCSolver::Stats::linesearch_failures)
        .def_readwrite("lbfgs_failures", &pyPANOCSolver::Stats::lbfgs_failures)
        .def_readwrite("lbfgs_rejected", &pyPANOCSolver::Stats::lbfgs_rejected)
        .def_readwrite("τ_1_accepted", &pyPANOCSolver::Stats::τ_1_accepted)
        .def_readwrite("count_τ", &pyPANOCSolver::Stats::count_τ)
        .def_readwrite("sum_τ", &pyPANOCSolver::Stats::sum_τ);

    py::class_<pyPANOCSolver>(m, "PANOCSolver")
        .def(py::init(pa::PolymorphicPANOCConstructor< //
                      pa::PolymorphicLBFGSDirection>()))
        .def(py::init(pa::PolymorphicPANOCConversion< //
                      pa::PolymorphicLBFGSDirection,
                      pa::LBFGSParams>()))
        .def(py::init(pa::PolymorphicPANOCConstructor< //
                      pa::PolymorphicPANOCDirectionTrampoline>()))
        .def(
            "__call__",
            [](pyPANOCSolver &solver, const pa::Problem &p, pa::crvec Σ,
               pa::real_t ε, bool always_overwrite_results, pa::vec x,
               pa::vec y)
                -> std::tuple<pa::vec, pa::vec, pa::vec, pyPANOCSolver::Stats> {
                pa::vec z(p.m);
                auto stats = solver(p, Σ, ε, always_overwrite_results, x, y, z);
                return std::make_tuple(std::move(x), std::move(y), std::move(z),
                                       std::move(stats));
            },
            py::call_guard<ostream_redirect>())
        .def("__str__", &pyPANOCSolver::get_name);
}