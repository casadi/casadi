#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <variant>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/interop/casadi/CasADiLoader.hpp>
#endif

template <class FuncProb, auto py_f, auto f, class Ret, class... Args>
void functional_setter_ret(FuncProb &p, std::optional<py::object> o) {
    if (o) {
        p.*py_f = *std::move(o);
        p.*f    = [&pf{p.*py_f}](Args... x) -> Ret { return py::cast<Ret>(pf(x...)); };
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
        p.*f    = [&pf{p.*py_f}](Args... x, Out r) -> void { r = py::cast<Ret>(pf(x...)); };
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
        .def("__copy__", [](const Box &self) { return Box{self}; })
        .def(
            "__deepcopy__", [](const Box &self, py::dict) { return Box{self}; }, "memo"_a)
        .def(py::pickle(
            [](const Box &b) { // __getstate__
                return py::make_tuple(b.upperbound, b.lowerbound);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");
                return Box{
                    py::cast<decltype(Box::upperbound)>(t[0]),
                    py::cast<decltype(Box::lowerbound)>(t[1]),
                };
            }))
        .def(py::init([](length_t n) {
                 return Box{vec::Constant(n, +alpaqa::inf<config_t>),
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

    using BoxConstrProblem = alpaqa::BoxConstrProblem<config_t>;
    py::class_<BoxConstrProblem>(m, "BoxConstrProblem",
                                 "C++ documentation: :cpp:class:`alpaqa::BoxConstrProblem`")
        .def(py::init<length_t, length_t>(), "n"_a, "m"_a,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints")
        .def("__copy__", [](const BoxConstrProblem &self) { return BoxConstrProblem{self}; })
        .def(
            "__deepcopy__",
            [](const BoxConstrProblem &self, py::dict) { return BoxConstrProblem{self}; }, "memo"_a)
        .def(py::pickle(
            [](const BoxConstrProblem &self) { // __getstate__
                self.check();
                return py::make_tuple(self.C, self.D);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 2)
                    throw std::runtime_error("Invalid state!");
                return BoxConstrProblem{
                    py::cast<Box>(t[0]),
                    py::cast<Box>(t[1]),
                };
            }))
        .def_readwrite("n", &BoxConstrProblem::n,
                       "Number of decision variables, dimension of :math:`x`")
        .def_readwrite("m", &BoxConstrProblem::m,
                       "Number of general constraints, dimension of :math:`g(x)`")
        .def_readwrite("C", &BoxConstrProblem::C, "Box constraints on :math:`x`")
        .def_readwrite("D", &BoxConstrProblem::D, "Box constraints on :math:`g(x)`")
        .def("eval_proj_diff_g", &BoxConstrProblem::eval_proj_diff_g, "z"_a, "p"_a)
        .def("eval_proj_multipliers", &BoxConstrProblem::eval_proj_multipliers, "y"_a, "M"_a,
             "penalty_alm_split"_a)
        .def("eval_prox_grad_step", &BoxConstrProblem::eval_prox_grad_step, "γ"_a, "x"_a,
             "grad_ψ"_a, "x̂"_a, "p"_a)
        .def(
            "eval_proj_diff_g",
            [](const BoxConstrProblem &prob, crvec z) {
                vec p(prob.get_m());
                prob.eval_proj_diff_g(z, p);
                return p;
            },
            "z"_a)
        .def(
            "eval_prox_grad_step",
            [](const BoxConstrProblem &prob, real_t γ, crvec x, crvec grad_ψ) {
                vec x̂(prob.get_n());
                vec p(prob.get_n());
                prob.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
                return std::make_tuple(std::move(x̂), std::move(p));
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def("get_box_C", &BoxConstrProblem::get_box_C)
        .def("get_box_D", &BoxConstrProblem::get_box_D);

    struct PyProblem {
        USING_ALPAQA_CONFIG(Conf);
        py::object o;

        PyProblem(py::object o) : o{std::move(o)} {}

        // clang-format off
        void eval_proj_diff_g(crvec z, rvec p) const { o.attr("eval_proj_diff_g")(z, p); }
        void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const { o.attr("eval_proj_multipliers")(y, M, penalty_alm_split); }
        void eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { o.attr("eval_prox_grad_step")(γ, x, grad_ψ, x̂, p); }
        real_t eval_f(crvec x) const { return py::cast<real_t>(o.attr("eval_f")(x)); }
        void eval_grad_f(crvec x, rvec grad_fx) const { o.attr("eval_grad_f")(x, grad_fx); }
        void eval_g(crvec x, rvec gx) const { o.attr("eval_g")(x, gx); }
        void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { o.attr("eval_grad_g_prod")(x, y, grad_gxy); }
        void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const { o.attr("eval_grad_gi")(x, i, grad_gi); }
        void eval_hess_L_prod(crvec x, crvec y, crvec v, rvec Hv) const { o.attr("eval_hess_L_prod")(x, y, v, Hv); }
        void eval_hess_L(crvec x, crvec y, rmat H) const { o.attr("eval_hess_L")(x, y, H); }
        real_t eval_f_grad_f(crvec x, rvec grad_fx) const { return py::cast<real_t>(o.attr("eval_f_grad_f")(x, grad_fx)); }
        real_t eval_f_g(crvec x, rvec g) const { return py::cast<real_t>(o.attr("eval_f_g")(x, g)); }
        real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const { return py::cast<real_t>(o.attr("eval_f_grad_f_g")(x, grad_fx, g)); }
        void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const { o.attr("eval_grad_f_grad_g_prod")(x, y, grad_f, grad_gxy); }
        void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const { o.attr("eval_grad_L")(x, y, grad_L, work_n); }
        real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const { return py::cast<real_t>(o.attr("eval_ψ")(x, y, Σ, ŷ)); }
        void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const { o.attr("eval_grad_ψ_from_ŷ")(x, ŷ, grad_ψ, work_n); }
        void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { o.attr("eval_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m); }
        real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { return py::cast<real_t>(o.attr("eval_ψ_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m)); }
        void check() const { if (auto ch = py::getattr(o, "check", py::none()); !ch.is_none()) ch(); }
        const Box &get_box_C() const { C = o.attr("get_box_C")(); return py::cast<const Box &>(C); }
        const Box &get_box_D() const { D = o.attr("get_box_D")(); return py::cast<const Box &>(D); }

        bool provides_eval_grad_gi() const { return py::hasattr(o, "eval_grad_gi"); }
        bool provides_eval_hess_L_prod() const { return py::hasattr(o, "eval_hess_L_prod"); }
        bool provides_eval_hess_L() const { return py::hasattr(o, "eval_hess_L"); }
        bool provides_eval_f_grad_f() const { return py::hasattr(o, "eval_f_grad_f"); }
        bool provides_eval_f_g() const { return py::hasattr(o, "eval_f_g"); }
        bool provides_eval_f_grad_f_g() const { return py::hasattr(o, "eval_f_grad_f_g"); }
        bool provides_eval_grad_f_grad_g_prod() const { return py::hasattr(o, "eval_grad_f_grad_g_prod"); }
        bool provides_eval_grad_L() const { return py::hasattr(o, "eval_grad_L"); }
        bool provides_eval_ψ() const { return py::hasattr(o, "eval_ψ"); }
        bool provides_eval_grad_ψ_from_ŷ() const { return py::hasattr(o, "eval_grad_ψ_from_ŷ"); }
        bool provides_eval_grad_ψ() const { return py::hasattr(o, "eval_grad_ψ"); }
        bool provides_eval_ψ_grad_ψ() const { return py::hasattr(o, "eval_ψ_grad_ψ"); }
        bool provides_get_box_C() const { return py::hasattr(o, "get_box_C"); }
        bool provides_get_box_D() const { return py::hasattr(o, "get_box_D"); }
        // clang-format on

        length_t get_n() const { return py::cast<length_t>(o.attr("n")); }
        length_t get_m() const { return py::cast<length_t>(o.attr("m")); }

        // To keep the references to the boxes alive
        mutable py::object C;
        mutable py::object D;
    };

    using TEProblem = alpaqa::TypeErasedProblem<config_t>;
    py::class_<TEProblem> te_problem(m, "TEProblem",
                                     "C++ documentation: :cpp:class:`alpaqa::TypeErasedProblem`");
    te_problem //
        .def(py::init<const TEProblem &>())
        .def("__copy__", [](const TEProblem &self) { return TEProblem{self}; })
        .def(
            "__deepcopy__", [](const TEProblem &self, py::dict) { return TEProblem{self}; },
            "memo"_a)
        // clang-format off
        .def("eval_proj_diff_g", &TEProblem::eval_proj_diff_g, "z"_a, "p"_a)
        .def("eval_proj_multipliers", &TEProblem::eval_proj_multipliers, "y"_a, "M"_a, "penalty_alm_split"_a)
        .def("eval_prox_grad_step", &TEProblem::eval_prox_grad_step, "γ"_a, "x"_a, "grad_ψ"_a, "x̂"_a, "p"_a)
        .def("eval_f", &TEProblem::eval_f, "x"_a)
        .def("eval_grad_f", &TEProblem::eval_grad_f, "x"_a, "grad_fx"_a)
        .def("eval_g", &TEProblem::eval_g, "x"_a, "gx"_a)
        .def("eval_grad_g_prod", &TEProblem::eval_grad_g_prod, "x"_a, "y"_a, "grad_gxy"_a)
        .def("eval_grad_gi", &TEProblem::eval_grad_gi, "x"_a, "i"_a, "grad_gi"_a)
        .def("eval_hess_L_prod", &TEProblem::eval_hess_L_prod, "x"_a, "y"_a, "v"_a, "Hv"_a)
        .def("eval_hess_L", &TEProblem::eval_hess_L, "x"_a, "y"_a, "H"_a)
        .def("eval_f_grad_f", &TEProblem::eval_f_grad_f, "x"_a, "grad_fx"_a)
        .def("eval_f_g", &TEProblem::eval_f_g, "x"_a, "g"_a)
        .def("eval_f_grad_f_g", &TEProblem::eval_f_grad_f_g, "x"_a, "grad_fx"_a, "g"_a)
        .def("eval_grad_f_grad_g_prod", &TEProblem::eval_grad_f_grad_g_prod, "x"_a, "y"_a, "grad_f"_a, "grad_gxy"_a)
        .def("eval_grad_L", &TEProblem::eval_grad_L, "x"_a, "y"_a, "grad_L"_a, "work_n"_a)
        .def("eval_ψ", &TEProblem::eval_ψ, "x"_a, "y"_a, "Σ"_a, "ŷ"_a)
        .def("eval_grad_ψ_from_ŷ", &TEProblem::eval_grad_ψ_from_ŷ, "x"_a, "ŷ"_a, "grad_ψ"_a, "work_n"_a)
        .def("eval_grad_ψ", &TEProblem::eval_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a)
        .def("eval_ψ_grad_ψ", &TEProblem::eval_ψ_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a)
        .def("get_box_C", &TEProblem::get_box_C)
        .def("get_box_D", &TEProblem::get_box_D)

        .def("provides_eval_grad_gi", &TEProblem::provides_eval_grad_gi)
        .def("provides_eval_hess_L_prod", &TEProblem::provides_eval_hess_L_prod)
        .def("provides_eval_hess_L", &TEProblem::provides_eval_hess_L)
        .def("provides_eval_f_grad_f", &TEProblem::provides_eval_f_grad_f)
        .def("provides_eval_f_g", &TEProblem::provides_eval_f_g)
        .def("provides_eval_f_grad_f_g", &TEProblem::provides_eval_f_grad_f_g)
        .def("provides_eval_grad_f_grad_g_prod", &TEProblem::provides_eval_grad_f_grad_g_prod)
        .def("provides_eval_grad_L", &TEProblem::provides_eval_grad_L)
        .def("provides_eval_ψ", &TEProblem::provides_eval_ψ)
        .def("provides_eval_grad_ψ_from_ŷ", &TEProblem::provides_eval_grad_ψ_from_ŷ)
        .def("provides_eval_grad_ψ", &TEProblem::provides_eval_grad_ψ)
        .def("provides_eval_ψ_grad_ψ", &TEProblem::provides_eval_ψ_grad_ψ)
        .def("provides_get_box_C", &TEProblem::provides_get_box_C)
        .def("provides_get_box_D", &TEProblem::provides_get_box_D)
        // clang-format on
        .def(
            "eval_proj_diff_g",
            [](const TEProblem &prob, crvec z) {
                vec p(prob.get_m());
                prob.eval_proj_diff_g(z, p);
                return p;
            },
            "z"_a)
        .def(
            "eval_prox_grad_step",
            [](const TEProblem &prob, real_t γ, crvec x, crvec grad_ψ) {
                vec x̂(prob.get_n());
                vec p(prob.get_n());
                prob.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
                return std::make_tuple(std::move(x̂), std::move(p));
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def(
            "eval_grad_f",
            [](const TEProblem &p, crvec x) {
                vec g(p.get_n());
                p.eval_grad_f(x, g);
                return g;
            },
            "x"_a)
        .def(
            "eval_g",
            [](const TEProblem &p, crvec x) {
                vec g(p.get_m());
                p.eval_g(x, g);
                return g;
            },
            "x"_a)
        .def(
            "eval_grad_g_prod",
            [](const TEProblem &p, crvec x, crvec y) {
                vec g(p.get_n());
                p.eval_grad_g_prod(x, y, g);
                return g;
            },
            "x"_a, "y"_a)
        .def(
            "eval_ψ",
            [](const TEProblem &p, crvec x, crvec y, crvec Σ) {
                vec ŷ(p.get_m());
                auto ψ = p.eval_ψ(x, y, Σ, ŷ);
                return std::make_tuple(std::move(ψ), std::move(ŷ));
            },
            "x"_a, "y"_a, "Σ"_a)
        .def(
            "eval_grad_ψ_from_ŷ",
            [](const TEProblem &p, crvec x, crvec ŷ) {
                vec grad_ψ(p.get_n()), work(p.get_n());
                p.eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work);
                return grad_ψ;
            },
            "x"_a, "ŷ"_a)
        .def(
            "eval_grad_ψ",
            [](const TEProblem &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                p.eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return grad_ψ;
            },
            "x"_a, "y"_a, "Σ"_a)
        .def(
            "eval_ψ_grad_ψ",
            [](const TEProblem &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                auto ψ = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return std::make_tuple(std::move(ψ), std::move(grad_ψ));
            },
            "x"_a, "y"_a, "Σ"_a);

    using TEControlProblem = alpaqa::TypeErasedControlProblem<config_t>;
    py::class_<TEControlProblem> te_control_problem(
        m, "TEControlProblem", "C++ documentation: :cpp:class:`alpaqa::TypeErasedControlProblem`");
    te_control_problem //
        .def(py::init<const TEControlProblem &>())
        .def("__copy__", [](const TEControlProblem &self) { return TEControlProblem{self}; })
        .def(
            "__deepcopy__",
            [](const TEControlProblem &self, py::dict) { return TEControlProblem{self}; },
            "memo"_a);

    if constexpr (std::is_same_v<typename Conf::real_t, double>) {
#if ALPAQA_HAVE_CASADI
        using CasADiProblem      = alpaqa::CasADiProblem<config_t>;
        auto load_CasADi_problem = [](const char *so_name, unsigned n, unsigned m, unsigned p,
                                      bool second_order) {
            return std::make_unique<CasADiProblem>(so_name, n, m, p, second_order);
        };
#else
        class CasADiProblem : BoxConstrProblem {};
        auto load_CasADi_problem = [](const char *, unsigned, unsigned, unsigned,
                                      bool) -> std::unique_ptr<CasADiProblem> {
            throw std::runtime_error("This version of alpaqa was compiled without CasADi support");
        };
#endif

        py::class_<CasADiProblem, BoxConstrProblem>(
            m, "CasADiProblem",
            "C++ documentation: :cpp:class:`alpaqa::CasADiProblem`\n\n"
            "See :py:class:`alpaqa._alpaqa.float64.TEProblem` for the full documentation.")
            .def("__copy__", [](const CasADiProblem &self) { return CasADiProblem{self}; })
            .def(
                "__deepcopy__",
                [](const CasADiProblem &self, py::dict) { return CasADiProblem{self}; }, "memo"_a)
#if ALPAQA_HAVE_CASADI
            // clang-format off
            .def("eval_f", &CasADiProblem::eval_f, "x"_a)
            .def("eval_grad_f", &CasADiProblem::eval_grad_f, "x"_a, "grad_fx"_a)
            .def("eval_g", &CasADiProblem::eval_g, "x"_a, "gx"_a)
            .def("eval_grad_g_prod", &CasADiProblem::eval_grad_g_prod, "x"_a, "y"_a, "grad_gxy"_a)
            .def("eval_grad_gi", &CasADiProblem::eval_grad_gi, "x"_a, "i"_a, "grad_gi"_a)
            .def("eval_hess_L_prod", &CasADiProblem::eval_hess_L_prod, "x"_a, "y"_a, "v"_a, "Hv"_a)
            .def("eval_hess_L", &CasADiProblem::eval_hess_L, "x"_a, "y"_a, "H"_a)
            .def("eval_grad_L", &CasADiProblem::eval_grad_L, "x"_a, "y"_a, "grad_L"_a, "work_n"_a)
            .def("eval_ψ", &CasADiProblem::eval_ψ, "x"_a, "y"_a, "Σ"_a, "ŷ"_a)
            .def("eval_grad_ψ_from_ŷ", &CasADiProblem::eval_grad_ψ_from_ŷ, "x"_a, "ŷ"_a, "grad_ψ"_a, "work_n"_a)
            .def("eval_grad_ψ", &CasADiProblem::eval_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a)
            .def("eval_ψ_grad_ψ", &CasADiProblem::eval_ψ_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a)
            // clang-format on
            .def(
                "eval_grad_f",
                [](const CasADiProblem &p, crvec x) {
                    vec g(p.get_n());
                    p.eval_grad_f(x, g);
                    return g;
                },
                "x"_a)
            .def(
                "eval_g",
                [](const CasADiProblem &p, crvec x) {
                    vec g(p.get_m());
                    p.eval_g(x, g);
                    return g;
                },
                "x"_a)
            .def(
                "eval_grad_g_prod",
                [](const CasADiProblem &p, crvec x, crvec y) {
                    vec g(p.get_n());
                    p.eval_grad_g_prod(x, y, g);
                    return g;
                },
                "x"_a, "y"_a)
            .def(
                "eval_ψ",
                [](const CasADiProblem &p, crvec x, crvec y, crvec Σ) {
                    vec ŷ(p.get_m());
                    auto ψ = p.eval_ψ(x, y, Σ, ŷ);
                    return std::make_tuple(std::move(ψ), std::move(ŷ));
                },
                "x"_a, "y"_a, "Σ"_a)
            .def(
                "eval_grad_ψ_from_ŷ",
                [](const CasADiProblem &p, crvec x, crvec ŷ) {
                    vec grad_ψ(p.get_n()), work(p.get_n());
                    p.eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work);
                    return grad_ψ;
                },
                "x"_a, "ŷ"_a)
            .def(
                "eval_grad_ψ",
                [](const CasADiProblem &p, crvec x, crvec y, crvec Σ) {
                    vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                    p.eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                    return grad_ψ;
                },
                "x"_a, "y"_a, "Σ"_a)
            .def(
                "eval_ψ_grad_ψ",
                [](const CasADiProblem &p, crvec x, crvec y, crvec Σ) {
                    vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                    auto ψ = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                    return std::make_tuple(std::move(ψ), std::move(grad_ψ));
                },
                "x"_a, "y"_a, "Σ"_a)
            .def_property(
                "param", [](CasADiProblem &p) -> rvec { return p.param; },
                [](CasADiProblem &p, crvec param) {
                    if (param.size() != p.param.size())
                        throw std::invalid_argument("Invalid parameter dimension: got " +
                                                    std::to_string(param.size()) + ", should be " +
                                                    std::to_string(p.param.size()) + ".");
                    p.param = param;
                },
                "Parameter vector :math:`p` of the problem");
#endif
        ;
#if ALPAQA_HAVE_CASADI
        te_problem.def(py::init<const CasADiProblem &>());
        py::implicitly_convertible<CasADiProblem, TEProblem>();
#endif

        m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a, "n"_a = 0, "m"_a = 0,
              "p"_a = 0, "second_order"_a = false, "Load a compiled CasADi problem.\n\n");

#if ALPAQA_HAVE_CASADI
        using CasADiControlProblem       = alpaqa::CasADiControlProblem<config_t>;
        auto load_CasADi_control_problem = [](const char *so_name, unsigned N, unsigned nx,
                                              unsigned nu, unsigned p) {
            return std::make_unique<CasADiControlProblem>(so_name, N, nx, nu, p);
        };
#else
        class CasADiControlProblem {};
        auto load_CasADi_control_problem = [](const char *so_name, unsigned N, unsigned nx,
                                              unsigned nu,
                                              unsigned p) -> std::unique_ptr<CasADiControlProblem> {
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
#if ALPAQA_HAVE_CASADI
            .def_readonly("N", &CasADiControlProblem::N)
            .def_readonly("nx", &CasADiControlProblem::nx)
            .def_readonly("nu", &CasADiControlProblem::nu)
            .def_readwrite("U", &CasADiControlProblem::U)
            .def_readwrite("cost_structure", &CasADiControlProblem::cost_structure)
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
                "Parameter vector :math:`p` of the problem");
#endif
        ;
#if ALPAQA_HAVE_CASADI
        te_control_problem.def(py::init<const CasADiControlProblem &>());
        py::implicitly_convertible<CasADiControlProblem, TEControlProblem>();
#endif

        m.def("load_casadi_control_problem", load_CasADi_control_problem, "so_name"_a, "N"_a,
              "nx"_a = 0, "nu"_a = 0, "p"_a = 0,
              "Load a compiled CasADi optimal control problem.\n\n");

        static constexpr auto te_pwc = []<class P>(P &&p) {
            using PwC = alpaqa::ProblemWithCounters<std::remove_cvref_t<P>>;
            auto te_p = TEProblem::template make<PwC>(std::forward<P>(p));
            auto eval = te_p.template as<PwC>().evaluations;
            return std::make_tuple(std::move(te_p), std::move(eval));
        };
        m.def(
            "problem_with_counters", [](const CasADiProblem &p) { return te_pwc(p); }, "problem"_a,
            "Wrap the problem to count all function evaluations.\n\n"
            ":param problem: The original problem to wrap. Copied.\n"
            ":return: * Wrapped problem.\n"
            "         * Counters for wrapped problem.\n\n");
        m.def(
            "problem_with_counters", [](py::object p) { return te_pwc(PyProblem{std::move(p)}); },
            "problem"_a);
    }
    // Must be last
    te_problem.def(py::init([](py::object o) { return TEProblem::template make<PyProblem>(o); }));
}

template void register_problems<alpaqa::EigenConfigd>(py::module_ &);
template void register_problems<alpaqa::EigenConfigf>(py::module_ &);
template void register_problems<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_problems<alpaqa::EigenConfigq>(py::module_ &);
#endif