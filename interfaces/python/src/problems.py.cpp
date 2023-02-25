#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <variant>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/casadi/CasADiProblem.hpp>
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
    using alpaqa::util::check_dim;

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
                return Box::from_lower_upper(py::cast<decltype(Box::lowerbound)>(t[1]),
                                             py::cast<decltype(Box::upperbound)>(t[0]));
            }))
        .def(py::init<length_t>(), "n"_a,
             "Create an :math:`n`-dimensional box at with bounds at "
             ":math:`\\pm\\infty` (no constraints).")
        .def(py::init([](vec lower, vec upper) {
                 if (lower.size() != upper.size())
                     throw std::invalid_argument("Upper and lower bound dimensions do not match");
                 return Box::from_lower_upper(std::move(lower), std::move(upper));
             }),
             py::kw_only(), "lower"_a, "upper"_a, "Create a box with the given bounds.")
        .def_readwrite("lowerbound", &Box::lowerbound)
        .def_readwrite("upperbound", &Box::upperbound);

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
                real_t hx̂ = prob.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
                return std::make_tuple(std::move(x̂), std::move(p), hx̂);
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def("get_box_C", &BoxConstrProblem::get_box_C)
        .def("get_box_D", &BoxConstrProblem::get_box_D);

    struct PyProblem {
        USING_ALPAQA_CONFIG(Conf);
        py::object o;

        PyProblem(py::object o) : o{std::move(o)} {}

        // clang-format off
        void eval_proj_diff_g(crvec z, rvec p) const { py::gil_scoped_acquire gil; o.attr("eval_proj_diff_g")(z, p); }
        void eval_proj_multipliers(rvec y, real_t M, index_t penalty_alm_split) const { py::gil_scoped_acquire gil; o.attr("eval_proj_multipliers")(y, M, penalty_alm_split); }
        real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_prox_grad_step")(γ, x, grad_ψ, x̂, p)); }
        real_t eval_f(crvec x) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_f")(x)); }
        void eval_grad_f(crvec x, rvec grad_fx) const { py::gil_scoped_acquire gil; o.attr("eval_grad_f")(x, grad_fx); }
        void eval_g(crvec x, rvec gx) const { py::gil_scoped_acquire gil; o.attr("eval_g")(x, gx); }
        void eval_grad_g_prod(crvec x, crvec y, rvec grad_gxy) const { py::gil_scoped_acquire gil; o.attr("eval_grad_g_prod")(x, y, grad_gxy); }
        void eval_grad_gi(crvec x, index_t i, rvec grad_gi) const { py::gil_scoped_acquire gil; o.attr("eval_grad_gi")(x, i, grad_gi); }
        void eval_hess_L_prod(crvec x, crvec y, real_t scale, crvec v, rvec Hv) const { py::gil_scoped_acquire gil; o.attr("eval_hess_L_prod")(x, y, scale, v, Hv); }
        // void eval_hess_L(crvec x, crvec y, rmat H) const { py::gil_scoped_acquire gil; o.attr("eval_hess_L")(x, y, H); } // TODO
        void eval_hess_ψ_prod(crvec x, crvec y, crvec Σ, real_t scale, crvec v, rvec Hv) const { py::gil_scoped_acquire gil; o.attr("eval_hess_ψ_prod")(x, y, Σ, scale, v, Hv); }
        // void eval_hess_ψ(crvec x, crvec y, crvec Σ, rmat H) const { py::gil_scoped_acquire gil; o.attr("eval_hess_ψ")(x, y, Σ, H); } // TODO
        real_t eval_f_grad_f(crvec x, rvec grad_fx) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_f_grad_f")(x, grad_fx)); }
        real_t eval_f_g(crvec x, rvec g) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_f_g")(x, g)); }
        real_t eval_f_grad_f_g(crvec x, rvec grad_fx, rvec g) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_f_grad_f_g")(x, grad_fx, g)); }
        void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const { py::gil_scoped_acquire gil; o.attr("eval_grad_f_grad_g_prod")(x, y, grad_f, grad_gxy); }
        void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const { py::gil_scoped_acquire gil; o.attr("eval_grad_L")(x, y, grad_L, work_n); }
        real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_ψ")(x, y, Σ, ŷ)); }
        void eval_grad_ψ_from_ŷ(crvec x, crvec ŷ, rvec grad_ψ, rvec work_n) const { py::gil_scoped_acquire gil; o.attr("eval_grad_ψ_from_ŷ")(x, ŷ, grad_ψ, work_n); }
        void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { py::gil_scoped_acquire gil; o.attr("eval_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m); }
        real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_ψ_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m)); }
        void check() const { py::gil_scoped_acquire gil; if (auto ch = py::getattr(o, "check", py::none()); !ch.is_none()) ch(); }
        const Box &get_box_C() const { py::gil_scoped_acquire gil; alpaqa::ScopedMallocAllower ma; C = py::cast<Box>(o.attr("get_box_C")()); return C; }
        const Box &get_box_D() const { py::gil_scoped_acquire gil; alpaqa::ScopedMallocAllower ma; D = py::cast<Box>(o.attr("get_box_D")()); return D; }

        bool provides_eval_grad_gi() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_gi"); }
        bool provides_eval_hess_L_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_L_prod"); }
        // bool provides_eval_hess_L() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_L"); }
        bool provides_eval_hess_ψ_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_ψ_prod"); }
        // bool provides_eval_hess_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_ψ"); }
        bool provides_eval_f_grad_f() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_f_grad_f"); }
        bool provides_eval_f_g() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_f_g"); }
        bool provides_eval_f_grad_f_g() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_f_grad_f_g"); }
        bool provides_eval_grad_f_grad_g_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_f_grad_g_prod"); }
        bool provides_eval_grad_L() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_L"); }
        bool provides_eval_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_ψ"); }
        bool provides_eval_grad_ψ_from_ŷ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_ψ_from_ŷ"); }
        bool provides_eval_grad_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_ψ"); }
        bool provides_eval_ψ_grad_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_ψ_grad_ψ"); }
        bool provides_get_box_C() const { py::gil_scoped_acquire gil; return py::hasattr(o, "get_box_C"); }
        bool provides_get_box_D() const { py::gil_scoped_acquire gil; return py::hasattr(o, "get_box_D"); }

        length_t get_n() const { py::gil_scoped_acquire gil; return py::cast<length_t>(o.attr("n")); }
        length_t get_m() const { py::gil_scoped_acquire gil; return py::cast<length_t>(o.attr("m")); }
        // clang-format on

        // To keep the references to the boxes alive
        mutable Box C;
        mutable Box D;
    };

    using TEProblem = alpaqa::TypeErasedProblem<config_t>;
    py::class_<TEProblem> te_problem(m, "Problem",
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
        .def("eval_hess_L_prod", &TEProblem::eval_hess_L_prod, "x"_a, "y"_a, "scale"_a, "v"_a, "Hv"_a)
        // .def("eval_hess_L", &TEProblem::eval_hess_L, "x"_a, "y"_a, "H"_a) // TODO
        .def("eval_hess_ψ_prod", &TEProblem::eval_hess_ψ_prod, "x"_a, "y"_a, "Σ"_a, "scale"_a, "v"_a, "Hv"_a)
        // .def("eval_hess_ψ", &TEProblem::eval_hess_ψ, "x"_a, "y"_a, "Σ"_a, "H"_a) // TODO
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
        // .def("provides_eval_hess_L", &TEProblem::provides_eval_hess_L) // TODO
        .def("provides_eval_hess_ψ_prod", &TEProblem::provides_eval_hess_ψ_prod)
        // .def("provides_eval_hess_ψ", &TEProblem::provides_eval_hess_ψ) // TODO
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
                real_t hx̂ = prob.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
                return std::make_tuple(std::move(x̂), std::move(p), hx̂);
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

    // ProblemWithCounters
    static constexpr auto te_pwc = []<class P>(P &&p) {
        using PwC = alpaqa::ProblemWithCounters<std::remove_cvref_t<P>>;
        auto te_p = TEProblem::template make<PwC>(std::forward<P>(p));
        auto eval = te_p.template as<PwC>().evaluations;
        return std::make_tuple(std::move(te_p), std::move(eval));
    };

    if constexpr (std::is_same_v<typename Conf::real_t, double>) {
#if ALPAQA_HAVE_CASADI
        using CasADiProblem      = alpaqa::CasADiProblem<config_t>;
        auto load_CasADi_problem = [](const char *so_name) {
            return std::make_unique<CasADiProblem>(so_name);
        };
#else
        struct CasADiProblem : BoxConstrProblem {};
        auto load_CasADi_problem = [](const char *) -> std::unique_ptr<CasADiProblem> {
            throw std::runtime_error("This version of alpaqa was compiled without CasADi support");
        };
#endif

        py::class_<CasADiProblem, BoxConstrProblem>(
            m, "CasADiProblem",
            "C++ documentation: :cpp:class:`alpaqa::CasADiProblem`\n\n"
            "See :py:class:`alpaqa._alpaqa.float64.Problem` for the full documentation.")
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
            .def("eval_hess_L_prod", &CasADiProblem::eval_hess_L_prod, "x"_a, "y"_a, "scale"_a, "v"_a, "Hv"_a)
            // .def("eval_hess_L", &CasADiProblem::eval_hess_L, "x"_a, "y"_a, "H"_a) // TODO
            .def("eval_hess_ψ_prod", &CasADiProblem::eval_hess_ψ_prod, "x"_a, "y"_a, "Σ"_a, "scale"_a, "v"_a, "Hv"_a)
            // .def("eval_hess_ψ", &CasADiProblem::eval_hess_ψ, "x"_a, "y"_a, "Σ"_a, "H"_a) // TODO
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

        m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a,
              "Load a compiled CasADi problem.\n\n");

#if ALPAQA_HAVE_CASADI
        m.def(
            "problem_with_counters", [](const CasADiProblem &p) { return te_pwc(p); }, "problem"_a,
            "Wrap the problem to count all function evaluations.\n\n"
            ":param problem: The original problem to wrap. Copied.\n"
            ":return: * Wrapped problem.\n"
            "         * Counters for wrapped problem.\n\n");
#endif
    }
    m.def(
        "problem_with_counters", [](py::object p) { return te_pwc(PyProblem{std::move(p)}); },
        "problem"_a);

    m.def("provided_functions", [](const TEProblem &problem) {
        std::ostringstream os;
        alpaqa::print_provided_functions(os, problem);
        return os.str();
    });

    // Must be last
    te_problem.def(py::init([](py::object o) { return TEProblem::template make<PyProblem>(o); }));
}

template void register_problems<alpaqa::EigenConfigd>(py::module_ &);
template void register_problems<alpaqa::EigenConfigf>(py::module_ &);
template void register_problems<alpaqa::EigenConfigl>(py::module_ &);
#ifdef ALPAQA_WITH_QUAD_PRECISION
template void register_problems<alpaqa::EigenConfigq>(py::module_ &);
#endif