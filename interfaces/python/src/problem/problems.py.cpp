#include <alpaqa/config/config.hpp>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>
#include <sstream>
#include <variant>

namespace py = pybind11;
using namespace py::literals;

#include <alpaqa/problem/box-constr-problem.hpp>
#include <alpaqa/problem/problem-with-counters.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/check-dim.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/casadi/CasADiProblem.hpp>
#endif

#include <util/copy.hpp>
#include <util/member.hpp>

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

template <class T, class... Args>
void problem_methods(py::class_<T, Args...> &cls) {
    USING_ALPAQA_CONFIG_TEMPLATE(T::config_t);
    // clang-format off
    cls.def("eval_proj_diff_g", &T::eval_proj_diff_g, "z"_a, "e"_a);
    cls.def("eval_proj_multipliers", &T::eval_proj_multipliers, "y"_a, "M"_a);
    cls.def("eval_prox_grad_step", &T::eval_prox_grad_step, "γ"_a, "x"_a, "grad_ψ"_a, "x̂"_a, "p"_a);
    cls.def("eval_inactive_indices_res_lna", &T::eval_inactive_indices_res_lna, "γ"_a, "x"_a, "grad_ψ"_a, "J"_a);
    cls.def("eval_f", &T::eval_f, "x"_a);
    cls.def("eval_grad_f", &T::eval_grad_f, "x"_a, "grad_fx"_a);
    cls.def("eval_g", &T::eval_g, "x"_a, "gx"_a);
    cls.def("eval_grad_g_prod", &T::eval_grad_g_prod, "x"_a, "y"_a, "grad_gxy"_a);
    cls.def("eval_grad_gi", &T::eval_grad_gi, "x"_a, "i"_a, "grad_gi"_a);
    cls.def("eval_hess_L_prod", &T::eval_hess_L_prod, "x"_a, "y"_a, "scale"_a, "v"_a, "Hv"_a);
    // cls.def("eval_hess_L", &T::eval_hess_L, "x"_a, "y"_a, "H"_a); // TODO
    cls.def("eval_hess_ψ_prod", &T::eval_hess_ψ_prod, "x"_a, "y"_a, "Σ"_a, "scale"_a, "v"_a, "Hv"_a);
    // cls.def("eval_hess_ψ", &T::eval_hess_ψ, "x"_a, "y"_a, "Σ"_a, "H"_a); // TODO
    cls.def("eval_f_grad_f", &T::eval_f_grad_f, "x"_a, "grad_fx"_a);
    if constexpr (requires { &T::eval_f_g; })
        cls.def("eval_f_g", &T::eval_f_g, "x"_a, "g"_a);
    if constexpr (requires { &T::eval_grad_f_grad_g_prod; })
        cls.def("eval_grad_f_grad_g_prod", &T::eval_grad_f_grad_g_prod, "x"_a, "y"_a, "grad_f"_a, "grad_gxy"_a);
    cls.def("eval_grad_L", &T::eval_grad_L, "x"_a, "y"_a, "grad_L"_a, "work_n"_a);
    cls.def("eval_ψ", &T::eval_ψ, "x"_a, "y"_a, "Σ"_a, "ŷ"_a);
    cls.def("eval_grad_ψ", &T::eval_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a);
    cls.def("eval_ψ_grad_ψ", &T::eval_ψ_grad_ψ, "x"_a, "y"_a, "Σ"_a, "grad_ψ"_a, "work_n"_a, "work_m"_a);
    cls.def("get_box_C", &T::get_box_C);
    cls.def("get_box_D", &T::get_box_D);

    if constexpr (requires { &T::provides_eval_inactive_indices_res_lna; })
        cls.def("provides_eval_inactive_indices_res_lna", &T::provides_eval_inactive_indices_res_lna);
    cls.def("provides_eval_grad_gi", &T::provides_eval_grad_gi);
    cls.def("provides_eval_hess_L_prod", &T::provides_eval_hess_L_prod);
    // cls.def("provides_eval_hess_L", &T::provides_eval_hess_L); // TODO
    cls.def("provides_eval_hess_ψ_prod", &T::provides_eval_hess_ψ_prod);
    // cls.def("provides_eval_hess_ψ", &T::provides_eval_hess_ψ); // TODO
    if constexpr (requires { &T::provides_eval_f_grad_f; })
        cls.def("provides_eval_f_grad_f", &T::provides_eval_f_grad_f);
    if constexpr (requires { &T::provides_eval_f_g; })
        cls.def("provides_eval_f_g", &T::provides_eval_f_g);
    if constexpr (requires { &T::provides_eval_grad_f_grad_g_prod; })
        cls.def("provides_eval_grad_f_grad_g_prod", &T::provides_eval_grad_f_grad_g_prod);
    if constexpr (requires { &T::provides_eval_grad_L; })
        cls.def("provides_eval_grad_L", &T::provides_eval_grad_L);
    if constexpr (requires { &T::provides_eval_ψ; })
        cls.def("provides_eval_ψ", &T::provides_eval_ψ);
    if constexpr (requires { &T::provides_eval_grad_ψ; })
        cls.def("provides_eval_grad_ψ", &T::provides_eval_grad_ψ);
    if constexpr (requires { &T::provides_eval_ψ_grad_ψ; })
        cls.def("provides_eval_ψ_grad_ψ", &T::provides_eval_ψ_grad_ψ);
    if constexpr (requires { &T::provides_get_box_C; })
        cls.def("provides_get_box_C", &T::provides_get_box_C);
    if constexpr (requires { &T::provides_get_box_D; })
        cls.def("provides_get_box_D", &T::provides_get_box_D);
    // clang-format on
    cls.def(
           "eval_proj_diff_g",
           [](const T &prob, crvec z) {
               vec e(prob.get_m());
               prob.eval_proj_diff_g(z, e);
               return e;
           },
           "z"_a)
        .def(
            "eval_prox_grad_step",
            [](const T &prob, real_t γ, crvec x, crvec grad_ψ) {
                vec x̂(prob.get_n());
                vec p(prob.get_n());
                real_t hx̂ = prob.eval_prox_grad_step(γ, x, grad_ψ, x̂, p);
                return std::make_tuple(std::move(x̂), std::move(p), hx̂);
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def(
            "eval_inactive_indices_res_lna",
            [](const T &prob, real_t γ, crvec x, crvec grad_ψ) {
                indexvec J_sto(prob.get_n());
                index_t nJ = prob.eval_inactive_indices_res_lna(γ, x, grad_ψ, J_sto);
                return indexvec{J_sto.topRows(nJ)};
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def(
            "eval_grad_f",
            [](const T &p, crvec x) {
                vec g(p.get_n());
                p.eval_grad_f(x, g);
                return g;
            },
            "x"_a)
        .def(
            "eval_g",
            [](const T &p, crvec x) {
                vec g(p.get_m());
                p.eval_g(x, g);
                return g;
            },
            "x"_a)
        .def(
            "eval_grad_g_prod",
            [](const T &p, crvec x, crvec y) {
                vec g(p.get_n());
                p.eval_grad_g_prod(x, y, g);
                return g;
            },
            "x"_a, "y"_a)
        .def(
            "eval_ψ",
            [](const T &p, crvec x, crvec y, crvec Σ) {
                vec ŷ(p.get_m());
                auto ψ = p.eval_ψ(x, y, Σ, ŷ);
                return std::make_tuple(std::move(ψ), std::move(ŷ));
            },
            "x"_a, "y"_a, "Σ"_a)
        .def(
            "eval_grad_ψ",
            [](const T &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                p.eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return grad_ψ;
            },
            "x"_a, "y"_a, "Σ"_a)
        .def(
            "eval_ψ_grad_ψ",
            [](const T &p, crvec x, crvec y, crvec Σ) {
                vec grad_ψ(p.get_n()), work_n(p.get_n()), work_m(p.get_m());
                auto ψ = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                return std::make_tuple(std::move(ψ), std::move(grad_ψ));
            },
            "x"_a, "y"_a, "Σ"_a);
}

template <alpaqa::Config Conf>
void register_problems(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);
    using alpaqa::util::check_dim;

    using Box = alpaqa::Box<config_t>;
    py::class_<Box> box(m, "Box", "C++ documentation: :cpp:class:`alpaqa::Box`");
    default_copy_methods(box);
    box //
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
        .def_property("lowerbound", vector_getter<&Box::lowerbound>(),
                      vector_setter<&Box::lowerbound>("lowerbound"))
        .def_property("upperbound", vector_getter<&Box::upperbound>(),
                      vector_setter<&Box::upperbound>("upperbound"));

    using BoxConstrProblem = alpaqa::BoxConstrProblem<config_t>;
    py::class_<BoxConstrProblem> box_constr_problem(
        m, "BoxConstrProblem", "C++ documentation: :cpp:class:`alpaqa::BoxConstrProblem`");
    default_copy_methods(box_constr_problem);
    box_constr_problem //
        .def(py::init<length_t, length_t>(), "n"_a, "m"_a,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints")
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
        .def_readwrite("l1_reg", &BoxConstrProblem::l1_reg,
                       py::return_value_policy::reference_internal,
                       ":math:`\\ell_1` regularization on :math:`x`")
        .def_readwrite("penalty_alm_split", &BoxConstrProblem::penalty_alm_split,
                       py::return_value_policy::reference_internal,
                       "Index between quadratic penalty and augmented Lagrangian constraints")
        .def("eval_proj_diff_g", &BoxConstrProblem::eval_proj_diff_g, "z"_a, "e"_a)
        .def("eval_proj_multipliers", &BoxConstrProblem::eval_proj_multipliers, "y"_a, "M"_a)
        .def("eval_prox_grad_step", &BoxConstrProblem::eval_prox_grad_step, "γ"_a, "x"_a,
             "grad_ψ"_a, "x̂"_a, "p"_a)
        .def("eval_inactive_indices_res_lna", &BoxConstrProblem::eval_inactive_indices_res_lna,
             "γ"_a, "x"_a, "grad_ψ"_a, "J"_a)
        .def(
            "eval_proj_diff_g",
            [](const BoxConstrProblem &prob, crvec z) {
                vec e(prob.get_m());
                prob.eval_proj_diff_g(z, e);
                return e;
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
        .def(
            "eval_inactive_indices_res_lna",
            [](const BoxConstrProblem &prob, real_t γ, crvec x, crvec grad_ψ) {
                indexvec J_sto(prob.get_n());
                index_t nJ = prob.eval_inactive_indices_res_lna(γ, x, grad_ψ, J_sto);
                return indexvec{J_sto.topRows(nJ)};
            },
            "γ"_a, "x"_a, "grad_ψ"_a)
        .def("get_box_C", &BoxConstrProblem::get_box_C)
        .def("get_box_D", &BoxConstrProblem::get_box_D);

    struct PyProblem {
        USING_ALPAQA_CONFIG(Conf);
        py::object o;

        PyProblem(py::object o) : o{std::move(o)} {}

        // clang-format off
        void eval_proj_diff_g(crvec z, rvec e) const { py::gil_scoped_acquire gil; o.attr("eval_proj_diff_g")(z, e); }
        void eval_proj_multipliers(rvec y, real_t M) const { py::gil_scoped_acquire gil; o.attr("eval_proj_multipliers")(y, M); }
        real_t eval_prox_grad_step(real_t γ, crvec x, crvec grad_ψ, rvec x̂, rvec p) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_prox_grad_step")(γ, x, grad_ψ, x̂, p)); }
        index_t eval_inactive_indices_res_lna(real_t γ, crvec x, crvec grad_ψ, rindexvec J) const { py::gil_scoped_acquire gil; return py::cast<index_t>(o.attr("eval_inactive_indices_res_lna")(γ, x, grad_ψ, J)); }
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
        void eval_grad_f_grad_g_prod(crvec x, crvec y, rvec grad_f, rvec grad_gxy) const { py::gil_scoped_acquire gil; o.attr("eval_grad_f_grad_g_prod")(x, y, grad_f, grad_gxy); }
        void eval_grad_L(crvec x, crvec y, rvec grad_L, rvec work_n) const { py::gil_scoped_acquire gil; o.attr("eval_grad_L")(x, y, grad_L, work_n); }
        real_t eval_ψ(crvec x, crvec y, crvec Σ, rvec ŷ) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_ψ")(x, y, Σ, ŷ)); }
        void eval_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { py::gil_scoped_acquire gil; o.attr("eval_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m); }
        real_t eval_ψ_grad_ψ(crvec x, crvec y, crvec Σ, rvec grad_ψ, rvec work_n, rvec work_m) const { py::gil_scoped_acquire gil; return py::cast<real_t>(o.attr("eval_ψ_grad_ψ")(x, y, Σ, grad_ψ, work_n, work_m)); }
        void check() const { py::gil_scoped_acquire gil; if (auto ch = py::getattr(o, "check", py::none()); !ch.is_none()) ch(); }
        const Box &get_box_C() const { py::gil_scoped_acquire gil; alpaqa::ScopedMallocAllower ma; C = py::cast<Box>(o.attr("get_box_C")()); return C; }
        const Box &get_box_D() const { py::gil_scoped_acquire gil; alpaqa::ScopedMallocAllower ma; D = py::cast<Box>(o.attr("get_box_D")()); return D; }

        bool provides_eval_inactive_indices_res_lna() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_inactive_indices_res_lna"); }
        bool provides_eval_grad_gi() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_gi"); }
        bool provides_eval_hess_L_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_L_prod"); }
        // bool provides_eval_hess_L() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_L"); }
        bool provides_eval_hess_ψ_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_ψ_prod"); }
        // bool provides_eval_hess_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_hess_ψ"); }
        bool provides_eval_f_grad_f() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_f_grad_f"); }
        bool provides_eval_f_g() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_f_g"); }
        bool provides_eval_grad_f_grad_g_prod() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_f_grad_g_prod"); }
        bool provides_eval_grad_L() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_grad_L"); }
        bool provides_eval_ψ() const { py::gil_scoped_acquire gil; return py::hasattr(o, "eval_ψ"); }
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
    default_copy_methods(te_problem);
    problem_methods(te_problem);

    // ProblemWithCounters
    struct ProblemWithCounters {
        TEProblem problem;
        std::shared_ptr<alpaqa::EvalCounter> evaluations;
    };
    py::class_<ProblemWithCounters>(m, "ProblemWithCounters")
        .def_readonly("problem", &ProblemWithCounters::problem)
        .def_readonly("evaluations", &ProblemWithCounters::evaluations);
    static constexpr auto te_pwc = []<class P>(P &&p) -> ProblemWithCounters {
        using PwC = alpaqa::ProblemWithCounters<P>;
        auto te_p = TEProblem::template make<PwC>(std::forward<P>(p));
        auto eval = te_p.template as<PwC>().evaluations;
        return {std::move(te_p), std::move(eval)};
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

        py::class_<CasADiProblem, BoxConstrProblem> casadi_problem(
            m, "CasADiProblem",
            "C++ documentation: :cpp:class:`alpaqa::CasADiProblem`\n\n"
            "See :py:class:`alpaqa._alpaqa.float64.Problem` for the full documentation.");
        default_copy_methods(casadi_problem);
#if ALPAQA_HAVE_CASADI
        problem_methods(casadi_problem);
        casadi_problem.def_property(
            "param", [](CasADiProblem &p) -> rvec { return p.param; },
            [](CasADiProblem &p, crvec param) {
                alpaqa::util::check_dim_msg<config_t>(param, p.param.size(),
                                                      "Invalid parameter size");
                p.param = param;
            },
            "Parameter vector :math:`p` of the problem");
#endif
#if ALPAQA_HAVE_CASADI
        te_problem.def(py::init<const CasADiProblem &>());
        py::implicitly_convertible<CasADiProblem, TEProblem>();
#endif

        m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a,
              "Load a compiled CasADi problem.\n\n");

#if ALPAQA_HAVE_CASADI
        m.def(
            "problem_with_counters", [](CasADiProblem &p) { return te_pwc(p); },
            py::keep_alive<0, 1>(), "problem"_a,
            "Wrap the problem to count all function evaluations.\n\n"
            ":param problem: The original problem to wrap. Copied.\n"
            ":return: * Wrapped problem.\n"
            "         * Counters for wrapped problem.\n\n");
#endif
    }
    m.def(
        "problem_with_counters", [](py::object p) { return te_pwc(PyProblem{std::move(p)}); },
        py::keep_alive<0, 1>(), "problem"_a);

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