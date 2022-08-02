#include <alpaqa/accelerators/lbfgs.hpp>
#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <alpaqa/problem/problem.hpp>
#include <alpaqa/problem/wrapped-problem-with-counters.hpp>
#if ALPAQA_HAVE_CASADI
#include <alpaqa/interop/casadi/CasADiLoader.hpp>
#endif

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/eval.h>
#include <pybind11/gil.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
using namespace py::literals;
constexpr auto ret_ref_internal = py::return_value_policy::reference_internal;

#include <limits>
#include <string>
#include <string_view>

template <class T>
constexpr auto copy() {
    return [](const T &self) { return T{self}; };
};
template <class T>
constexpr auto deepcopy() {
    return [](const T &self, py::dict) { return T{self}; };
};
template <class T, class M>
constexpr auto wrap_ref(M T::*m) {
    return py::cpp_function([=](T &self) { return Eigen::Ref<M>(self.*m); }, ret_ref_internal);
}
template <class T, class M>
constexpr auto wrap_set(M T::*m) {
    return [=](T &self, const M &v) { self.*m = v; };
}

#include "kwargs-to-struct.hpp"
#include "trampolines.hpp"
#include "type-erased-inner-solver.hpp"
#include "type-erased-panoc-direction.hpp"

#include <alpaqa/inner/src/panoc.tpp>
#include <alpaqa/outer/src/alm.tpp>

void register_common_classes(py::module_ &m) {
    // ----------------------------------------------------------------------------------------- //
    py::class_<alpaqa::EvalCounter> evalcounter(m, "EvalCounter",
                                                "C++ documentation: "
                                                ":cpp:class:`alpaqa::EvalCounter`\n\n");
    py::class_<alpaqa::EvalCounter::EvalTimer>(evalcounter, "EvalTimer",
                                               "C++ documentation: "
                                               ":cpp:class:`alpaqa::EvalCounter::EvalTimer`\n\n")
        .def(py::pickle(
            [](const alpaqa::EvalCounter::EvalTimer &p) {   // __getstate__
                return py::make_tuple(p.f,                  //
                                      p.grad_f,             //
                                      p.f_grad_f,           //
                                      p.f_g,                //
                                      p.f_grad_f_g,         //
                                      p.grad_f_grad_g_prod, //
                                      p.g,                  //
                                      p.grad_g_prod,        //
                                      p.grad_gi,            //
                                      p.grad_L,             //
                                      p.hess_L_prod,        //
                                      p.hess_L,             //
                                      p.ψ,                  //
                                      p.grad_ψ,             //
                                      p.grad_ψ_from_ŷ,      //
                                      p.ψ_grad_ψ            //
                );
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 16)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::EvalCounter::EvalTimer;
                return T{
                    py::cast<decltype(T::f)>(t[0]),
                    py::cast<decltype(T::grad_f)>(t[1]),
                    py::cast<decltype(T::f_grad_f)>(t[2]),
                    py::cast<decltype(T::f_g)>(t[3]),
                    py::cast<decltype(T::f_grad_f_g)>(t[4]),
                    py::cast<decltype(T::grad_f_grad_g_prod)>(t[5]),
                    py::cast<decltype(T::g)>(t[6]),
                    py::cast<decltype(T::grad_g_prod)>(t[7]),
                    py::cast<decltype(T::grad_gi)>(t[8]),
                    py::cast<decltype(T::grad_L)>(t[9]),
                    py::cast<decltype(T::hess_L_prod)>(t[10]),
                    py::cast<decltype(T::hess_L)>(t[11]),
                    py::cast<decltype(T::ψ)>(t[12]),
                    py::cast<decltype(T::grad_ψ)>(t[13]),
                    py::cast<decltype(T::grad_ψ_from_ŷ)>(t[14]),
                    py::cast<decltype(T::ψ_grad_ψ)>(t[15]),
                };
            }))
        .def_readwrite("f", &alpaqa::EvalCounter::EvalTimer::f)
        .def_readwrite("grad_f", &alpaqa::EvalCounter::EvalTimer::grad_f)
        .def_readwrite("f_grad_f", &alpaqa::EvalCounter::EvalTimer::f_grad_f)
        .def_readwrite("f_g", &alpaqa::EvalCounter::EvalTimer::f_g)
        .def_readwrite("f_grad_f_g", &alpaqa::EvalCounter::EvalTimer::f_grad_f_g)
        .def_readwrite("grad_f_grad_g_prod", &alpaqa::EvalCounter::EvalTimer::grad_f_grad_g_prod)
        .def_readwrite("g", &alpaqa::EvalCounter::EvalTimer::g)
        .def_readwrite("grad_g_prod", &alpaqa::EvalCounter::EvalTimer::grad_g_prod)
        .def_readwrite("grad_gi", &alpaqa::EvalCounter::EvalTimer::grad_gi)
        .def_readwrite("grad_L", &alpaqa::EvalCounter::EvalTimer::grad_L)
        .def_readwrite("hess_L_prod", &alpaqa::EvalCounter::EvalTimer::hess_L_prod)
        .def_readwrite("hess_L", &alpaqa::EvalCounter::EvalTimer::hess_L)
        .def_readwrite("ψ", &alpaqa::EvalCounter::EvalTimer::ψ)
        .def_readwrite("grad_ψ", &alpaqa::EvalCounter::EvalTimer::grad_ψ)
        .def_readwrite("grad_ψ_from_ŷ", &alpaqa::EvalCounter::EvalTimer::grad_ψ_from_ŷ)
        .def_readwrite("ψ_grad_ψ", &alpaqa::EvalCounter::EvalTimer::ψ_grad_ψ);

    evalcounter
        .def(py::pickle(
            [](const alpaqa::EvalCounter &p) {              // __getstate__
                return py::make_tuple(p.f,                  //
                                      p.grad_f,             //
                                      p.f_grad_f,           //
                                      p.f_g,                //
                                      p.f_grad_f_g,         //
                                      p.grad_f_grad_g_prod, //
                                      p.g,                  //
                                      p.grad_g_prod,        //
                                      p.grad_gi,            //
                                      p.grad_L,             //
                                      p.hess_L_prod,        //
                                      p.hess_L,             //
                                      p.ψ,                  //
                                      p.grad_ψ,             //
                                      p.grad_ψ_from_ŷ,      //
                                      p.ψ_grad_ψ,           //
                                      p.time                //
                );
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 17)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::EvalCounter;
                return T{
                    py::cast<decltype(T::f)>(t[0]),
                    py::cast<decltype(T::grad_f)>(t[1]),
                    py::cast<decltype(T::f_grad_f)>(t[2]),
                    py::cast<decltype(T::f_g)>(t[3]),
                    py::cast<decltype(T::f_grad_f_g)>(t[4]),
                    py::cast<decltype(T::grad_f_grad_g_prod)>(t[5]),
                    py::cast<decltype(T::g)>(t[6]),
                    py::cast<decltype(T::grad_g_prod)>(t[7]),
                    py::cast<decltype(T::grad_gi)>(t[8]),
                    py::cast<decltype(T::grad_L)>(t[9]),
                    py::cast<decltype(T::hess_L_prod)>(t[10]),
                    py::cast<decltype(T::hess_L)>(t[11]),
                    py::cast<decltype(T::ψ)>(t[12]),
                    py::cast<decltype(T::grad_ψ)>(t[13]),
                    py::cast<decltype(T::grad_ψ_from_ŷ)>(t[14]),
                    py::cast<decltype(T::ψ_grad_ψ)>(t[15]),
                    py::cast<decltype(T::time)>(t[16]),
                };
            }))
        .def_readwrite("f", &alpaqa::EvalCounter::f)
        .def_readwrite("grad_f", &alpaqa::EvalCounter::grad_f)
        .def_readwrite("f_grad_f", &alpaqa::EvalCounter::f_grad_f)
        .def_readwrite("f_g", &alpaqa::EvalCounter::f_g)
        .def_readwrite("f_grad_f_g", &alpaqa::EvalCounter::f_grad_f_g)
        .def_readwrite("grad_f_grad_g_prod", &alpaqa::EvalCounter::grad_f_grad_g_prod)
        .def_readwrite("g", &alpaqa::EvalCounter::g)
        .def_readwrite("grad_g_prod", &alpaqa::EvalCounter::grad_g_prod)
        .def_readwrite("grad_gi", &alpaqa::EvalCounter::grad_gi)
        .def_readwrite("grad_L", &alpaqa::EvalCounter::grad_L)
        .def_readwrite("hess_L_prod", &alpaqa::EvalCounter::hess_L_prod)
        .def_readwrite("hess_L", &alpaqa::EvalCounter::hess_L)
        .def_readwrite("ψ", &alpaqa::EvalCounter::ψ)
        .def_readwrite("grad_ψ", &alpaqa::EvalCounter::grad_ψ)
        .def_readwrite("grad_ψ_from_ŷ", &alpaqa::EvalCounter::grad_ψ_from_ŷ)
        .def_readwrite("ψ_grad_ψ", &alpaqa::EvalCounter::ψ_grad_ψ)
        .def_readwrite("time", &alpaqa::EvalCounter::time)
        .def("__str__", [](const alpaqa::EvalCounter &c) {
            std::ostringstream os;
            os << c;
            return os.str();
        });

    py::enum_<alpaqa::LBFGSStepSize>(m, "LBFGSStepsize",
                                     "C++ documentation: :cpp:enum:`alpaqa::LBFGSStepSize`")
        .value("BasedOnGradientStepSize", alpaqa::LBFGSStepSize::BasedOnGradientStepSize)
        .value("BasedOnCurvature", alpaqa::LBFGSStepSize::BasedOnCurvature)
        .export_values();

    using SolverStatus = alpaqa::SolverStatus;
    py::enum_<SolverStatus>(m, "SolverStatus", py::arithmetic(),
                            "C++ documentation: :cpp:enum:`alpaqa::SolverStatus`")
        .value("Busy", SolverStatus::Busy, "In progress.")
        .value("Converged", SolverStatus::Converged, "Converged and reached given tolerance")
        .value("MaxTime", SolverStatus::MaxTime, "Maximum allowed execution time exceeded")
        .value("MaxIter", SolverStatus::MaxIter, "Maximum number of iterations exceeded")
        .value("NotFinite", SolverStatus::NotFinite, "Intermediate results were infinite or NaN")
        .value("NoProgress", SolverStatus::NoProgress, "No progress was made in the last iteration")
        .value("Interrupted", SolverStatus::Interrupted, "Solver was interrupted by the user")
        .export_values();

    py::enum_<alpaqa::PANOCStopCrit>(m, "PANOCStopCrit",
                                     "C++ documentation: :cpp:enum:`alpaqa::PANOCStopCrit`")
        .value("ApproxKKT", alpaqa::PANOCStopCrit::ApproxKKT)
        .value("ApproxKKT2", alpaqa::PANOCStopCrit::ApproxKKT2)
        .value("ProjGradNorm", alpaqa::PANOCStopCrit::ProjGradNorm)
        .value("ProjGradNorm2", alpaqa::PANOCStopCrit::ProjGradNorm2)
        .value("ProjGradUnitNorm", alpaqa::PANOCStopCrit::ProjGradUnitNorm)
        .value("ProjGradUnitNorm2", alpaqa::PANOCStopCrit::ProjGradUnitNorm2)
        .value("FPRNorm", alpaqa::PANOCStopCrit::FPRNorm)
        .value("FPRNorm2", alpaqa::PANOCStopCrit::FPRNorm2)
        .value("Ipopt", alpaqa::PANOCStopCrit::Ipopt)
        .value("LBFGSBpp", alpaqa::PANOCStopCrit::LBFGSBpp)
        .export_values();
}

template <alpaqa::Config Conf>
void register_classes_for(py::module_ &m) {
    USING_ALPAQA_CONFIG(Conf);

    // ----------------------------------------------------------------------------------------- //
    using Box = alpaqa::Box<config_t>;
    py::class_<Box>(m, "Box", "C++ documentation: :cpp:class:`alpaqa::Box`")
        .def(py::init([](length_t n) {
                 return Box{vec::Constant(n, alpaqa::inf<config_t>),
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

    // ----------------------------------------------------------------------------------------- //
    using ProblemBase = alpaqa::ProblemBase<config_t>;
    py::class_<ProblemBase, std::shared_ptr<ProblemBase>>(m, "ProblemBase")
        .def("__copy__", [](const ProblemBase &p) { return p.clone(); })
        .def("__deepcopy__", [](const ProblemBase &p, py::dict) { return p.clone(); })
        .def_readwrite("n", &ProblemBase::n, "Number of unknowns, dimension of :math:`x`")
        .def_readwrite("m", &ProblemBase::m,
                       "Number of general constraints, dimension of :math:`g(x)`")
        .def("eval_f", &ProblemBase::eval_f)
        .def("eval_grad_f", &ProblemBase::eval_grad_f)
        .def("eval_grad_f",
             [](const ProblemBase &p, crvec x) {
                 vec g(p.n);
                 p.eval_grad_f(x, g);
                 return g;
             })
        .def("eval_g", &ProblemBase::eval_g)
        .def("eval_g",
             [](const ProblemBase &p, crvec x) {
                 vec g(p.m);
                 p.eval_g(x, g);
                 return g;
             })
        .def("eval_grad_g_prod", &ProblemBase::eval_grad_g_prod)
        .def("eval_grad_g_prod",
             [](const ProblemBase &p, crvec x, crvec y) {
                 vec g(p.n);
                 p.eval_grad_g_prod(x, y, g);
                 return g;
             })
        .def("eval_ψ_ŷ", &ProblemBase::eval_ψ_ŷ)
        .def("eval_ψ_ŷ",
             [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
                 vec ŷ(p.m);
                 auto ψ = p.eval_ψ_ŷ(x, y, Σ, ŷ);
                 return std::make_tuple(ψ, ŷ);
             })
        .def("eval_grad_ψ_from_ŷ", &ProblemBase::eval_grad_ψ_from_ŷ)
        .def("eval_grad_ψ_from_ŷ",
             [](const ProblemBase &p, crvec x, crvec ŷ) {
                 vec grad_ψ(p.n), work(p.n);
                 p.eval_grad_ψ_from_ŷ(x, ŷ, grad_ψ, work);
                 return grad_ψ;
             })
        .def("eval_grad_ψ", &ProblemBase::eval_grad_ψ)
        .def("eval_grad_ψ",
             [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
                 vec grad_ψ(p.n), work_n(p.n), work_m(p.m);
                 p.eval_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
                 return grad_ψ;
             })
        .def("eval_ψ_grad_ψ", &ProblemBase::eval_ψ_grad_ψ)
        .def("eval_ψ_grad_ψ", [](const ProblemBase &p, crvec x, crvec y, crvec Σ) {
            vec grad_ψ(p.n), work_n(p.n), work_m(p.m);
            auto ψ = p.eval_ψ_grad_ψ(x, y, Σ, grad_ψ, work_n, work_m);
            return std::make_tuple(ψ, grad_ψ);
        });
    using Problem = alpaqa::Problem<config_t>;
    py::class_<Problem, ProblemBase, ProblemTrampoline<Problem>, std::shared_ptr<Problem>>(
        m, "Problem", "C++ documentation: :cpp:class:`Problem`")
        // .def(py::init())
        .def(py::init<length_t, length_t, length_t>(), "n"_a, "m"_a, "p"_a = 0,
             ":param n: Number of unknowns\n"
             ":param m: Number of constraints\n"
             ":param p: Number of parameters")
        .def_readwrite("C", &Problem::C, "Box constraints on :math:`x`")
        .def_readwrite("D", &Problem::D, "Box constraints on :math:`g(x)`")
        .def_property(
            "param", py::overload_cast<>(&Problem::get_param),
            [](Problem &p, crvec param) {
                if (param.size() != p.get_param().size())
                    throw std::invalid_argument("Invalid parameter dimension: got " +
                                                std::to_string(param.size()) + ", should be " +
                                                std::to_string(p.get_param().size()) + ".");
                p.set_param(param);
            },
            "Parameter vector :math:`p` of the problem");
    using CountedProblem = alpaqa::WrappedProblemWithCounters<config_t, std::shared_ptr>;
    py::class_<CountedProblem, ProblemBase, ProblemTrampoline<CountedProblem>,
               std::shared_ptr<CountedProblem>>(m, "CountedProblem")
        .def_readwrite("evaluations", &CountedProblem::evaluations);

    m.def(
        "with_counters",
        [](std::shared_ptr<ProblemBase> prob) { return CountedProblem{std::move(prob)}; }, "prob"_a,
        "Return a counted version of the given problem.");

    // ----------------------------------------------------------------------------------------- //

#if ALPAQA_HAVE_CASADI
    using CasADiProblem = alpaqa::CasADiProblem<config_t>;
#else
    class CasADiProblem : public Problem {};
#endif
    py::class_<CasADiProblem, Problem, ProblemTrampoline<CasADiProblem>,
               std::shared_ptr<CasADiProblem>>(
        m, "CasADiProblem",
        "C++ documentation: "
        ":cpp:class:`alpaqa::CasADiProblem`\n\n"
        "See :py:class:`alpaqa._alpaqa.Problem` for the full documentation.");

    // ----------------------------------------------------------------------------------------- //
    using TypeErasedPANOCDirection = alpaqa::TypeErasedPANOCDirection<Conf>;
    py::class_<TypeErasedPANOCDirection>(m, "PANOCDirection")
        .def("__str__", &TypeErasedPANOCDirection::template get_name<>);

    // ----------------------------------------------------------------------------------------- //
    using LBFGS = alpaqa::LBFGS<config_t>;
    py::class_<LBFGS> lbfgs(m, "LBFGS");
    using LBFGSParams = typename LBFGS::Params;
    py::class_<LBFGSParams> lbfgsparams(lbfgs, "Params");
    using CBFGS = alpaqa::CBFGSParams<config_t>;
    py::class_<CBFGS> cbfgs(lbfgsparams, "CBFGS");
    using LBFGSSign = typename LBFGS::Sign;
    py::enum_<LBFGSSign> lbfgssign(lbfgs, "Sign",
                                   "C++ documentation :cpp:enum:`alpaqa::LBFGS::Sign`");
    cbfgs //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<CBFGS>))
        .def("to_dict", &struct_to_dict<CBFGS>)
        .def_readwrite("α", &CBFGS::α)
        .def_readwrite("ϵ", &CBFGS::ϵ);
    lbfgsparams //
        .def(py::init())
        .def(py::init(&kwargs_to_struct<LBFGSParams>))
        .def("to_dict", &struct_to_dict<LBFGSParams>)
        .def_readwrite("memory", &LBFGSParams::memory)
        .def_readwrite("cbfgs", &LBFGSParams::cbfgs);
    lbfgssign //
        .value("Positive", LBFGSSign::Positive)
        .value("Negative", LBFGSSign::Negative)
        .export_values();
    lbfgs //
        .def(py::init<LBFGSParams>(), "params"_a)
        .def(py::init(
                 [](py::dict params) -> LBFGS { return {kwargs_to_struct<LBFGSParams>(params)}; }),
             "params"_a)
        .def(py::init<LBFGSParams, length_t>(), "params"_a, "n"_a)
        .def(py::init([](py::dict params, length_t n) -> LBFGS {
                 return {kwargs_to_struct<LBFGSParams>(params), n};
             }),
             "params"_a, "n"_a)
        .def_static("update_valid", LBFGS::update_valid, "params"_a, "yᵀs"_a, "sᵀs"_a, "pᵀp"_a)
        .def(
            "update",
            [](LBFGS &self, crvec xk, crvec xkp1, crvec pk, crvec pkp1, LBFGSSign sign,
               bool forced) {
                if (xk.size() != self.n())
                    throw std::invalid_argument("xk dimension mismatch");
                if (xkp1.size() != self.n())
                    throw std::invalid_argument("xkp1 dimension mismatch");
                if (pk.size() != self.n())
                    throw std::invalid_argument("pk dimension mismatch");
                if (pkp1.size() != self.n())
                    throw std::invalid_argument("pkp1 dimension mismatch");
                return self.update(xk, xkp1, pk, pkp1, sign, forced);
            },
            "xk"_a, "xkp1"_a, "pk"_a, "pkp1"_a, "sign"_a = LBFGS::Sign::Positive,
            "forced"_a = false)
        .def(
            "update_sy",
            [](LBFGS &self, crvec sk, crvec yk, real_t pkp1Tpkp1, bool forced) {
                if (sk.size() != self.n())
                    throw std::invalid_argument("sk dimension mismatch");
                if (yk.size() != self.n())
                    throw std::invalid_argument("yk dimension mismatch");
                return self.update_sy(sk, yk, pkp1Tpkp1, forced);
            },
            "sk"_a, "yk"_a, "pkp1Tpkp1"_a, "forced"_a = false)
        .def(
            "apply",
            [](LBFGS &self, rvec q, real_t γ) {
                if (q.size() != self.n())
                    throw std::invalid_argument("q dimension mismatch");
                return self.apply(q, γ);
            },
            "q"_a, "γ"_a)
        .def(
            "apply_masked",
            [](LBFGS &self, rvec q, real_t γ, const std::vector<index_t> &J) {
                return self.apply_masked(q, γ, J);
            },
            "q"_a, "γ"_a, "J"_a)
        .def("reset", &LBFGS::reset)
        .def("current_history", &LBFGS::current_history)
        .def("resize", &LBFGS::resize, "n"_a)
        .def("scale_y", &LBFGS::scale_y, "factor"_a)
        .def_property_readonly("n", &LBFGS::n)
        .def("s", [](LBFGS &self, index_t i) -> rvec { return self.s(i); })
        .def("y", [](LBFGS &self, index_t i) -> rvec { return self.y(i); })
        .def("ρ", [](LBFGS &self, index_t i) -> real_t & { return self.ρ(i); })
        .def("α", [](LBFGS &self, index_t i) -> real_t & { return self.α(i); })
        .def_property_readonly("params", &LBFGS::get_params)
        .def("__str__", &LBFGS::get_name);

    // ----------------------------------------------------------------------------------------- //
    using LipschitzEstimateParams = alpaqa::LipschitzEstimateParams<config_t>;
    py::class_<LipschitzEstimateParams>(
        m, "LipschitzEstimateParams",
        "C++ documentation: :cpp:class:`alpaqa::LipschitzEstimateParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<LipschitzEstimateParams>))
        .def("to_dict", &struct_to_dict<LipschitzEstimateParams>)
        .def_readwrite("L_0", &LipschitzEstimateParams::L_0)
        .def_readwrite("ε", &LipschitzEstimateParams::ε)
        .def_readwrite("δ", &LipschitzEstimateParams::δ)
        .def_readwrite("Lγ_factor", &LipschitzEstimateParams::Lγ_factor);

    using PANOCParams = alpaqa::PANOCParams<config_t>;
    py::class_<PANOCParams>(m, "PANOCParams", "C++ documentation: :cpp:class:`alpaqa::PANOCParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<PANOCParams>))
        .def("to_dict", &struct_to_dict<PANOCParams>)
        .def_readwrite("Lipschitz", &PANOCParams::Lipschitz)
        .def_readwrite("max_iter", &PANOCParams::max_iter)
        .def_readwrite("max_time", &PANOCParams::max_time)
        .def_readwrite("τ_min", &PANOCParams::τ_min)
        .def_readwrite("L_min", &PANOCParams::L_min)
        .def_readwrite("L_max", &PANOCParams::L_max)
        .def_readwrite("max_no_progress", &PANOCParams::max_no_progress)
        .def_readwrite("print_interval", &PANOCParams::print_interval)
        .def_readwrite("quadratic_upperbound_tolerance_factor",
                       &PANOCParams::quadratic_upperbound_tolerance_factor)
        .def_readwrite("update_lipschitz_in_linesearch",
                       &PANOCParams::update_lipschitz_in_linesearch)
        .def_readwrite("alternative_linesearch_cond", &PANOCParams::alternative_linesearch_cond)
        .def_readwrite("lbfgs_stepsize", &PANOCParams::lbfgs_stepsize);

    // ----------------------------------------------------------------------------------------- //
    using PANOCProgressInfo = alpaqa::PANOCProgressInfo<config_t>;
    py::class_<PANOCProgressInfo>(m, "PANOCProgressInfo",
                                  "Data passed to the PANOC progress callback.\n\n"
                                  "C++ documentation: :cpp:class:`alpaqa::PANOCProgressInfo`")
        .def_readonly("k", &PANOCProgressInfo::k, "Iteration")
        .def_readonly("x", &PANOCProgressInfo::x, "Decision variable :math:`x`")
        .def_readonly("p", &PANOCProgressInfo::p, "Projected gradient step :math:`p`")
        .def_readonly("norm_sq_p", &PANOCProgressInfo::norm_sq_p, ":math:`\\left\\|p\\right\\|^2`")
        .def_readonly("x̂", &PANOCProgressInfo::x̂,
                      "Decision variable after projected gradient step :math:`\\hat x`")
        .def_readonly("φγ", &PANOCProgressInfo::φγ,
                      "Forward-backward envelope :math:`\\varphi_\\gamma(x)`")
        .def_readonly("ψ", &PANOCProgressInfo::ψ, "Objective value :math:`\\psi(x)`")
        .def_readonly("grad_ψ", &PANOCProgressInfo::grad_ψ,
                      "Gradient of objective :math:`\\nabla\\psi(x)`")
        .def_readonly("ψ_hat", &PANOCProgressInfo::ψ_hat)
        .def_readonly("grad_ψ_hat", &PANOCProgressInfo::grad_ψ_hat)
        .def_readonly("L", &PANOCProgressInfo::L,
                      "Estimate of Lipschitz constant of objective :math:`L`")
        .def_readonly("γ", &PANOCProgressInfo::γ, "Step size :math:`\\gamma`")
        .def_readonly("τ", &PANOCProgressInfo::τ, "Line search parameter :math:`\\tau`")
        .def_readonly("ε", &PANOCProgressInfo::ε, "Tolerance reached :math:`\\varepsilon_k`")
        .def_readonly("Σ", &PANOCProgressInfo::Σ, "Penalty factor :math:`\\Sigma`")
        .def_readonly("y", &PANOCProgressInfo::y, "Lagrange multipliers :math:`y`")
        .def_property_readonly(
            "fpr", [](const PANOCProgressInfo &p) { return std::sqrt(p.norm_sq_p) / p.γ; },
            "Fixed-point residual :math:`\\left\\|p\\right\\| / \\gamma`");

    using PANOCSolver = alpaqa::PANOCSolver<TypeErasedPANOCDirection>;
    py::class_<PANOCSolver>(m, "PANOCSolver", "C++ documentation: :cpp:class:`alpaqa::PANOCSolver`")
        .def(py::init([](const PANOCParams &params, const LBFGS &lbfgs) {
            return PANOCSolver{params, alpaqa::erase_direction<LBFGS>(lbfgs)};
        }));

    using InnerSolver = alpaqa::TypeErasedInnerSolver<config_t>;
    py::class_<InnerSolver>(m, "InnerSolver")
        .def(py::init<PANOCSolver>())
        .def("__call__",
             [](InnerSolver &self, const alpaqa::ProblemBase<config_t> &p, crvec Σ, real_t ε,
                bool a, rvec x, vec y, rvec e) { return self(p, Σ, ε, a, x, y, e).to_dict(); })
        .def_property_readonly("name", &InnerSolver::template get_name<>);

    using ALMSolver = alpaqa::ALMSolver<InnerSolver>;
    using ALMParams = typename ALMSolver::Params;
    py::class_<ALMParams>(m, "ALMParams", "C++ documentation: :cpp:class:`alpaqa::ALMParams`")
        .def(py::init())
        .def(py::init(&kwargs_to_struct<ALMParams>))
        .def("to_dict", &struct_to_dict<ALMParams>)
        .def_readwrite("ε", &ALMParams::ε)
        .def_readwrite("δ", &ALMParams::δ)
        .def_readwrite("Δ", &ALMParams::Δ)
        .def_readwrite("Δ_lower", &ALMParams::Δ_lower)
        .def_readwrite("Σ_0", &ALMParams::Σ_0)
        .def_readwrite("σ_0", &ALMParams::σ_0)
        .def_readwrite("Σ_0_lower", &ALMParams::Σ_0_lower)
        .def_readwrite("ε_0", &ALMParams::ε_0)
        .def_readwrite("ε_0_increase", &ALMParams::ε_0_increase)
        .def_readwrite("ρ", &ALMParams::ρ)
        .def_readwrite("ρ_increase", &ALMParams::ρ_increase)
        .def_readwrite("θ", &ALMParams::θ)
        .def_readwrite("M", &ALMParams::M)
        .def_readwrite("Σ_max", &ALMParams::Σ_max)
        .def_readwrite("Σ_min", &ALMParams::Σ_min)
        .def_readwrite("max_iter", &ALMParams::max_iter)
        .def_readwrite("max_time", &ALMParams::max_time)
        .def_readwrite("max_num_initial_retries", &ALMParams::max_num_initial_retries)
        .def_readwrite("max_num_retries", &ALMParams::max_num_retries)
        .def_readwrite("max_total_num_retries", &ALMParams::max_total_num_retries)
        .def_readwrite("print_interval", &ALMParams::print_interval)
        .def_readwrite("single_penalty_factor", &ALMParams::single_penalty_factor);

    py::class_<ALMSolver>(m, "ALMSolver",
                          "Main augmented Lagrangian solver.\n\n"
                          "C++ documentation: :cpp:class:`alpaqa::ALMSolver`")
        // Default constructor
        .def(py::init(
                 []() -> ALMSolver { throw alpaqa::not_implemented_error("ALMSolver.__init__"); }),
             "Build an ALM solver using Structured PANOC as inner solver.")
        // Params and solver
        .def(py::init([](params_or_dict<ALMParams> params, const PANOCSolver &inner) {
                 return std::make_unique<ALMSolver>(var_kwargs_to_struct(params),
                                                    InnerSolver{inner});
             }),
             "alm_params"_a, "panoc_solver"_a, "Build an ALM solver using PANOC as inner solver.")
        // .def(
        //     py::init(PolymorphicALMConstructor<alpaqa::PolymorphicPGASolver>()),
        //     "alm_params"_a, "pga_solver"_a,
        //     "Build an ALM solver using the projected gradient algorithm as "
        //     "inner solver.")
        // .def(py::init(PolymorphicALMConstructor<
        //               alpaqa::PolymorphicStructuredPANOCLBFGSSolver>()),
        //      "alm_params"_a, "structuredpanoc_solver"_a,
        //      "Build an ALM solver using Structured PANOC as inner solver.")
        // .def(py::init(PolymorphicALMConstructor<
        //               alpaqa::PolymorphicInnerSolverTrampoline>()),
        //      "alm_params"_a, "inner_solver"_a,
        //      "Build an ALM solver using a custom inner solver.")
        // // Only solver (default params)
        // .def(py::init(PolymorphicALMConstructorDefaultParams<
        //               alpaqa::PolymorphicPANOCSolver>()),
        //      "panoc_solver"_a,
        //      "Build an ALM solver using PANOC as inner solver.")
        // .def(py::init(PolymorphicALMConstructorDefaultParams<
        //               alpaqa::PolymorphicPGASolver>()),
        //      "pga_solver"_a,
        //      "Build an ALM solver using the projected gradient algorithm as "
        //      "inner solver.")
        // .def(py::init(PolymorphicALMConstructorDefaultParams<
        //               alpaqa::PolymorphicStructuredPANOCLBFGSSolver>()),
        //      "structuredpanoc_solver"_a,
        //      "Build an ALM solver using Structured PANOC as inner solver.")
        // .def(py::init(PolymorphicALMConstructorDefaultParams<
        //               alpaqa::PolymorphicInnerSolverTrampoline>()),
        //      "inner_solver"_a,
        //      "Build an ALM solver using a custom inner solver.")
        // // Only params (default solver)
        // .def(py::init([](const alpaqa::ALMParams &params) {
        //          return alpaqa::PolymorphicALMSolver{
        //              params,
        //              std::static_pointer_cast<
        //                  alpaqa::PolymorphicInnerSolverBase>(
        //                  std::make_shared<
        //                      alpaqa::PolymorphicStructuredPANOCLBFGSSolver>(
        //                      alpaqa::StructuredPANOCLBFGSParams{},
        //                      alpaqa::LBFGSParams{})),
        //          };
        //      }),
        //      "alm_params"_a,
        //      "Build an ALM solver using Structured PANOC as inner solver.")
        // Other functions and properties
        .def_property_readonly(
            "inner_solver",
            py::cpp_function([](ALMSolver &self) -> InnerSolver & { return self.inner_solver; },
                             ret_ref_internal))
        .def(
            "__call__",
            [](ALMSolver &solver, const alpaqa::ProblemBase<config_t> &p, std::optional<vec> x,
               std::optional<vec> y) -> std::tuple<vec, vec, py::dict> {
                if (!x)
                    x = vec::Zero(p.n);
                else if (x->size() != p.n)
                    throw std::invalid_argument(
                        "Length of x does not match problem size problem.n");
                if (!y)
                    y = vec::Zero(p.m);
                else if (y->size() != p.m)
                    throw std::invalid_argument(
                        "Length of y does not match problem size problem.m");
                if (p.get_C().lowerbound.size() != p.n)
                    throw std::invalid_argument(
                        "Length of problem.C.lowerbound does not match problem "
                        "size problem.n");
                if (p.get_C().upperbound.size() != p.n)
                    throw std::invalid_argument(
                        "Length of problem.C.upperbound does not match problem "
                        "size problem.n");
                if (p.get_D().lowerbound.size() != p.m)
                    throw std::invalid_argument(
                        "Length of problem.D.lowerbound does not match problem "
                        "size problem.m");
                if (p.get_D().upperbound.size() != p.m)
                    throw std::invalid_argument(
                        "Length of problem.D.upperbound does not match problem "
                        "size problem.m");

                auto stats = solver(p, *y, *x);
                return std::make_tuple(std::move(*x), std::move(*y),
                                       alpaqa::conv::stats_to_dict<InnerSolver>(stats));
            },
            "problem"_a, "x"_a = std::nullopt, "y"_a = std::nullopt,
            py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>(),
            "Solve.\n\n"
            ":param problem: Problem to solve.\n"
            ":param y: Initial guess for Lagrange multipliers :math:`y`\n"
            ":param x: Initial guess for decision variables :math:`x`\n\n"
            ":return: * Lagrange multipliers :math:`y` at the solution\n"
            "         * Solution :math:`x`\n"
            "         * Statistics\n\n")
        .def("__str__", &ALMSolver::get_name)
        .def_property_readonly("params", &ALMSolver::get_params);

    if constexpr (std::is_base_of_v<alpaqa::EigenConfigd, config_t>) {

#if !ALPAQA_HAVE_CASADI
        constexpr static auto load_CasADi_problem = [](const char *, unsigned, unsigned, unsigned,
                                                       bool,
                                                       bool) -> std::unique_ptr<CasADiProblem> {
            throw std::runtime_error("This version of alpaqa was compiled without CasADi support");
        };
#else
        using CountedCasADiProblem = alpaqa::WrappedProblemWithCounters<config_t, std::unique_ptr>;
        constexpr static auto load_CasADi_problem = [](const char *so_name, unsigned n, unsigned m,
                                                       unsigned p, bool second_order) {
            return std::make_shared<CasADiProblem>(so_name, n, m, p, second_order);
        };
#endif
        m.def("load_casadi_problem", load_CasADi_problem, "so_name"_a, "n"_a = 0, "m"_a = 0,
              "p"_a = 0, "second_order"_a = false, "Load a compiled CasADi problem.\n\n");
    }
}

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc()               = "C++ implementation of alpaqa";
    m.attr("__version__") = VERSION_INFO;
#if ALPAQA_HAVE_CASADI
    m.attr("with_casadi") = true;
#else
    m.attr("with_casadi") = false;
#endif

    py::register_exception<alpaqa::not_implemented_error>(m, "not_implemented_error",
                                                          PyExc_NotImplementedError);

    register_common_classes(m);

    auto m_single = m.def_submodule("float32", "Single precision");
    register_classes_for<alpaqa::EigenConfigf>(m_single);
    auto m_double = m.def_submodule("float64", "Double precision");
    register_classes_for<alpaqa::EigenConfigd>(m_double);
    auto m_long_double = m.def_submodule("longdouble", "Long double precision");
    register_classes_for<alpaqa::EigenConfigl>(m_long_double);
#ifdef ALPAQA_WITH_QUAD_PRECISION
    auto m_quad = m.def_submodule("float128", "Quadruple precision (f128)");
    register_classes_for<alpaqa::EigenConfigq>(m_quad);
    // Note: this is usually disabled because NumPy doesn't support it.
#endif
}

// TODO: Why is this needed?
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::DefaultConfig>>;
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigf>>;
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigd>>;
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigl>>;
#ifdef ALPAQA_WITH_QUAD_PRECISION
template class alpaqa::PANOCSolver<alpaqa::TypeErasedPANOCDirection<alpaqa::EigenConfigq>>;
#endif