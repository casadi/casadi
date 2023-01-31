#include <alpaqa/problem/ocproblem-counters.hpp>
#include <alpaqa/problem/problem-counters.hpp>
#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>
#include <memory>
#include <sstream>

namespace py = pybind11;

void register_counters(py::module_ &m) {
    // ----------------------------------------------------------------------------------------- //
    py::class_<alpaqa::EvalCounter, std::shared_ptr<alpaqa::EvalCounter>> evalcounter(
        m, "EvalCounter",
        "C++ documentation: "
        ":cpp:class:`alpaqa::EvalCounter`\n\n");
    py::class_<alpaqa::EvalCounter::EvalTimer>(evalcounter, "EvalTimer",
                                               "C++ documentation: "
                                               ":cpp:class:`alpaqa::EvalCounter::EvalTimer`\n\n")
        .def(py::pickle(
            [](const alpaqa::EvalCounter::EvalTimer &p) { // __getstate__
                return py::make_tuple(p.proj_diff_g, p.proj_multipliers, p.prox_grad_step, p.f,
                                      p.grad_f, p.f_grad_f, p.f_g, p.f_grad_f_g,
                                      p.grad_f_grad_g_prod, p.g, p.grad_g_prod, p.grad_gi, p.grad_L,
                                      p.hess_L_prod, p.hess_L, p.hess_ψ, p.ψ, p.grad_ψ,
                                      p.grad_ψ_from_ŷ, p.ψ_grad_ψ);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 20)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::EvalCounter::EvalTimer;
                return T{
                    .proj_diff_g        = py::cast<decltype(T::proj_diff_g)>(t[0]),
                    .proj_multipliers   = py::cast<decltype(T::proj_multipliers)>(t[1]),
                    .prox_grad_step     = py::cast<decltype(T::prox_grad_step)>(t[2]),
                    .f                  = py::cast<decltype(T::f)>(t[3]),
                    .grad_f             = py::cast<decltype(T::grad_f)>(t[4]),
                    .f_grad_f           = py::cast<decltype(T::f_grad_f)>(t[5]),
                    .f_g                = py::cast<decltype(T::f_g)>(t[6]),
                    .f_grad_f_g         = py::cast<decltype(T::f_grad_f_g)>(t[7]),
                    .grad_f_grad_g_prod = py::cast<decltype(T::grad_f_grad_g_prod)>(t[8]),
                    .g                  = py::cast<decltype(T::g)>(t[9]),
                    .grad_g_prod        = py::cast<decltype(T::grad_g_prod)>(t[10]),
                    .grad_gi            = py::cast<decltype(T::grad_gi)>(t[11]),
                    .grad_L             = py::cast<decltype(T::grad_L)>(t[12]),
                    .hess_L_prod        = py::cast<decltype(T::hess_L_prod)>(t[13]),
                    .hess_L             = py::cast<decltype(T::hess_L)>(t[14]),
                    .hess_ψ             = py::cast<decltype(T::hess_ψ)>(t[15]),
                    .ψ                  = py::cast<decltype(T::ψ)>(t[16]),
                    .grad_ψ             = py::cast<decltype(T::grad_ψ)>(t[17]),
                    .grad_ψ_from_ŷ      = py::cast<decltype(T::grad_ψ_from_ŷ)>(t[18]),
                    .ψ_grad_ψ           = py::cast<decltype(T::ψ_grad_ψ)>(t[19]),
                };
            }))
        .def_readwrite("proj_diff_g", &alpaqa::EvalCounter::EvalTimer::proj_diff_g)
        .def_readwrite("proj_multipliers", &alpaqa::EvalCounter::EvalTimer::proj_multipliers)
        .def_readwrite("prox_grad_step", &alpaqa::EvalCounter::EvalTimer::prox_grad_step)
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
        .def_readwrite("hess_ψ", &alpaqa::EvalCounter::EvalTimer::hess_ψ)
        .def_readwrite("ψ", &alpaqa::EvalCounter::EvalTimer::ψ)
        .def_readwrite("grad_ψ", &alpaqa::EvalCounter::EvalTimer::grad_ψ)
        .def_readwrite("grad_ψ_from_ŷ", &alpaqa::EvalCounter::EvalTimer::grad_ψ_from_ŷ)
        .def_readwrite("ψ_grad_ψ", &alpaqa::EvalCounter::EvalTimer::ψ_grad_ψ);

    evalcounter
        .def(py::pickle(
            [](const alpaqa::EvalCounter &p) { // __getstate__
                return py::make_tuple(p.proj_diff_g, p.proj_multipliers, p.prox_grad_step, p.f,
                                      p.grad_f, p.f_grad_f, p.f_g, p.f_grad_f_g,
                                      p.grad_f_grad_g_prod, p.g, p.grad_g_prod, p.grad_gi, p.grad_L,
                                      p.hess_L_prod, p.hess_L, p.hess_ψ, p.ψ, p.grad_ψ,
                                      p.grad_ψ_from_ŷ, p.ψ_grad_ψ, p.time);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 21)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::EvalCounter;
                return T{
                    .proj_diff_g        = py::cast<decltype(T::proj_diff_g)>(t[0]),
                    .proj_multipliers   = py::cast<decltype(T::proj_multipliers)>(t[1]),
                    .prox_grad_step     = py::cast<decltype(T::prox_grad_step)>(t[2]),
                    .f                  = py::cast<decltype(T::f)>(t[3]),
                    .grad_f             = py::cast<decltype(T::grad_f)>(t[4]),
                    .f_grad_f           = py::cast<decltype(T::f_grad_f)>(t[5]),
                    .f_g                = py::cast<decltype(T::f_g)>(t[6]),
                    .f_grad_f_g         = py::cast<decltype(T::f_grad_f_g)>(t[7]),
                    .grad_f_grad_g_prod = py::cast<decltype(T::grad_f_grad_g_prod)>(t[8]),
                    .g                  = py::cast<decltype(T::g)>(t[9]),
                    .grad_g_prod        = py::cast<decltype(T::grad_g_prod)>(t[10]),
                    .grad_gi            = py::cast<decltype(T::grad_gi)>(t[11]),
                    .grad_L             = py::cast<decltype(T::grad_L)>(t[12]),
                    .hess_L_prod        = py::cast<decltype(T::hess_L_prod)>(t[13]),
                    .hess_L             = py::cast<decltype(T::hess_L)>(t[14]),
                    .hess_ψ             = py::cast<decltype(T::hess_ψ)>(t[15]),
                    .ψ                  = py::cast<decltype(T::ψ)>(t[16]),
                    .grad_ψ             = py::cast<decltype(T::grad_ψ)>(t[17]),
                    .grad_ψ_from_ŷ      = py::cast<decltype(T::grad_ψ_from_ŷ)>(t[18]),
                    .ψ_grad_ψ           = py::cast<decltype(T::ψ_grad_ψ)>(t[19]),
                    .time               = py::cast<decltype(T::time)>(t[20]),
                };
            }))
        .def_readwrite("proj_diff_g", &alpaqa::EvalCounter::proj_diff_g)
        .def_readwrite("proj_multipliers", &alpaqa::EvalCounter::proj_multipliers)
        .def_readwrite("prox_grad_step", &alpaqa::EvalCounter::prox_grad_step)
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
        .def_readwrite("hess_ψ", &alpaqa::EvalCounter::hess_ψ)
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

    py::class_<alpaqa::OCPEvalCounter, std::shared_ptr<alpaqa::OCPEvalCounter>> ocpevalcounter(
        m, "OCPEvalCounter",
        "C++ documentation: "
        ":cpp:class:`alpaqa::OCPEvalCounter`\n\n");
    py::class_<alpaqa::OCPEvalCounter::OCPEvalTimer>(
        ocpevalcounter, "OCPEvalTimer",
        "C++ documentation: "
        ":cpp:class:`alpaqa::OCPEvalCounter::OCPEvalTimer`\n\n")
        .def(py::pickle(
            [](const alpaqa::OCPEvalCounter::OCPEvalTimer &p) { // __getstate__
                return py::make_tuple(p.f, p.jac_f, p.grad_f_prod, p.h, p.h_N, p.l, p.l_N, p.qr,
                                      p.q_N, p.add_Q, p.add_Q_N, p.add_R_masked, p.add_S_masked,
                                      p.add_R_prod_masked, p.add_S_prod_masked, p.constr,
                                      p.constr_N, p.grad_constr_prod, p.grad_constr_prod_N,
                                      p.add_gn_hess_constr, p.add_gn_hess_constr_N);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 21)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::OCPEvalCounter::OCPEvalTimer;
                return T{
                    py::cast<decltype(T::f)>(t[0]),
                    py::cast<decltype(T::jac_f)>(t[1]),
                    py::cast<decltype(T::grad_f_prod)>(t[2]),
                    py::cast<decltype(T::h)>(t[3]),
                    py::cast<decltype(T::h_N)>(t[4]),
                    py::cast<decltype(T::l)>(t[5]),
                    py::cast<decltype(T::l_N)>(t[6]),
                    py::cast<decltype(T::qr)>(t[7]),
                    py::cast<decltype(T::q_N)>(t[8]),
                    py::cast<decltype(T::add_Q)>(t[9]),
                    py::cast<decltype(T::add_Q_N)>(t[10]),
                    py::cast<decltype(T::add_R_masked)>(t[11]),
                    py::cast<decltype(T::add_S_masked)>(t[12]),
                    py::cast<decltype(T::add_R_prod_masked)>(t[13]),
                    py::cast<decltype(T::add_S_prod_masked)>(t[14]),
                    py::cast<decltype(T::constr)>(t[15]),
                    py::cast<decltype(T::constr_N)>(t[16]),
                    py::cast<decltype(T::grad_constr_prod)>(t[17]),
                    py::cast<decltype(T::grad_constr_prod_N)>(t[18]),
                    py::cast<decltype(T::add_gn_hess_constr)>(t[19]),
                    py::cast<decltype(T::add_gn_hess_constr_N)>(t[20]),
                };
            }))
        // clang-format off
        .def_readwrite("f", &alpaqa::OCPEvalCounter::OCPEvalTimer::f)
        .def_readwrite("jac_f", &alpaqa::OCPEvalCounter::OCPEvalTimer::jac_f)
        .def_readwrite("grad_f_prod", &alpaqa::OCPEvalCounter::OCPEvalTimer::grad_f_prod)
        .def_readwrite("h", &alpaqa::OCPEvalCounter::OCPEvalTimer::h)
        .def_readwrite("h_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::h_N)
        .def_readwrite("l", &alpaqa::OCPEvalCounter::OCPEvalTimer::l)
        .def_readwrite("l_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::l_N)
        .def_readwrite("qr", &alpaqa::OCPEvalCounter::OCPEvalTimer::qr)
        .def_readwrite("q_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::q_N)
        .def_readwrite("add_Q", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_Q)
        .def_readwrite("add_Q_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_Q_N)
        .def_readwrite("add_R_masked", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_R_masked)
        .def_readwrite("add_S_masked", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_S_masked)
        .def_readwrite("add_R_prod_masked", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_R_prod_masked)
        .def_readwrite("add_S_prod_masked", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_S_prod_masked)
        .def_readwrite("constr", &alpaqa::OCPEvalCounter::OCPEvalTimer::constr)
        .def_readwrite("constr_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::constr_N)
        .def_readwrite("grad_constr_prod", &alpaqa::OCPEvalCounter::OCPEvalTimer::grad_constr_prod)
        .def_readwrite("grad_constr_prod_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::grad_constr_prod_N)
        .def_readwrite("add_gn_hess_constr", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_gn_hess_constr)
        .def_readwrite("add_gn_hess_constr_N", &alpaqa::OCPEvalCounter::OCPEvalTimer::add_gn_hess_constr_N)
        // clang-format on
        ;

    ocpevalcounter
        .def(py::pickle(
            [](const alpaqa::OCPEvalCounter &p) { // __getstate__
                return py::make_tuple(p.f, p.jac_f, p.grad_f_prod, p.h, p.h_N, p.l, p.l_N, p.qr,
                                      p.q_N, p.add_Q, p.add_Q_N, p.add_R_masked, p.add_S_masked,
                                      p.add_R_prod_masked, p.add_S_prod_masked, p.constr,
                                      p.constr_N, p.grad_constr_prod, p.grad_constr_prod_N,
                                      p.add_gn_hess_constr, p.add_gn_hess_constr_N, p.time);
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 22)
                    throw std::runtime_error("Invalid state!");
                using T = alpaqa::OCPEvalCounter;
                return T{
                    py::cast<decltype(T::f)>(t[0]),
                    py::cast<decltype(T::jac_f)>(t[1]),
                    py::cast<decltype(T::grad_f_prod)>(t[2]),
                    py::cast<decltype(T::h)>(t[3]),
                    py::cast<decltype(T::h_N)>(t[4]),
                    py::cast<decltype(T::l)>(t[5]),
                    py::cast<decltype(T::l_N)>(t[6]),
                    py::cast<decltype(T::qr)>(t[7]),
                    py::cast<decltype(T::q_N)>(t[8]),
                    py::cast<decltype(T::add_Q)>(t[9]),
                    py::cast<decltype(T::add_Q_N)>(t[10]),
                    py::cast<decltype(T::add_R_masked)>(t[11]),
                    py::cast<decltype(T::add_S_masked)>(t[12]),
                    py::cast<decltype(T::add_R_prod_masked)>(t[13]),
                    py::cast<decltype(T::add_S_prod_masked)>(t[14]),
                    py::cast<decltype(T::constr)>(t[15]),
                    py::cast<decltype(T::constr_N)>(t[16]),
                    py::cast<decltype(T::grad_constr_prod)>(t[17]),
                    py::cast<decltype(T::grad_constr_prod_N)>(t[18]),
                    py::cast<decltype(T::add_gn_hess_constr)>(t[19]),
                    py::cast<decltype(T::add_gn_hess_constr_N)>(t[20]),
                    py::cast<decltype(T::time)>(t[21]),
                };
            }))
        // clang-format off
        .def_readwrite("f", &alpaqa::OCPEvalCounter::f)
        .def_readwrite("jac_f", &alpaqa::OCPEvalCounter::jac_f)
        .def_readwrite("grad_f_prod", &alpaqa::OCPEvalCounter::grad_f_prod)
        .def_readwrite("h", &alpaqa::OCPEvalCounter::h)
        .def_readwrite("h_N", &alpaqa::OCPEvalCounter::h_N)
        .def_readwrite("l", &alpaqa::OCPEvalCounter::l)
        .def_readwrite("l_N", &alpaqa::OCPEvalCounter::l_N)
        .def_readwrite("qr", &alpaqa::OCPEvalCounter::qr)
        .def_readwrite("q_N", &alpaqa::OCPEvalCounter::q_N)
        .def_readwrite("add_Q", &alpaqa::OCPEvalCounter::add_Q)
        .def_readwrite("add_Q_N", &alpaqa::OCPEvalCounter::add_Q_N)
        .def_readwrite("add_R_masked", &alpaqa::OCPEvalCounter::add_R_masked)
        .def_readwrite("add_S_masked", &alpaqa::OCPEvalCounter::add_S_masked)
        .def_readwrite("add_R_prod_masked", &alpaqa::OCPEvalCounter::add_R_prod_masked)
        .def_readwrite("add_S_prod_masked", &alpaqa::OCPEvalCounter::add_S_prod_masked)
        .def_readwrite("constr", &alpaqa::OCPEvalCounter::constr)
        .def_readwrite("constr_N", &alpaqa::OCPEvalCounter::constr_N)
        .def_readwrite("grad_constr_prod", &alpaqa::OCPEvalCounter::grad_constr_prod)
        .def_readwrite("grad_constr_prod_N", &alpaqa::OCPEvalCounter::grad_constr_prod_N)
        .def_readwrite("add_gn_hess_constr", &alpaqa::OCPEvalCounter::add_gn_hess_constr)
        .def_readwrite("add_gn_hess_constr_N", &alpaqa::OCPEvalCounter::add_gn_hess_constr_N)
        .def_readwrite("time", &alpaqa::OCPEvalCounter::time)
        // clang-format on
        .def("__str__", [](const alpaqa::OCPEvalCounter &c) {
            std::ostringstream os;
            os << c;
            return os.str();
        });
}
