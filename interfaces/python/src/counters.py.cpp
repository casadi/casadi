#include <alpaqa/problem/problem-counters.hpp>
#include <pybind11/chrono.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void register_counters(py::module_ &m) {
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
}
