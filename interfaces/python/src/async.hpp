#pragma once

#include <pybind11/gil.h>
namespace py = pybind11;

#include <chrono>
#include <exception>
#include <future>
#include <tuple>
#include <utility>
using namespace std::chrono_literals;

#include "stream-replacer.hpp"
#include "thread-checker.hpp"

template <class Solver, class Invoker, class... CheckedArgs>
auto async_solve(bool async, Solver &solver, Invoker &invoke_solver, CheckedArgs &...checked_args) {
    if (!async) {
        // Replace the output stream
        StreamReplacer stream{&solver};
        // Invoke the solver synchronously
        auto &&stats = invoke_solver();
        return stats;
    } else {
        // Check that the user doesn't use the same solver/problem in multiple threads
        ThreadChecker solver_checker{&solver};
        auto checkers = std::make_tuple(ThreadChecker{&checked_args}...);
        // Replace the output stream
        StreamReplacer stream{&solver};
        // Invoke the solver asynchronously
        auto stats = std::async(std::launch::async, invoke_solver);
        {
            py::gil_scoped_release gil;
            while (stats.wait_for(50ms) != std::future_status::ready) {
                py::gil_scoped_acquire gil;
                // Check if Python received a signal (e.g. Ctrl+C)
                if (PyErr_CheckSignals() != 0) {
                    // Nicely ask the solver to stop
                    solver.stop();
                    // It should return a result soon
                    if (py::gil_scoped_release gil;
                        stats.wait_for(15s) != std::future_status::ready) {
                        // If it doesn't, we terminate the entire program,
                        // because the solver uses variables local to this
                        // function, so we cannot safely return without
                        // waiting for the solver to finish.
                        std::terminate();
                    }
                    if (PyErr_Occurred())
                        throw py::error_already_set();
                    break;
                }
            }
        }
        return std::move(stats.get());
    }
}