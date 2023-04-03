#pragma once

#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/problem/ocproblem.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#include <functional>
#include <optional>
#include <stdexcept>
#include <variant>

#include <dict/stats-to-dict.hpp>
#include <inner/type-erased-inner-solver.hpp>
#include <inner/type-erased-solver-stats.hpp>
#include <util/async.hpp>

#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace alpaqa {

/**
 * Type erasing the @ref ALMSolver class, so that it can be instantiated with
 * inner solvers of different problem types (e.g. @ref TypeErasedProblem and
 * @ref TypeErasedControlProblem).
 * 
 * To this end, it represents the possible problems as a variant, and then
 * std::visits the variant in the @ref call method, throwing an exception if
 * the given problem type does not match the inner solver's problem type.
 */
template <Config Conf>
struct ALMSolverVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Stats = typename ALMSolver<TypeErasedInnerSolver<Conf, TypeErasedProblem<Conf>>>::Stats;
    using Problem =
        std::variant<const TypeErasedProblem<Conf> *, const TypeErasedControlProblem<Conf> *>;

    // clang-format off
    required_function_t<py::tuple(const Problem &, std::optional<vec> x, std::optional<vec> y, bool async)>
        call = nullptr;
    required_function_t<void()>
        stop = nullptr;
    required_const_function_t<std::string()>
        get_name = nullptr;
    required_const_function_t<py::object()>
        get_params = nullptr;
    required_const_function_t<py::object()>
        get_inner_solver = nullptr;
    // clang-format on

    template <class T>
    ALMSolverVTable(util::VTableTypeTag<T> t) : util::BasicVTable{t} {
        stop       = util::type_erased_wrapped<&T::stop>();
        get_name   = util::type_erased_wrapped<&T::get_name>();
        get_params = [](const void *self_) -> py::object {
            auto &self = *std::launder(reinterpret_cast<const T *>(self_));
            return py::cast(self.get_params());
        };
        get_inner_solver = [](const void *self_) -> py::object {
            auto &self = *std::launder(reinterpret_cast<const T *>(self_));
            return py::cast(self.inner_solver);
        };
        call = [](void *self_, const Problem &p, std::optional<vec> x, std::optional<vec> y,
                  bool async) {
            auto &self       = *std::launder(reinterpret_cast<T *>(self_));
            auto call_solver = [&]<class P>(const P *p) -> py::tuple {
                if constexpr (!std::is_same_v<P, typename T::Problem>)
                    throw std::invalid_argument(
                        "Unsupported problem type (got '" + demangled_typename(typeid(P)) +
                        "', expected '" + demangled_typename(typeid(typename T::Problem)) + "')");
                else
                    return safe_call_solver(self, p, x, y, async);
            };
            return std::visit(call_solver, p);
        };
    }
    ALMSolverVTable() = default;

    template <class T>
    static decltype(auto) safe_call_solver(T &self, const auto &p, std::optional<vec> &x,
                                           std::optional<vec> &y, bool async) {
        using InnerSolver = typename T::InnerSolver;
        alpaqa::util::check_dim_msg<config_t>(x, p->get_n(),
                                              "Length of x does not match problem size problem.n");
        alpaqa::util::check_dim_msg<config_t>(y, p->get_m(),
                                              "Length of y does not match problem size problem.m");
        auto invoke_solver = [&] { return self(*p, *x, *y); };
        auto stats         = async_solve(async, self, invoke_solver, *p);
        return py::make_tuple(std::move(*x), std::move(*y),
                              alpaqa::conv::stats_to_dict<InnerSolver>(std::move(stats)));
    }
};

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedALMSolver : public util::TypeErased<ALMSolverVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = ALMSolverVTable<Conf>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using Stats          = typename VTable::Stats;
    using Problem        = typename VTable::Problem;
    using TypeErased::TypeErased;

  protected:
    using TypeErased::call;
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TypeErasedALMSolver make(Args &&...args) {
        return TypeErased::template make<TypeErasedALMSolver, T>(std::forward<Args>(args)...);
    }

    decltype(auto) operator()(const Problem &p, std::optional<vec> x, std::optional<vec> y,
                              bool async) {
        return call(vtable.call, p, x, y, async);
    }
    decltype(auto) stop() { return call(vtable.stop); }
    decltype(auto) get_name() const { return call(vtable.get_name); }
    decltype(auto) get_params() const { return call(vtable.get_params); }
    decltype(auto) get_inner_solver() const { return call(vtable.get_inner_solver); }
};

} // namespace alpaqa