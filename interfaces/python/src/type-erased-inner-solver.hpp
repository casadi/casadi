#pragma once

#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include <optional>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "type-erased-solver-stats.hpp"

namespace alpaqa {

template <Config Conf, class ProblemT>
struct InnerSolverVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Stats        = TypeErasedInnerSolverStats<Conf>;
    using Problem      = ProblemT;
    using SolveOptions = InnerSolveOptions<config_t>;

    // clang-format off
    required_function_t<Stats(const Problem &, const SolveOptions &, rvec, rvec, crvec, rvec)>
        call = nullptr;
    required_function_t<void()>
        stop = nullptr;
    required_const_function_t<std::string()>
        get_name = nullptr;
    // clang-format on

    template <class T>
    InnerSolverVTable(util::VTableTypeTag<T> t) : util::BasicVTable{t} {
        stop     = util::type_erased_wrapped<&T::stop>();
        get_name = util::type_erased_wrapped<&T::get_name>();
        call     = []<class... Args>(void *self_, const Problem &p, Args... args) {
            auto &self = *std::launder(reinterpret_cast<T *>(self_));
            return Stats{self(p, std::forward<Args>(args)...)};
        };
    }
    InnerSolverVTable() = default;
};

template <Config Conf, class ProblemT, class Allocator = std::allocator<std::byte>>
class TypeErasedInnerSolver
    : public util::TypeErased<InnerSolverVTable<Conf, ProblemT>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = InnerSolverVTable<Conf, ProblemT>;
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
    static TypeErasedInnerSolver make(Args &&...args) {
        return TypeErased::template make<TypeErasedInnerSolver, T>(std::forward<Args>(args)...);
    }

    template <class... Args>
    decltype(auto) operator()(Args &&...args) {
        return call(vtable.call, std::forward<Args>(args)...);
    }
    decltype(auto) stop() { return call(vtable.stop); }
    decltype(auto) get_name() const { return call(vtable.get_name); }
};

} // namespace alpaqa

template <class InnerSolverT>
struct InnerSolverConversion {
    using InnerSolver = InnerSolverT;
    std::optional<py::class_<InnerSolver>> cls;
    void initialize(py::class_<InnerSolver> &&cls) {
        assert(!this->cls);
        this->cls.emplace(std::move(cls));
    }
    template <class T>
    void implicitly_convertible_to() {
        assert(this->cls);
        cls->def(py::init([](const T &t) { return std::make_unique<InnerSolver>(t); }));
        py::implicitly_convertible<T, InnerSolver>();
    }
};

/// Global instance of the py::class_<InnerSolverT> binding, for registering
/// converting constructors from concrete inner solvers later.
template <class InnerSolverT>
inline InnerSolverConversion<InnerSolverT> inner_solver_class;
