#pragma once

#include <alpaqa/inner/inner-solve-options.hpp>
#include <alpaqa/problem/type-erased-problem.hpp>
#include "type-erased-solver-stats.hpp"

namespace alpaqa {

template <Config Conf>
struct InnerSolverVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Stats        = TypeErasedInnerSolverStats<Conf>;
    using Problem      = TypeErasedProblem<Conf>;
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

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedInnerSolver : public util::TypeErased<InnerSolverVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = InnerSolverVTable<Conf>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using Stats          = typename VTable::Stats;
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
    template <class... Args>
    decltype(auto) stop(Args &&...args) {
        return call(vtable.stop, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) get_name(Args &&...args) const {
        return call(vtable.get_name, std::forward<Args>(args)...);
    }
};

} // namespace alpaqa