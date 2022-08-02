#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/inner/directions/panoc-direction-update.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/util/type-erasure.hpp>

#include <new>
#include <string>
#include <utility>

namespace alpaqa {

template <Config Conf>
struct PANOCDirectionVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);

    void (*initialize)(void *self, crvec x_0, crvec x̂_0, crvec p_0, crvec grad_0)        = nullptr;
    bool (*update)(void *self, crvec xₖ, crvec xₙₑₓₜ, crvec pₖ, crvec pₙₑₓₜ, crvec grad_new,
                   const Box<config_t> &C, real_t γ_new)                             = nullptr;
    bool (*apply)(const void *self, crvec xₖ, crvec x̂ₖ, crvec pₖ, real_t γ, rvec qₖ) = nullptr;
    void (*changed_γ)(void *self, real_t γₖ, real_t old_γₖ)                          = nullptr;
    void (*reset)(void *self)                                                        = nullptr;
    std::string (*get_name)(const void *self)                                        = nullptr;

    template <class T>
    PANOCDirectionVTable(util::VTableTypeTag<T> t) : util::BasicVTable{t} {
        initialize = util::type_erased_wrapped<&T::initialize>();
        update     = util::type_erased_wrapped<&T::update>();
        apply      = util::type_erased_wrapped<&T::apply>();
        changed_γ  = util::type_erased_wrapped<&T::changed_γ>();
        reset      = util::type_erased_wrapped<&T::reset>();
        get_name   = util::type_erased_wrapped<&T::get_name>();
    }
    PANOCDirectionVTable() = default;
};

template <Config Conf>
constexpr size_t te_pd_buff_size = util::required_te_buffer_size_for<PANOCDirection<LBFGS<Conf>>>();

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedPANOCDirection
    : public util::TypeErased<PANOCDirectionVTable<Conf>, Allocator, te_pd_buff_size<Conf>> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = PANOCDirectionVTable<Conf>;
    using allocator_type = Allocator;
    using TypeErased     = util::TypeErased<VTable, allocator_type, te_pd_buff_size<Conf>>;
    using TypeErased::TypeErased;

  private:
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TypeErasedPANOCDirection make(Args &&...args) {
        return TypeErased::template make<TypeErasedPANOCDirection, T>(std::forward<Args>(args)...);
    }

    template <class... Args>
    decltype(auto) initialize(Args &&...args) {
        return vtable.initialize(self, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) update(Args &&...args) {
        return vtable.update(self, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) apply(Args &&...args) const {
        return vtable.apply(self, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) changed_γ(Args &&...args) {
        return vtable.changed_γ(self, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) reset(Args &&...args) {
        return vtable.reset(self, std::forward<Args>(args)...);
    }
    template <class... Args>
    decltype(auto) get_name(Args &&...args) const {
        return vtable.get_name(self, std::forward<Args>(args)...);
    }
};

template <Config Conf>
struct PANOCDirection<TypeErasedPANOCDirection<Conf>> : TypeErasedPANOCDirection<Conf> {
    PANOCDirection(TypeErasedPANOCDirection<Conf> &&o)
        : TypeErasedPANOCDirection<Conf>(std::move(o)) {}
};

template <class T, class... Args>
auto erase_direction(Args &&...args) {
    return TypeErasedPANOCDirection<typename T::config_t>::template make<PANOCDirection<T>>(
        std::forward<Args>(args)...);
}

} // namespace alpaqa