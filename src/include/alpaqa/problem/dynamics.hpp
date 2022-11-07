#pragma once

#include <alpaqa/config/config.hpp>
#include <alpaqa/problem/box.hpp>
#include <alpaqa/problem/structure.hpp>
#include <alpaqa/util/type-erasure.hpp>

namespace alpaqa {

template <Config Conf>
struct ControlProblem {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    ControlProblem(length_t N, length_t nx, length_t nu) : N{N}, nx{nx}, nu{nu} {}

    void get_U(Box &U) const;
    void get_x_init(rvec x_init) const;
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    real_t eval_l(index_t timestep, crvec h) const;
    real_t eval_l_N(crvec h) const;
    void eval_grad_l(index_t timestep, crvec h, rvec grad_lh) const;
    void eval_grad_l_N(crvec h, rvec grad_lh) const;
    void eval_hess_l(index_t timestep, crvec h, rmat hess_lh) const;
    void eval_hess_l_N(crvec h, rmat hess_lh) const;
    CostStructure get_l_structure() const;
    void check() const;

    length_t N, nx, nu;
};

template <Config Conf>
struct ControlProblemVTable : util::BasicVTable {
    USING_ALPAQA_CONFIG(Conf);
    using Box = alpaqa::Box<config_t>;

    // clang-format off
    required_const_function_t<void(Box &U)>
        get_U;
    required_const_function_t<void(rvec x_init)>
        get_x_init;
    required_const_function_t<void(index_t timestep, crvec x, crvec u, rvec fxu)>
        eval_f;
    required_const_function_t<void(index_t timestep, crvec x, crvec u, rmat J_fxu)>
        eval_jac_f;
    required_const_function_t<real_t(index_t timestep, crvec h)>
        eval_l;
    required_const_function_t<real_t(crvec h)>
        eval_l_N;
    required_const_function_t<void(index_t timestep, crvec h, rvec grad_lh)>
        eval_grad_l;
    required_const_function_t<void(crvec h, rvec grad_lh)>
        eval_grad_l_N;
    required_const_function_t<void(index_t timestep, crvec h, rmat hess_lh)>
        eval_hess_l;
    required_const_function_t<void(crvec h, rmat hess_lh)>
        eval_hess_l_N;
    required_const_function_t<CostStructure()>
        get_l_structure;
    required_const_function_t<void()>
        check;
    // clang-format on

    length_t N, nu, nx, nh;

    template <class P>
    ControlProblemVTable(util::VTableTypeTag<P> t) : util::BasicVTable{t} {
        get_U           = util::type_erased_wrapped<&P::get_U>();
        get_x_init      = util::type_erased_wrapped<&P::get_x_init>();
        eval_f          = util::type_erased_wrapped<&P::eval_f>();
        eval_jac_f      = util::type_erased_wrapped<&P::eval_jac_f>();
        eval_l          = util::type_erased_wrapped<&P::eval_l>();
        eval_l_N        = util::type_erased_wrapped<&P::eval_l_N>();
        eval_grad_l     = util::type_erased_wrapped<&P::eval_grad_l>();
        eval_grad_l_N   = util::type_erased_wrapped<&P::eval_grad_l_N>();
        eval_hess_l     = util::type_erased_wrapped<&P::eval_hess_l>();
        eval_hess_l_N   = util::type_erased_wrapped<&P::eval_hess_l_N>();
        get_l_structure = util::type_erased_wrapped<&P::get_l_structure>();
        check           = util::type_erased_wrapped<&P::check>();
        N               = t.t->get_N();
        nu              = t.t->get_nu();
        nx              = t.t->get_nx();
        nh              = t.t->get_nh();
    }
    ControlProblemVTable() = default;
};

template <Config Conf = DefaultConfig, class Allocator = std::allocator<std::byte>>
class TypeErasedControlProblem : public util::TypeErased<ControlProblemVTable<Conf>, Allocator> {
  public:
    USING_ALPAQA_CONFIG(Conf);
    using VTable         = ControlProblemVTable<config_t>;
    using allocator_type = Allocator;
    using Box            = typename VTable::Box;
    using TypeErased     = util::TypeErased<VTable, allocator_type>;
    using TypeErased::TypeErased;

  protected:
    using TypeErased::call;
    using TypeErased::self;
    using TypeErased::vtable;

  public:
    template <class T, class... Args>
    static TypeErasedControlProblem make(Args &&...args) {
        return TypeErased::template make<TypeErasedControlProblem, T>(std::forward<Args>(args)...);
    }

    /// @name Problem dimensions
    /// @{

    /// Horizon length.
    length_t get_N() const { return vtable.N; }
    /// Number of inputs.
    length_t get_nu() const { return vtable.nu; }
    /// Number of states.
    length_t get_nx() const { return vtable.nx; }
    /// Number of outputs.
    length_t get_ny() const { return vtable.ny; }

    /// @}

    void get_U(Box &U) const;
    void get_x_init(rvec x_init) const;
    void eval_f(index_t timestep, crvec x, crvec u, rvec fxu) const;
    void eval_jac_f(index_t timestep, crvec x, crvec u, rmat J_fxu) const;
    real_t eval_l(index_t timestep, crvec h) const;
    real_t eval_l_N(crvec h) const;
    void eval_grad_l(index_t timestep, crvec h, rvec grad_lh) const;
    void eval_grad_l_N(crvec h, rvec grad_lh) const;
    void eval_hess_l(index_t timestep, crvec h, rmat hess_lh) const;
    void eval_hess_l_N(crvec h, rmat hess_lh) const;
    CostStructure get_l_structure() const;
    void check() const;
};

template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::get_U(Box &U) const {
    return call(vtable.get_U, U);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::get_x_init(rvec x_init) const {
    return call(vtable.get_x_init, x_init);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_f(index_t timestep, crvec x, crvec u,
                                                       rvec fxu) const {
    return call(vtable.eval_f, timestep, x, u, fxu);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_jac_f(index_t timestep, crvec x, crvec u,
                                                           rmat J_fxu) const {
    return call(vtable.eval_jac_f, timestep, x, u, J_fxu);
}
template <Config Conf, class Allocator>
auto TypeErasedControlProblem<Conf, Allocator>::eval_l(index_t timestep, crvec h) const -> real_t {
    return call(vtable.eval_l, timestep, h);
}
template <Config Conf, class Allocator>
auto TypeErasedControlProblem<Conf, Allocator>::eval_l_N(crvec h) const -> real_t {
    return call(vtable.eval_l_N, h);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_grad_l(index_t timestep, crvec h,
                                                            rvec grad_lh) const {
    return call(vtable.eval_grad_l, timestep, h, grad_lh);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_grad_l_N(crvec h, rvec grad_lh) const {
    return call(vtable.eval_grad_l_N, h, grad_lh);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_hess_l(index_t timestep, crvec h,
                                                            rmat hess_lh) const {
    return call(vtable.eval_hess_l, timestep, h, hess_lh);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::eval_hess_l_N(crvec h, rmat hess_lh) const {
    return call(vtable.eval_hess_l_N, h, hess_lh);
}
template <Config Conf, class Allocator>
CostStructure TypeErasedControlProblem<Conf, Allocator>::get_l_structure() const {
    return call(vtable.get_l_structure);
}
template <Config Conf, class Allocator>
void TypeErasedControlProblem<Conf, Allocator>::check() const {
    return call(vtable.check);
}

} // namespace alpaqa