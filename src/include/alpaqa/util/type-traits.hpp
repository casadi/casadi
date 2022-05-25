#pragma once

#include <type_traits>

namespace alpaqa::util {

template <class M>
struct class_from_member_ptr_impl {};

template <class C, class Ret, class... Args>
struct class_from_member_ptr_impl<Ret (C::*)(Args...)> {
    using type = C;
};

template <class C, class Ret, class... Args>
struct class_from_member_ptr_impl<Ret (C::*)(Args...) const> {
    using type = std::add_const_t<C>;
};

template <class M>
using class_from_member_ptr_impl_t =
    typename class_from_member_ptr_impl<M>::type;

template <auto M>
using class_from_member_ptr_t = class_from_member_ptr_impl_t<decltype(M)>;

} // namespace alpaqa::util