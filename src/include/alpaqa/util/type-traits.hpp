#pragma once

#include <concepts>
#include <memory>
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

template <class First, class... Pack>
struct last_type {
    using type = typename last_type<Pack...>::type;
};
template <class Only>
struct last_type<Only> {
    using type = Only;
};
template <class... Pack>
using last_type_t = typename last_type<Pack...>::type;

template <class... Pack>
struct first_type_or_void;
template <class First, class... Pack>
struct first_type_or_void<First, Pack...> {
    using type = First;
};
template <>
struct first_type_or_void<> {
    using type = void;
};
template <class... Pack>
using first_type_or_void_t = typename first_type_or_void<Pack...>::type;

template <class... Pack>
concept no_leading_allocator = !
std::is_same_v<std::remove_cvref_t<first_type_or_void_t<Pack...>>,
               std::allocator_arg_t>;

} // namespace alpaqa::util