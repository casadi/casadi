#pragma once

#include <alpaqa/util/type-traits.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#ifndef NDEBUG
#include <typeinfo>
#endif

namespace alpaqa::util {

template <class T>
struct VTableTypeTag {
    T *t = nullptr;
};

/// Struct that stores the size of a polymorphic object, as well as pointers to
/// functions to copy, move or destroy the object.
/// Inherit from this struct to add useful functions.
struct BasicVTable {

    template <class>
    struct required_function;
    template <class R, class... Args>
    struct required_function<R(Args...)> {
        using type = R (*)(void *self, Args...);
    };
    template <class>
    struct required_const_function;
    template <class R, class... Args>
    struct required_const_function<R(Args...)> {
        using type = R (*)(const void *self, Args...);
    };
    template <class F>
    using required_function_t = typename required_function<F>::type;
    template <class F>
    using required_const_function_t = typename required_const_function<F>::type;

    template <class, class VTable = BasicVTable>
    struct optional_function;
    template <class R, class... Args, class VTable>
    struct optional_function<R(Args...), VTable> {
        using type = R (*)(void *self, Args..., const VTable &);
    };
    template <class, class VTable = BasicVTable>
    struct optional_const_function;
    template <class R, class... Args, class VTable>
    struct optional_const_function<R(Args...), VTable> {
        using type = R (*)(const void *self, Args..., const VTable &);
    };
    template <class F, class VTable = BasicVTable>
    using optional_function_t = typename optional_function<F, VTable>::type;
    template <class F, class VTable = BasicVTable>
    using optional_const_function_t =
        typename optional_const_function<F, VTable>::type;

    size_t size = static_cast<size_t>(0xDEADBEEFDEADBEEF);
    required_const_function_t<void(void *storage)> copy = nullptr;
    required_function_t<void(void *storage)> move       = nullptr;
    required_function_t<void()> destroy                 = nullptr;
#ifndef NDEBUG
    const std::type_info *type = &typeid(void);
#endif

    BasicVTable() = default;

    template <class T>
    BasicVTable(VTableTypeTag<T>) {
        size = sizeof(T);
        copy = [](const void *self, void *storage) {
            new (storage) T(*std::launder(reinterpret_cast<const T *>(self)));
        };
        move = [](void *self, void *storage) {
            new (storage)
                T(std::move(*std::launder(reinterpret_cast<T *>(self))));
        };
        destroy = [](void *self) {
            std::launder(reinterpret_cast<T *>(self))->~T();
        };
#ifndef NDEBUG
        type = &typeid(T);
#endif
    }
};

namespace detail {
template <class... ExtraArgs>
struct Launderer {
  private:
    template <auto M, class V, class T, class R, class... Args>
    static constexpr auto do_invoke(V *self, Args... args, ExtraArgs...) -> R {
        return std::invoke(M, *std::launder(reinterpret_cast<T *>(self)),
                           std::forward<Args>(args)...);
    }
    template <auto M, class T, class R, class... Args>
    static constexpr auto invoke_ovl(R (T::*)(Args...) const) {
        return do_invoke<M, const void, const T, R, Args...>;
    }
    template <auto M, class T, class R, class... Args>
    static constexpr auto invoke_ovl(R (T::*)(Args...)) {
        return do_invoke<M, void, T, R, Args...>;
    }

  public:
    template <auto M>
    static constexpr auto invoke() {
        return invoke_ovl<M>(M);
    }
};
} // namespace detail

template <auto Method, class... ExtraArgs>
constexpr auto type_erased_wrapped() {
    return detail::Launderer<ExtraArgs...>::template invoke<Method>();
}

template <class VTable, class Allocator>
inline constexpr size_t default_te_buffer_size() {
    struct S {
        [[no_unique_address]] Allocator allocator;
        void *self = nullptr;
        VTable vtable;
    };
    const size_t max_size = 128;
    return max_size - std::min(max_size, sizeof(S));
}

template <class... Types>
inline constexpr size_t required_te_buffer_size_for() {
    constexpr size_t sizes[] = {sizeof(Types)...};
    return *std::max_element(std::begin(sizes), std::end(sizes));
}

/// Class for polymorphism through type erasure. Saves the entire vtable, and
/// uses small buffer optimization.
template <class VTable           = BasicVTable,
          class Allocator        = std::allocator<std::byte>,
          size_t SmallBufferSize = default_te_buffer_size<VTable, Allocator>()>
class TypeErased {
  public:
    static constexpr size_t small_buffer_size = SmallBufferSize;
    using allocator_type                      = Allocator;

  private:
    struct EmptyBuffer {
        static constexpr std::byte *data() { return nullptr; }
    };
    using allocator_traits  = std::allocator_traits<allocator_type>;
    using buffer_array_type = std::array<std::byte, small_buffer_size>;
    using buffer_type = std::conditional_t<small_buffer_size == 0, EmptyBuffer,
                                           buffer_array_type>;
    [[no_unique_address]] alignas(std::max_align_t) buffer_type small_buffer;
    [[no_unique_address]] allocator_type allocator;

  protected:
    void *self = nullptr;
    VTable vtable;

  public:
    /// Default constructor.
    TypeErased() = default;
    /// Default constructor (allocator aware).
    TypeErased(std::allocator_arg_t, allocator_type alloc)
        : allocator{std::move(alloc)} {}
    /// Copy constructor.
    TypeErased(const TypeErased &other)
        : allocator{allocator_traits::select_on_container_copy_construction(
              other.allocator)} {
        do_copy_assign<false>(other);
    }
    /// Copy constructor (allocator aware).
    TypeErased(const TypeErased &other, allocator_type alloc)
        : allocator{std::move(alloc)} {
        do_copy_assign<false>(other);
    }
    /// Copy assignment.
    TypeErased &operator=(const TypeErased &other) {
        if (&other == this)
            return *this;
        cleanup();
        do_copy_assign<true>(other);
        return *this;
    }

    /// Move constructor.
    TypeErased(TypeErased &&other) : allocator{std::move(other.allocator)} {
        do_move_assign<false>(std::move(other));
    }
    /// Move constructor (allocator aware).
    TypeErased(TypeErased &&other, allocator_type alloc)
        : allocator{std::move(alloc)} {
        do_move_assign<false>(std::move(other));
    }
    /// Move assignment.
    TypeErased &operator=(TypeErased &&other) {
        if (&other == this)
            return *this;
        cleanup();
        do_move_assign<true>(std::move(other));
        return *this;
    }

    /// Destructor.
    ~TypeErased() { cleanup(); }

    /// Main constructor that type-erases the given argument. Requirement
    /// prevents this constructor from taking precedence over the copy and move
    /// constructors.
    template <class Tref>
    explicit TypeErased(Tref &&d, allocator_type alloc = {}) //
        requires(!std::is_base_of_v<TypeErased, std::remove_cvref_t<Tref>>)
        : allocator{std::move(alloc)} {
        using T = std::remove_cvref_t<Tref>;
        construct_inplace<T>(std::forward<Tref>(d));
    }

    /// Construct a type-erased wrapper of type Ret for an object of type T,
    /// initialized in-place with the given arguments.
    template <class Ret, class T, class... Args>
    requires std::is_base_of_v<TypeErased, Ret>
    static Ret make(Args &&...args) {
        Ret r{};
        r.template construct_inplace<T>(std::forward<Args>(args)...);
        return r;
    }
    /// Construct a type-erased wrapper of type Ret for an object of type T,
    /// initialized in-place with the given arguments.
    template <class Ret, class T, class... Args>
    requires std::is_base_of_v<TypeErased, Ret>
    static Ret make(std::allocator_arg_t, allocator_type alloc,
                    Args &&...args) {
        Ret r{std::move(alloc)};
        r.template construct_inplace<T>(std::forward<Args>(args)...);
        return r;
    }

    /// Check if this wrapper wraps an object. False for default-constructed
    /// objects.
    explicit operator bool() const { return self; }

    /// Get a copy of the allocator.
    allocator_type get_allocator() const { return allocator; }

    /// Convert the type-erased object to the given type. The type is checked
    /// in debug builds only, use with caution.
    template <class T>
    T &as() {
#ifndef NDEBUG
        assert(typeid(T) == *vtable.type);
#endif
        return *reinterpret_cast<T *>(self);
    }
    /// @copydoc as()
    template <class T>
    const T &as() const {
#ifndef NDEBUG
        assert(typeid(T) == *vtable.type);
#endif
        return *reinterpret_cast<const T *>(self);
    }

  private:
    /// Ensure that storage is available, either by using the small buffer if
    /// it is large enough, or by calling the allocator.
    void ensure_storage(size_t size) {
        assert(!self);
        self = size <= small_buffer_size ? small_buffer.data()
                                         : allocator.allocate(size);
    }
    void ensure_storage() { ensure_storage(vtable.size); }

    /// Destroy the type-erased object (if not empty), and deallocate the memory
    /// if necessary.
    void cleanup() {
        if (self) {
            vtable.destroy(self);
            if (vtable.size > small_buffer_size)
                allocator.deallocate(reinterpret_cast<std::byte *>(self),
                                     vtable.size);
            self = nullptr;
        }
    }

    template <bool CopyAllocator>
    void do_copy_assign(const TypeErased &other) {
        constexpr bool prop_alloc =
            allocator_traits::propagate_on_container_copy_assignment::value;
        if constexpr (CopyAllocator && prop_alloc)
            allocator = other.allocator;
        vtable = other.vtable;
        ensure_storage();
        vtable.copy(other.self, self);
    }

    template <bool MoveAllocator>
    void do_move_assign(TypeErased &&other) {
        constexpr bool prop_alloc =
            allocator_traits::propagate_on_container_move_assignment::value;
        if constexpr (MoveAllocator && prop_alloc)
            allocator = std::move(other.allocator);
        vtable = other.vtable;
        // If dynamically allocated, simply steal storage
        if (vtable.size > small_buffer_size) {
            // If we can steal the storage because of equal allocators, do so
            if (allocator == other.allocator) {
                self = std::exchange(other.self, nullptr);
            }
            // If the allocators are not the same, we cannot steal the storage,
            // so do an explicit move
            else {
                self = allocator.allocate(vtable.size);
                vtable.move(other.self, self);
                other.cleanup();
            }
        }
        // Otherwise, use the small buffer and do an explicit move
        else if (other.self) {
            self = small_buffer.data();
            vtable.move(other.self, self);
            other.cleanup();
        }
    }

  protected:
    /// Ensure storage and construct the type-erased object of type T in-place.
    template <class T, class... Args>
    void construct_inplace(Args &&...args) {
        ensure_storage(sizeof(T));
        T *t   = new (self) T{std::forward<Args>(args)...};
        vtable = VTable{VTableTypeTag<T>{t}};
    }

    /// Call the vtable function @p f with the given arguments @p args,
    /// implicitly passing the @ref self pointer and @ref vtable reference if
    /// necessary.
    template <class Ret, class... FArgs, class... Args>
    decltype(auto) call(Ret (*f)(const void *, FArgs...),
                        Args &&...args) const {
        assert(f);
        assert(self);
        using LastArg = util::last_type_t<FArgs...>;
        if constexpr (std::is_same_v<LastArg, const VTable &>)
            return f(self, std::forward<Args>(args)..., vtable);
        else
            return f(self, std::forward<Args>(args)...);
    }
    /// @copydoc call
    template <class Ret, class... FArgs, class... Args>
    decltype(auto) call(Ret (*f)(void *, FArgs...), Args &&...args) {
        assert(f);
        assert(self);
        using LastArg = util::last_type_t<FArgs...>;
        if constexpr (std::is_same_v<LastArg, const VTable &>)
            return f(self, std::forward<Args>(args)..., vtable);
        else
            return f(self, std::forward<Args>(args)...);
    }
    /// @copydoc call
    template <class Ret>
    decltype(auto) call(Ret (*f)(const void *)) const {
        assert(f);
        assert(self);
        return f(self);
    }
    /// @copydoc call
    template <class Ret>
    decltype(auto) call(Ret (*f)(void *)) {
        assert(f);
        assert(self);
        return f(self);
    }
};

} // namespace alpaqa::util