#pragma once

#include <alpaqa/export.hpp>
#include <alpaqa/util/demangled-typename.hpp>
#include <alpaqa/util/noop-delete.hpp>
#include <alpaqa/util/type-traits.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace alpaqa::util {

class ALPAQA_EXPORT bad_type_erased_type : public std::logic_error {
  public:
    bad_type_erased_type(const std::type_info &actual_type,
                         const std::type_info &requested_type,
                         const std::string &message = "")
        : std::logic_error{message}, actual_type{actual_type},
          requested_type{requested_type} {}

    [[nodiscard]] const char *what() const noexcept override {
        message = "";
        if (const char *w = std::logic_error::what(); w && *w) {
            message += w;
            message += ": ";
        }
        message = "Type requested: " + demangled_typename(requested_type) +
                  ", type contained: " + demangled_typename(actual_type);
        return message.c_str();
    }

  private:
    const std::type_info &actual_type;
    const std::type_info &requested_type;
    mutable std::string message;
};

/// Struct that stores the size of a polymorphic object, as well as pointers to
/// functions to copy, move or destroy the object.
/// Inherit from this struct to add useful functions.
struct BasicVTable {

    template <class>
    struct required_function; // undefined
    template <class R, class... Args>
    struct required_function<R(Args...)> {
        using type = R (*)(void *self, Args...);
    };
    template <class>
    struct required_const_function; // undefined
    template <class R, class... Args>
    struct required_const_function<R(Args...)> {
        using type = R (*)(const void *self, Args...);
    };
    template <class, class VTable = BasicVTable>
    struct optional_function; // undefined
    template <class R, class... Args, class VTable>
    struct optional_function<R(Args...), VTable> {
        using type = R (*)(void *self, Args..., const VTable &);
    };
    template <class, class VTable = BasicVTable>
    struct optional_const_function; // undefined
    template <class R, class... Args, class VTable>
    struct optional_const_function<R(Args...), VTable> {
        using type = R (*)(const void *self, Args..., const VTable &);
    };
    /// A required function includes a void pointer to self, in addition to the
    /// arguments of @p F.
    template <class F>
    using required_function_t = typename required_function<F>::type;
    /// @copydoc required_function_t
    /// For const-qualified member functions.
    template <class F>
    using required_const_function_t = typename required_const_function<F>::type;
    /// An optional function includes a void pointer to self, the arguments of
    /// @p F, and an additional reference to the VTable, so that it can be
    /// implemented in terms of other functions.
    template <class F, class VTable = BasicVTable>
    using optional_function_t = typename optional_function<F, VTable>::type;
    /// @copydoc optional_function_t
    /// For const-qualified member functions.
    template <class F, class VTable = BasicVTable>
    using optional_const_function_t =
        typename optional_const_function<F, VTable>::type;

    /// Copy-construct a new instance into storage.
    required_const_function_t<void(void *storage)> copy = nullptr;
    /// Move-construct a new instance into storage.
    required_function_t<void(void *storage)> move = nullptr;
    /// Destruct the given instance.
    required_function_t<void()> destroy = nullptr;
    /// The original type of the stored object.
    const std::type_info *type = &typeid(void);

    BasicVTable() = default;

    template <class T>
    BasicVTable(std::in_place_t, T &) noexcept {
        copy = [](const void *self, void *storage) {
            new (storage) T(*std::launder(reinterpret_cast<const T *>(self)));
        };
        // TODO: require that move constructor is noexcept?
        move = [](void *self, void *storage) noexcept {
            new (storage)
                T(std::move(*std::launder(reinterpret_cast<T *>(self))));
        };
        destroy = [](void *self) {
            std::launder(reinterpret_cast<T *>(self))->~T();
        };
        type = &typeid(T);
    }
};

namespace detail {
template <class Class, class... ExtraArgs>
struct Launderer {
  private:
    template <auto M, class V, class C, class R, class... Args>
    [[gnu::always_inline]] static constexpr auto
    do_invoke(V *self, Args... args, ExtraArgs...) -> R {
        return std::invoke(M, *std::launder(reinterpret_cast<C *>(self)),
                           std::forward<Args>(args)...);
    }
    template <auto M, class T, class R, class... Args>
        requires std::is_base_of_v<T, Class>
    [[gnu::always_inline]] static constexpr auto invoker_ovl(R (T::*)(Args...)
                                                                 const) {
        return do_invoke<M, const void, const Class, R, Args...>;
    }
    template <auto M, class T, class R, class... Args>
        requires std::is_base_of_v<T, Class>
    [[gnu::always_inline]] static constexpr auto
    invoker_ovl(R (T::*)(Args...)) {
        return do_invoke<M, void, Class, R, Args...>;
    }

  public:
    /// Returns a function that accepts a void pointer, casts it to the class
    /// type of the member function @p Method, launders it, and then invokes
    /// @p Method with it, passing on the arguments to @p Method. The function
    /// can also accept additional arguments at the end, of type @p ExtraArgs.
    template <auto Method>
    [[gnu::always_inline]] static constexpr auto invoker() {
        return invoker_ovl<Method>(Method);
    }
};
} // namespace detail

/// @copydoc detail::Launderer::invoker
template <class Class, auto Method, class... ExtraArgs>
[[gnu::always_inline]] constexpr auto type_erased_wrapped() {
    return detail::Launderer<Class, ExtraArgs...>::template invoker<Method>();
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

/// Similar to `std::in_place_t`.
template <class T>
struct te_in_place_t {};
/// Convenience instance of @ref te_in_place_t.
template <class T>
inline constexpr te_in_place_t<T> te_in_place;

/// Class for polymorphism through type erasure. Saves the entire vtable, and
/// uses small buffer optimization.
///
/// @todo   Decouple allocation/small buffer optimization.
template <class VTable           = BasicVTable,
          class Allocator        = std::allocator<std::byte>,
          size_t SmallBufferSize = default_te_buffer_size<VTable, Allocator>()>
class TypeErased {
  public:
    static constexpr size_t small_buffer_size = SmallBufferSize;
    using allocator_type                      = Allocator;

  private:
    using allocator_traits = std::allocator_traits<allocator_type>;
    using buffer_type      = std::array<std::byte, small_buffer_size>;
    [[no_unique_address]] alignas(std::max_align_t) buffer_type small_buffer;
    [[no_unique_address]] allocator_type allocator;

  private:
    /// True if @p T is not a child class of @ref TypeErased.
    template <class T>
    static constexpr auto no_child_of_ours =
        !std::is_base_of_v<TypeErased, std::remove_cvref_t<T>>;

  protected:
    static constexpr size_t invalid_size =
        static_cast<size_t>(0xDEADBEEFDEADBEEF);
    /// Pointer to the stored object.
    void *self = nullptr;
    /// Size required to store the object.
    size_t size = invalid_size;
    VTable vtable;

  public:
    /// Default constructor.
    TypeErased() noexcept(noexcept(allocator_type())) = default;
    /// Default constructor (allocator aware).
    template <class Alloc>
    TypeErased(std::allocator_arg_t, const Alloc &alloc) : allocator{alloc} {}
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
        // Check for self-assignment
        if (&other == this)
            return *this;
        // Delete our own storage before assigning a new value
        cleanup();
        do_copy_assign<true>(other);
        return *this;
    }

    /// Move constructor.
    TypeErased(TypeErased &&other) noexcept
        : allocator{std::move(other.allocator)} {
        size   = other.size;
        vtable = std::move(other.vtable);
        // If dynamically allocated, simply steal storage
        if (size > small_buffer_size) {
            // We stole the allocator, so we can steal the storage as well
            self = std::exchange(other.self, nullptr);
        }
        // Otherwise, use the small buffer and do an explicit move
        else if (other.self) {
            self = small_buffer.data();
            vtable.move(other.self, self);
            vtable.destroy(other.self); // nothing to deallocate
            other.self = nullptr;
        }
    }
    /// Move constructor (allocator aware).
    TypeErased(TypeErased &&other, const allocator_type &alloc) noexcept
        : allocator{alloc} {
        // Only continue if other actually contains a value
        if (other.self == nullptr)
            return;
        size   = other.size;
        vtable = std::move(other.vtable);
        // If dynamically allocated, simply steal other's storage
        if (size > small_buffer_size) {
            // Can we steal the storage because of equal allocators?
            if (allocator == other.allocator) {
                self = std::exchange(other.self, nullptr);
            }
            // If the allocators are not the same, we cannot steal the
            // storage, so do an explicit move
            else {
                self = allocator.allocate(size);
                vtable.move(other.self, self);
                // Cannot call other.cleanup() here because we stole the vtable
                vtable.destroy(other.self);
                other.deallocate();
            }
        }
        // Otherwise, use the small buffer and do an explicit move
        else if (other.self) {
            self = small_buffer.data();
            vtable.move(other.self, self);
            // Cannot call other.cleanup() here because we stole the vtable
            vtable.destroy(other.self); // nothing to deallocate
            other.self = nullptr;
        }
    }
    /// Move assignment.
    TypeErased &operator=(TypeErased &&other) noexcept {
        // Check for self-assignment
        if (&other == this)
            return *this;
        // Delete our own storage before assigning a new value
        cleanup();
        // Check if we are allowed to steal the allocator
        static constexpr bool prop_alloc =
            allocator_traits::propagate_on_container_move_assignment::value;
        if constexpr (prop_alloc)
            allocator = std::move(other.allocator);
        // Only assign if other contains a value
        if (other.self == nullptr)
            return *this;

        size   = other.size;
        vtable = std::move(other.vtable);
        // If dynamically allocated, simply steal other's storage
        if (size > small_buffer_size) {
            // Can we steal the storage because of equal allocators?
            if (prop_alloc || allocator == other.allocator) {
                self = std::exchange(other.self, nullptr);
            }
            // If the allocators are not the same, we cannot steal the
            // storage, so do an explicit move
            else {
                self = allocator.allocate(size);
                vtable.move(other.self, self);
                vtable.destroy(other.self);
                // Careful, we might have moved other.allocator!
                auto &deallocator    = prop_alloc ? allocator : other.allocator;
                using pointer_t      = typename allocator_traits::pointer;
                auto &&other_pointer = static_cast<pointer_t>(other.self);
                deallocator.deallocate(other_pointer, size);
                other.self = nullptr;
            }
        }
        // Otherwise, use the small buffer and do an explicit move
        else if (other.self) {
            self = small_buffer.data();
            vtable.move(other.self, self);
            vtable.destroy(other.self); // nothing to deallocate
            other.self = nullptr;
        }
        return *this;
    }

    /// Destructor.
    ~TypeErased() { cleanup(); }

    /// Main constructor that type-erases the given argument.
    template <class T, class Alloc>
    explicit TypeErased(std::allocator_arg_t, const Alloc &alloc, T &&d)
        : allocator{alloc} {
        construct_inplace<std::remove_cvref_t<T>>(std::forward<T>(d));
    }
    /// Main constructor that type-erases the object constructed from the given
    /// argument.
    template <class T, class Alloc, class... Args>
    explicit TypeErased(std::allocator_arg_t, const Alloc &alloc,
                        te_in_place_t<T>, Args &&...args)
        : allocator{alloc} {
        construct_inplace<std::remove_cvref_t<T>>(std::forward<Args>(args)...);
    }
    /// @copydoc TypeErased(std::allocator_arg_t, const Alloc &, T &&)
    /// Requirement prevents this constructor from taking precedence over the
    /// copy and move constructors.
    template <class T>
        requires no_child_of_ours<T>
    explicit TypeErased(T &&d) {
        construct_inplace<std::remove_cvref_t<T>>(std::forward<T>(d));
    }
    /// Main constructor that type-erases the object constructed from the given
    /// argument.
    template <class T, class... Args>
    explicit TypeErased(te_in_place_t<T>, Args &&...args) {
        construct_inplace<std::remove_cvref_t<T>>(std::forward<Args>(args)...);
    }

    /// Construct a type-erased wrapper of type Ret for an object of type T,
    /// initialized in-place with the given arguments.
    template <class Ret, class T, class Alloc, class... Args>
        requires std::is_base_of_v<TypeErased, Ret>
    static Ret make(std::allocator_arg_t tag, const Alloc &alloc,
                    Args &&...args) {
        Ret r{tag, alloc};
        r.template construct_inplace<T>(std::forward<Args>(args)...);
        return r;
    }
    /// Construct a type-erased wrapper of type Ret for an object of type T,
    /// initialized in-place with the given arguments.
    template <class Ret, class T, class... Args>
        requires no_leading_allocator<Args...>
    static Ret make(Args &&...args) {
        return make<Ret, T>(std::allocator_arg, allocator_type{},
                            std::forward<Args>(args)...);
    }

    /// Check if this wrapper wraps an object. False for default-constructed
    /// objects.
    explicit operator bool() const noexcept { return self != nullptr; }

    /// Get a copy of the allocator.
    allocator_type get_allocator() const noexcept { return allocator; }

    /// Query the contained type.
    [[nodiscard]] const std::type_info &type() const noexcept {
        return *vtable.type;
    }

    /// Convert the type-erased object to the given type. The type is checked
    /// in debug builds only, use with caution.
    template <class T>
    T &as() & {
        if (typeid(T) != type())
            throw bad_type_erased_type(type(), typeid(T));
        return *reinterpret_cast<T *>(self);
    }
    /// @copydoc as()
    template <class T>
    const T &as() const & {
        if (typeid(T) != type())
            throw bad_type_erased_type(type(), typeid(T));
        return *reinterpret_cast<const T *>(self);
    }
    /// @copydoc as()
    template <class T>
    const T &&as() && {
        if (typeid(T) != type())
            throw bad_type_erased_type(type(), typeid(T));
        return std::move(*reinterpret_cast<T *>(self));
    }

  private:
    /// Deallocates the storage when destroyed.
    struct Deallocator {
        TypeErased *instance;
        Deallocator(TypeErased *instance) noexcept : instance{instance} {}
        Deallocator(const Deallocator &)            = delete;
        Deallocator &operator=(const Deallocator &) = delete;
        Deallocator(Deallocator &&o) noexcept
            : instance{std::exchange(o.instance, nullptr)} {}
        Deallocator &operator=(Deallocator &&) noexcept = delete;
        void release() noexcept { instance = nullptr; }
        ~Deallocator() { instance ? instance->deallocate() : void(); }
    };

    /// Ensure that storage is available, either by using the small buffer if
    /// it is large enough, or by calling the allocator.
    /// Returns a RAII wrapper that deallocates the storage unless released.
    Deallocator allocate(size_t size) {
        assert(!self);
        assert(size != invalid_size);
        self       = size <= small_buffer_size ? small_buffer.data()
                                               : allocator.allocate(size);
        this->size = size;
        return {this};
    }

    /// Deallocate the memory without invoking the destructor.
    void deallocate() {
        assert(size != invalid_size);
        using pointer_t = typename allocator_traits::pointer;
        if (size > small_buffer_size)
            allocator.deallocate(reinterpret_cast<pointer_t>(self), size);
        self = nullptr;
    }

    /// Destroy the type-erased object (if not empty), and deallocate the memory
    /// if necessary.
    void cleanup() {
        if (self) {
            vtable.destroy(self);
            deallocate();
        }
    }

    template <bool CopyAllocator>
    void do_copy_assign(const TypeErased &other) {
        constexpr bool prop_alloc =
            allocator_traits::propagate_on_container_copy_assignment::value;
        if constexpr (CopyAllocator && prop_alloc)
            allocator = other.allocator;
        if (!other)
            return;
        vtable             = other.vtable;
        auto storage_guard = allocate(other.size);
        // If copy constructor throws, storage should be deallocated and self
        // set to null, otherwise the TypeErased destructor will attempt to call
        // the contained object's destructor, which is undefined behavior if
        // construction failed.
        vtable.copy(other.self, self);
        storage_guard.release();
    }

  protected:
    /// Ensure storage and construct the type-erased object of type T in-place.
    template <class T, class... Args>
    void construct_inplace(Args &&...args) {
        static_assert(std::is_same_v<T, std::remove_cvref_t<T>>);
        // Allocate memory
        auto storage_guard = allocate(sizeof(T));
        // Construct the stored object
        using destroyer = std::unique_ptr<T, noop_delete<T>>;
        destroyer object_guard{new (self) T{std::forward<Args>(args)...}};
        vtable = VTable{std::in_place, static_cast<T &>(*object_guard.get())};
        object_guard.release();
        storage_guard.release();
    }

    /// Call the vtable function @p f with the given arguments @p args,
    /// implicitly passing the @ref self pointer and @ref vtable reference if
    /// necessary.
    template <class Ret, class... FArgs, class... Args>
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(const void *, FArgs...),
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
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(void *, FArgs...),
                                               Args &&...args) {
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
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(const void *)) const {
        assert(f);
        assert(self);
        return f(self);
    }
    /// @copydoc call
    template <class Ret>
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(void *)) {
        assert(f);
        assert(self);
        return f(self);
    }
    /// @copydoc call
    template <class Ret>
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(const void *,
                                                        const VTable &)) const {
        assert(f);
        assert(self);
        return f(self, vtable);
    }
    /// @copydoc call
    template <class Ret>
    [[gnu::always_inline]] decltype(auto) call(Ret (*f)(void *,
                                                        const VTable &)) {
        assert(f);
        assert(self);
        return f(self, vtable);
    }
};

} // namespace alpaqa::util