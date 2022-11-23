#include <gtest/gtest.h>

#include <alpaqa/util/type-erasure.hpp>

#include <cstdio>
#include <cstring>
#include <string_view>

#define PF() (std::printf("%s: %s\n", __PRETTY_FUNCTION__, msg))

struct Noisy {
    static unsigned created;
    static unsigned destroyed;
    Noisy(const char *msg = "default") {
        std::strcpy(this->msg, msg);
        PF();
        ++created;
    }
    Noisy(const Noisy &other) {
        std::strcpy(this->msg, other.msg);
        PF();
        ++created;
    }
    Noisy(Noisy &&other) {
        std::strcpy(this->msg, other.msg);
        std::strcpy(other.msg, "moved-from");
        PF();
        ++created;
    }
    Noisy &operator=(const Noisy &other) {
        std::strcpy(this->msg, other.msg);
        return PF(), *this;
    }
    Noisy &operator=(Noisy &&other) {
        std::strcpy(this->msg, other.msg);
        std::strcpy(other.msg, "moved-from");
        return PF(), *this;
    }
    ~Noisy() {
        PF();
        ++destroyed;
    }
    const char *get_msg() const { return PF(), msg; }
    void set_msg(const char *m) {
        PF();
        std::strcpy(msg, m);
    }
    char msg[16]{};
};

struct CustomVTable : alpaqa::util::BasicVTable {
    const char *(*get_msg)(const void *)  = nullptr;
    void (*set_msg)(void *, const char *) = nullptr;
    CustomVTable()                        = default;
    template <class T>
    CustomVTable(alpaqa::util::VTableTypeTag<T> t) : BasicVTable{t} {
        get_msg = [](const void *self) {
            return std::launder(reinterpret_cast<const T *>(self))->get_msg();
        };
        set_msg = [](void *self, const char *m) {
            std::launder(reinterpret_cast<T *>(self))->set_msg(m);
        };
    }
};

template <class Alloc = std::allocator<std::byte>>
struct CustomTypeErased : alpaqa::util::TypeErased<CustomVTable, Alloc, 0> {
    using TypeErased = alpaqa::util::TypeErased<CustomVTable, Alloc, 0>;
    using TypeErased::self;
    using TypeErased::TypeErased;
    using TypeErased::vtable;
    using typename TypeErased::allocator_type;
    template <class T, class... Args>
    static CustomTypeErased make(Args &&...args) {
        return TypeErased::template make<CustomTypeErased, T>(
            std::forward<Args>(args)...);
    }
    const char *get_msg() const { return vtable.get_msg(self); }
    void set_msg(const char *m) { vtable.set_msg(self, m); }
};

unsigned Noisy::created   = 0;
unsigned Noisy::destroyed = 0;

#include <array>
#include <cctype>
#include <iomanip>
#include <memory_resource>
#include <vector>

using PMRCTE = CustomTypeErased<std::pmr::polymorphic_allocator<std::byte>>;

void dump_buf(const auto &buffer) {
    constexpr size_t per_row = 16;
    size_t rows              = (buffer.size() + per_row - 1) / per_row;
    for (size_t r = 0; r < rows; ++r) {
        size_t rr   = r * per_row;
        size_t cols = std::min(per_row, buffer.size() - rr);
        std::cout << std::hex << std::setw(4) << rr << std::dec << ": ";
        for (size_t c = 0; c < cols; ++c) {
            auto ch = static_cast<uint8_t>(buffer[rr + c]);
            std::cout << std::hex << std::setw(2) << +ch << ' ' << std::dec;
        }
        std::cout << "  ";
        for (size_t c = 0; c < cols; ++c) {
            auto ch = static_cast<char>(buffer[rr + c]);
            if (std::isprint(ch)) {
                std::cout << "  " << ch;
            } else {
                std::cout << "  ·";
            }
        }
        std::cout << '\n';
    }
}

TEST(TypeErasure, TypeErased) {
    std::array<std::byte, 256> buffer{};
    {
        std::pmr::monotonic_buffer_resource mbr{
            buffer.data(), buffer.size(), std::pmr::null_memory_resource()};
        std::cout << "default: " << std::pmr::get_default_resource() << '\n';
        std::cout << "mbr:     " << &mbr << '\n';
        PMRCTE n{Noisy{"Message"}, &mbr};
        dump_buf(buffer);
        EXPECT_STREQ(n.get_msg(), "Message");
        EXPECT_TRUE(n);
        PMRCTE n2{std::move(n), &mbr};
        dump_buf(buffer);
        EXPECT_FALSE(n);
        EXPECT_TRUE(n2);
        EXPECT_STREQ(n2.get_msg(), "Message");
        PMRCTE n3{std::allocator_arg, &mbr};
        EXPECT_FALSE(n3);
        dump_buf(buffer);
        n3 = n2;
        EXPECT_TRUE(n2);
        EXPECT_TRUE(n3);
        dump_buf(buffer);
        EXPECT_STREQ(n2.get_msg(), "Message");
        EXPECT_STREQ(n3.get_msg(), "Message");
        n3.set_msg("foo");
        EXPECT_STREQ(n2.get_msg(), "Message");
        EXPECT_STREQ(n3.get_msg(), "foo");
        std::cout << sizeof(PMRCTE) << "\n";
        std::pmr::vector<PMRCTE> v{&mbr};
        v.reserve(2);
        dump_buf(buffer);
        v.emplace_back<Noisy>("vec 0");
        dump_buf(buffer);
        EXPECT_STREQ(v.back().get_msg(), "vec 0");
        v.emplace_back<Noisy>("vec 1");
        dump_buf(buffer);
        EXPECT_STREQ(v.front().get_msg(), "vec 0");
        EXPECT_STREQ(v.back().get_msg(), "vec 1");
    }
    std::cout << "Created:   " << Noisy::created << "\n"
              << "Destroyed: " << Noisy::destroyed << "\n";
    EXPECT_EQ(Noisy::created, Noisy::destroyed);
    dump_buf(buffer);
    auto buf_str = std::string_view{
        reinterpret_cast<const char *>(buffer.data()), buffer.size()};
    using namespace std::string_view_literals;
    EXPECT_NE(buf_str.find("Message"sv), buf_str.npos);
    EXPECT_NE(buf_str.find("foo\0age"sv), buf_str.npos);
    EXPECT_NE(buf_str.find("vec 0"sv), buf_str.npos);
    EXPECT_NE(buf_str.find("vec 1"sv), buf_str.npos);
}

TEST(TypeErasure, copyFromEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>::make<Noisy>();
    b      = a;
}

TEST(TypeErasure, copyToEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>::make<Noisy>();
    a      = b;
}

TEST(TypeErasure, copyEmptyToEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>();
    a      = b;
}

TEST(TypeErasure, copyConstructFromEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b{a};
}

TEST(TypeErasure, moveFromEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>::make<Noisy>();
    b      = std::move(a);
}

TEST(TypeErasure, moveToEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>::make<Noisy>();
    a      = std::move(b);
}

TEST(TypeErasure, moveEmptyToEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b = CustomTypeErased<>();
    a      = std::move(b);
}

TEST(TypeErasure, moveConstructFromEmpty) {
    struct test_exception {};
    auto a = CustomTypeErased<>();
    auto b{std::move(a)};
}

TEST(TypeErasure, throwingCopyCtor) {
    struct test_exception {};
    struct Throwing : Noisy {
        using Noisy::Noisy;
        Throwing(const Throwing &o) : Noisy{o} { throw test_exception(); }
    };
    auto a = CustomTypeErased<>::make<Throwing>();
    EXPECT_THROW({ auto b{a}; }, test_exception);
    auto c = CustomTypeErased<>();
    EXPECT_THROW({ c = a; }, test_exception);
}

TEST(TypeErasure, throwingCtor) {
    struct test_exception {};
    struct Throwing : Noisy {
        Throwing() { throw test_exception(); }
    };
    EXPECT_THROW(CustomTypeErased<>::make<Throwing>(), test_exception);
}
