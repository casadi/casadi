#pragma once

#include <alpaqa/util/quadmath/quadmath.hpp>

#include <limits>
#include <type_traits>

#include <Eigen/Core>

namespace alpaqa {

template <class T>
struct is_config : std::false_type {};
template <class T>
inline constexpr bool is_config_v = is_config<T>::value;

template <class Conf>
concept Config = is_config_v<Conf>;

struct DefaultConfig;
template <>
struct is_config<DefaultConfig> : std::true_type {};

#define USING_ALPAQA_CONFIG_NO_TYPENAME(Conf)                                  \
    using real_t [[maybe_unused]]     = Conf::real_t;                          \
    using vec [[maybe_unused]]        = Conf::vec;                             \
    using mvec [[maybe_unused]]       = Conf::mvec;                            \
    using cmvec [[maybe_unused]]      = Conf::cmvec;                           \
    using rvec [[maybe_unused]]       = Conf::rvec;                            \
    using crvec [[maybe_unused]]      = Conf::crvec;                           \
    using mat [[maybe_unused]]        = Conf::mat;                             \
    using mmat [[maybe_unused]]       = Conf::mmat;                            \
    using cmmat [[maybe_unused]]      = Conf::cmmat;                           \
    using rmat [[maybe_unused]]       = Conf::rmat;                            \
    using crmat [[maybe_unused]]      = Conf::crmat;                           \
    using length_t [[maybe_unused]]   = Conf::length_t;                        \
    using index_t [[maybe_unused]]    = Conf::index_t;                         \
    using indexvec [[maybe_unused]]   = Conf::indexvec;                        \
    using rindexvec [[maybe_unused]]  = Conf::rindexvec;                       \
    using crindexvec [[maybe_unused]] = Conf::crindexvec

#define USING_ALPAQA_CONFIG(Conf) /** @cond CONFIG_TYPES */                    \
    using config_t [[maybe_unused]] = Conf;                                    \
    USING_ALPAQA_CONFIG_NO_TYPENAME(typename Conf) /** @endcond */

#define USING_ALPAQA_CONFIG_TEMPLATE(Conf) /** @cond CONFIG_TYPES */           \
    using config_t [[maybe_unused]] = typename Conf;                           \
    USING_ALPAQA_CONFIG_NO_TYPENAME(typename Conf) /** @endcond */

// clang-format off
template <Config Conf = DefaultConfig> using real_t = typename Conf::real_t;
template <Config Conf = DefaultConfig> using vec = typename Conf::vec;
template <Config Conf = DefaultConfig> using mvec = typename Conf::mvec;
template <Config Conf = DefaultConfig> using cmvec = typename Conf::cmvec;
template <Config Conf = DefaultConfig> using rvec = typename Conf::rvec;
template <Config Conf = DefaultConfig> using crvec = typename Conf::crvec;
template <Config Conf = DefaultConfig> using mat = typename Conf::mat;
template <Config Conf = DefaultConfig> using mmat = typename Conf::mmat;
template <Config Conf = DefaultConfig> using cmmat = typename Conf::cmmat;
template <Config Conf = DefaultConfig> using rmat = typename Conf::rmat;
template <Config Conf = DefaultConfig> using crmat = typename Conf::crmat;
template <Config Conf = DefaultConfig> using length_t = typename Conf::length_t;
template <Config Conf = DefaultConfig> using index_t = typename Conf::index_t;
template <Config Conf = DefaultConfig> using indexvec = typename Conf::indexvec;
template <Config Conf = DefaultConfig> using rindexvec = typename Conf::rindexvec;
template <Config Conf = DefaultConfig> using crindexvec = typename Conf::crindexvec;

template <Config Conf>
constexpr const auto inf = std::numeric_limits<real_t<Conf>>::infinity();
template <Config Conf>
constexpr const auto NaN = std::numeric_limits<real_t<Conf>>::quiet_NaN();
// clang-format on

template <class RealT>
struct EigenConfig {
    /// Real scalar element type.
    using real_t = RealT;
    /// Dynamic vector type.
    using vec = Eigen::VectorX<real_t>;
    /// Map of vector type.
    using mvec = Eigen::Map<vec>;
    /// Immutable map of vector type.
    using cmvec = Eigen::Map<const vec>;
    /// Reference to mutable vector.
    using rvec = Eigen::Ref<vec>;
    /// Reference to immutable vector.
    using crvec = Eigen::Ref<const vec>;
    /// Dynamic matrix type.
    using mat = Eigen::MatrixX<real_t>;
    /// Map of matrix type.
    using mmat = Eigen::Map<mat>;
    /// Immutable map of matrix type.
    using cmmat = Eigen::Map<const mat>;
    /// Reference to mutable matrix.
    using rmat = Eigen::Ref<mat>;
    /// Reference to immutable matrix.
    using crmat = Eigen::Ref<const mat>;
    /// Type for lengths and sizes.
    using length_t = Eigen::Index;
    /// Type for vector and matrix indices.
    using index_t = Eigen::Index;
    /// Dynamic vector of indices.
    using indexvec = Eigen::VectorX<index_t>;
    /// Reference to mutable index vector.
    using rindexvec = Eigen::Ref<indexvec>;
    /// Reference to immutable index vector.
    using crindexvec = Eigen::Ref<const indexvec>;
};

/// Single-precision `float` configuration.
struct EigenConfigf : EigenConfig<float> {
    static constexpr const char *get_name() { return "EigenConfigf"; }
};
/// Double-precision `double` configuration.
struct EigenConfigd : EigenConfig<double> {
    static constexpr const char *get_name() { return "EigenConfigd"; }
};
/// `long double` configuration. (Quad precision on ARM64, 80-bit x87 floats
/// on Intel/AMD x86)
struct EigenConfigl : EigenConfig<long double> {
    static constexpr const char *get_name() { return "EigenConfigl"; }
};
#ifdef ALPAQA_WITH_QUAD_PRECISION
/// Quad-precision `__float128` configuration.
struct EigenConfigq : EigenConfig<__float128> {
    static constexpr const char *get_name() { return "EigenConfigq"; }
};
template <>
struct is_config<EigenConfigq> : std::true_type {};
#endif

struct DefaultConfig : EigenConfigd {};

template <>
struct is_config<EigenConfigf> : std::true_type {};
template <>
struct is_config<EigenConfigd> : std::true_type {};
template <>
struct is_config<EigenConfigl> : std::true_type {};

/// Global empty vector for convenience.
template <Config Conf>
inline const rvec<Conf> null_vec = mvec<Conf>{nullptr, 0};

namespace vec_util {

/// Get the maximum or infinity-norm of the given vector.
/// @returns @f$ \left\|v\right\|_\infty @f$
template <class Derived>
    requires(Derived::ColsAtCompileTime == 1)
auto norm_inf(const Eigen::MatrixBase<Derived> &v) {
    return v.template lpNorm<Eigen::Infinity>();
}

/// Get the 1-norm of the given vector.
/// @returns @f$ \left\|v\right\|_1 @f$
template <class Derived>
    requires(Derived::ColsAtCompileTime == 1)
auto norm_1(const Eigen::MatrixBase<Derived> &v) {
    return v.template lpNorm<1>();
}

} // namespace vec_util

} // namespace alpaqa
