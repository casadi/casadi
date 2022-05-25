#pragma once

#ifdef ALPAQA_WITH_QUAD_PRECISION
#pragma GCC system_header

#include <quadmath.h>

#include <alpaqa/export.h>

#include <cmath>
#include <iosfwd>
#include <limits>

namespace std {

// From GCC 11's <limits>

/// numeric_limits<__float128> specialization.
template <>
struct numeric_limits<__float128> {
    static constexpr bool is_specialized = true;

    static constexpr __float128 min() noexcept { return FLT128_MIN; }
    static constexpr __float128 max() noexcept { return FLT128_MAX; }
    static constexpr __float128 lowest() noexcept { return -FLT128_MAX; }

    static constexpr int digits       = FLT128_MANT_DIG;
    static constexpr int digits10     = FLT128_DIG;
    static constexpr int max_digits10 = (2 + FLT128_MANT_DIG * 643L / 2136);
    static constexpr bool is_signed   = true;
    static constexpr bool is_integer  = false;
    static constexpr bool is_exact    = false;
    static constexpr int radix        = 2;

    static constexpr __float128 epsilon() noexcept { return FLT128_EPSILON; }
    static constexpr __float128 round_error() noexcept { return 0.5Q; }

    static constexpr int min_exponent   = FLT128_MIN_EXP;
    static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
    static constexpr int max_exponent   = FLT128_MAX_EXP;
    static constexpr int max_exponent10 = FLT128_MAX_10_EXP;

    static constexpr bool has_infinity             = true;
    static constexpr bool has_quiet_NaN            = true;
    static constexpr bool has_signaling_NaN        = has_quiet_NaN;
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss          = false;

    static constexpr __float128 infinity() noexcept { return __builtin_huge_valq(); }
    static constexpr __float128 quiet_NaN() noexcept { return __builtin_nanq(""); }
    static constexpr __float128 signaling_NaN() noexcept { return __builtin_nansq(""); }
    static constexpr __float128 denorm_min() noexcept { return FLT128_DENORM_MIN; }

    static constexpr bool is_iec559 = has_infinity && has_quiet_NaN && has_denorm == denorm_present;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo  = false;

    static constexpr bool traps                    = false;
    static constexpr bool tinyness_before          = false;
    static constexpr float_round_style round_style = round_to_nearest;
};

// From GCC's quadmath.h

inline __float128 acos(__float128 x) { return ::acosq(x); }
inline __float128 acosh(__float128 x) { return ::acoshq(x); }
inline __float128 asin(__float128 x) { return ::asinq(x); }
inline __float128 asinh(__float128 x) { return ::asinhq(x); }
inline __float128 atan(__float128 x) { return ::atanq(x); }
inline __float128 atanh(__float128 x) { return ::atanhq(x); }
inline __float128 atan2(__float128 x, __float128 y) { return ::atan2q(x, y); }
inline __float128 cbrt(__float128 x) { return ::cbrtq(x); }
inline __float128 ceil(__float128 x) { return ::ceilq(x); }
inline __float128 copysign(__float128 x, __float128 y) { return ::copysignq(x, y); }
inline __float128 cosh(__float128 x) { return ::coshq(x); }
inline __float128 cos(__float128 x) { return ::cosq(x); }
inline __float128 erf(__float128 x) { return ::erfq(x); }
inline __float128 erfc(__float128 x) { return ::erfcq(x); }
inline __float128 exp2(__float128 x) { return ::exp2q(x); }
inline __float128 exp(__float128 x) { return ::expq(x); }
inline __float128 expm1(__float128 x) { return ::expm1q(x); }
inline __float128 fabs(__float128 x) { return ::fabsq(x); }
inline __float128 fdim(__float128 x, __float128 y) { return ::fdimq(x, y); }
inline int finite(__float128 x) { return ::finiteq(x); }
inline __float128 floor(__float128 x) { return ::floorq(x); }
inline __float128 fma(__float128 x, __float128 y, __float128 z) { return ::fmaq(x, y, z); }
inline __float128 fmax(__float128 x, __float128 y) { return ::fmaxq(x, y); }
inline __float128 fmin(__float128 x, __float128 y) { return ::fminq(x, y); }
inline __float128 fmod(__float128 x, __float128 y) { return ::fmodq(x, y); }
inline __float128 frexp(__float128 x, int *y) { return ::frexpq(x, y); }
inline __float128 hypot(__float128 x, __float128 y) { return ::hypotq(x, y); }
inline int isinf(__float128 x) { return ::isinfq(x); }
inline int ilogb(__float128 x) { return ::ilogbq(x); }
inline int isfinite(__float128 x) { return ::finiteq(x); }
inline int isnan(__float128 x) { return ::isnanq(x); }
inline int issignaling(__float128 x) { return ::issignalingq(x); }
inline __float128 j0(__float128 x) { return ::j0q(x); }
inline __float128 j1(__float128 x) { return ::j1q(x); }
inline __float128 jn(int n, __float128 x) { return ::jnq(n, x); }
inline __float128 ldexp(__float128 x, int y) { return ::ldexpq(x, y); }
inline __float128 lgamma(__float128 x) { return ::lgammaq(x); }
inline long long int llrint(__float128 x) { return ::llrintq(x); }
inline long long int llround(__float128 x) { return ::llroundq(x); }
inline __float128 logb(__float128 x) { return ::logbq(x); }
inline __float128 log(__float128 x) { return ::logq(x); }
inline __float128 log10(__float128 x) { return ::log10q(x); }
inline __float128 log2(__float128 x) { return ::log2q(x); }
inline __float128 log1p(__float128 x) { return ::log1pq(x); }
inline long int lrint(__float128 x) { return ::lrintq(x); }
inline long int lround(__float128 x) { return ::lroundq(x); }
inline __float128 modf(__float128 x, __float128 *y) { return ::modfq(x, y); }
// inline __float128 nan(const char *x) { return ::nanq(x); }
inline __float128 nearbyint(__float128 x) { return ::nearbyintq(x); }
inline __float128 nextafter(__float128 x, __float128 y) { return ::nextafterq(x, y); }
inline __float128 pow(__float128 x, __float128 y) { return ::powq(x, y); }
inline __float128 remainder(__float128 x, __float128 y) { return ::remainderq(x, y); }
inline __float128 remquo(__float128 x, __float128 y, int *z) { return ::remquoq(x, y, z); }
inline __float128 rint(__float128 x) { return ::rintq(x); }
inline __float128 round(__float128 x) { return ::roundq(x); }
inline __float128 scalbln(__float128 x, long int y) { return ::scalblnq(x, y); }
inline __float128 scalbn(__float128 x, int y) { return ::scalbnq(x, y); }
inline int signbit(__float128 x) { return ::signbitq(x); }
inline void sincosq(__float128 x, __float128 *y, __float128 *z) { return ::sincosq(x, y, z); }
inline __float128 sinh(__float128 x) { return ::sinhq(x); }
inline __float128 sin(__float128 x) { return ::sinq(x); }
inline __float128 sqrt(__float128 x) { return ::sqrtq(x); }
inline __float128 tan(__float128 x) { return ::tanq(x); }
inline __float128 tanh(__float128 x) { return ::tanhq(x); }
inline __float128 tgamma(__float128 x) { return ::tgammaq(x); }
inline __float128 trunc(__float128 x) { return ::truncq(x); }
inline __float128 y0(__float128 x) { return ::y0q(x); }
inline __float128 y1(__float128 x) { return ::y1q(x); }
inline __float128 yn(int n, __float128 x) { return ::ynq(n, x); }

} // namespace std

#endif