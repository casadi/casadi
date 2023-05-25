#include <alpaqa/inner/directions/panoc/anderson.hpp>
#include <alpaqa/inner/directions/panoc/lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-lbfgs.hpp>
#include <alpaqa/inner/directions/panoc/structured-newton.hpp>
#include <alpaqa/inner/directions/pantr/newton-tr.hpp>
#include <alpaqa/inner/internal/lipschitz.hpp>
#include <alpaqa/inner/internal/panoc-stop-crit.hpp>
#include <alpaqa/inner/panoc.hpp>
#include <alpaqa/inner/pantr.hpp>
#include <alpaqa/inner/zerofpr.hpp>
#include <alpaqa/outer/alm.hpp>
#include <alpaqa/util/io/csv.hpp>
#if ALPAQA_WITH_OCP
#include <alpaqa/inner/panoc-ocp.hpp>
#endif

#include <alpaqa/implementation/params/params.tpp>

#include <fstream>

#include "from_chars-compat.ipp"

namespace alpaqa::params {

template <>
void ALPAQA_EXPORT set_param(bool &b, ParamString s) {
    assert_key_empty<bool>(s);
    if (s.value == "0" || s.value == "false")
        b = false;
    else if (s.value == "1" || s.value == "true")
        b = true;
    else
        throw std::invalid_argument(
            "Invalid value '" + std::string(s.value) +
            "' for type 'bool' in '" + std::string(s.full_key) +
            "',\n  "
            "possible values are: '0', '1', 'true', 'false'");
}

template <>
void ALPAQA_EXPORT set_param(std::string_view &v, ParamString s) {
    assert_key_empty<bool>(s);
    v = s.value;
}

template <>
void ALPAQA_EXPORT set_param(std::string &v, ParamString s) {
    assert_key_empty<bool>(s);
    v = s.value;
}

template <class T>
    requires((std::floating_point<T> || std::integral<T>) && !std::is_enum_v<T>)
void set_param(T &f, ParamString s) {
    assert_key_empty<T>(s);
    const auto *val_end = s.value.data() + s.value.size();
    const auto *ptr     = set_param_float_int(f, s);
    if (ptr != val_end)
        throw std::invalid_argument("Invalid suffix '" +
                                    std::string(ptr, val_end) + "' for type '" +
                                    demangled_typename(typeid(T)) + "' in '" +
                                    std::string(s.full_key) + "'");
}

template <>
void ALPAQA_EXPORT set_param(alpaqa::vec<config_t> &v, ParamString s) {
    v.resize(std::count(s.value.begin(), s.value.end(), ',') + 1);
    std::string_view value, remainder = s.value;
    for (auto &e : v) {
        std::tie(value, remainder) = split_key(remainder, ',');
        set_param(e, {.full_key = s.full_key, .key = "", .value = value});
    }
}

template <>
void ALPAQA_EXPORT set_param(vec_from_file<config_t> &v, ParamString s) {
    assert_key_empty<vec_from_file<config_t>>(s);
    if (s.value.starts_with('@')) {
        std::string fpath{s.value.substr(1)};
        std::ifstream f(fpath);
        if (!f)
            throw std::invalid_argument("Unable to open file '" + fpath +
                                        "' in '" + std::string(s.full_key) +
                                        '\'');
        try {
            auto r      = alpaqa::csv::read_row_std_vector<real_t<config_t>>(f);
            auto r_size = static_cast<length_t<config_t>>(r.size());
            if (v.expected_size >= 0 && r_size != v.expected_size)
                throw std::invalid_argument(
                    "Incorrect size in '" + std::string(s.full_key) +
                    "' (got " + std::to_string(r.size()) + ", expected " +
                    std::to_string(v.expected_size) + ')');
            v.value.emplace(cmvec<config_t>{r.data(), r_size});
        } catch (alpaqa::csv::read_error &e) {
            throw std::invalid_argument(
                "Unable to read from file '" + fpath + "' in '" +
                std::string(s.full_key) +
                "': alpaqa::csv::read_error: " + e.what());
        }
    } else {
        alpaqa::params::set_param(v.value.emplace(), s);
        if (v.expected_size >= 0 && v.value->size() != v.expected_size)
            throw std::invalid_argument(
                "Incorrect size in '" + std::string(s.full_key) + "' (got " +
                std::to_string(v.value->size()) + ", expected " +
                std::to_string(v.expected_size) + ')');
    }
}

template <class Rep, class Period>
void set_param(std::chrono::duration<Rep, Period> &t, ParamString s) {
    using Duration = std::remove_cvref_t<decltype(t)>;
    assert_key_empty<Duration>(s);
    const auto *val_end = s.value.data() + s.value.size();
    double value;
#if ALPAQA_USE_FROM_CHARS_FLOAT
    auto [ptr, ec] = std::from_chars(s.value.data(), val_end, value);
    if (ec != std::errc())
        throw std::invalid_argument("Invalid value '" +
                                    std::string(ptr, val_end) + "' for type '" +
                                    demangled_typename(typeid(Duration)) +
                                    "' in '" + std::string(s.full_key) +
                                    "': " + std::make_error_code(ec).message());
#else
#pragma message "Using std::stod as a fallback to replace std::from_chars"
    size_t end_index;
    try {
        value = std::stod(std::string(s.value), &end_index);
    } catch (std::exception &e) {
        throw std::invalid_argument(
            "Invalid value '" + std::string(s.value) + "' for type '" +
            demangled_typename(typeid(Duration)) + "' in '" +
            std::string(s.full_key) + "': " + e.what());
    }
    const char *ptr = s.value.data() + end_index;
#endif
    std::string_view units{ptr, val_end};
    auto cast = [](auto t) { return std::chrono::duration_cast<Duration>(t); };
    if (units == "s" || units.empty())
        t = cast(std::chrono::duration<double, std::ratio<1, 1>>{value});
    else if (units == "ms")
        t = cast(std::chrono::duration<double, std::ratio<1, 1000>>{value});
    else if (units == "us" || units == "µs")
        t = cast(std::chrono::duration<double, std::ratio<1, 1000000>>{value});
    else if (units == "ns")
        t = cast(
            std::chrono::duration<double, std::ratio<1, 1000000000>>{value});
    else if (units == "min")
        t = cast(std::chrono::duration<double, std::ratio<60, 1>>{value});
    else
        throw std::invalid_argument("Invalid units '" + std::string(units) +
                                    "' in '" + std::string(s.full_key) + "'");
}

template <>
void ALPAQA_EXPORT set_param(LBFGSStepSize &t, ParamString s) {
    if (s.value == "BasedOnExternalStepSize")
        t = LBFGSStepSize::BasedOnExternalStepSize;
    else if (s.value == "BasedOnCurvature")
        t = LBFGSStepSize::BasedOnCurvature;
    else
        throw std::invalid_argument("Invalid value '" + std::string(s.value) +
                                    "' for type 'LBFGSStepSize' in '" +
                                    std::string(s.full_key) + "'");
}

template <>
void ALPAQA_EXPORT set_param(PANOCStopCrit &t, ParamString s) {
    if (s.value == "ApproxKKT")
        t = PANOCStopCrit::ApproxKKT;
    else if (s.value == "ApproxKKT2")
        t = PANOCStopCrit::ApproxKKT2;
    else if (s.value == "ProjGradNorm")
        t = PANOCStopCrit::ProjGradNorm;
    else if (s.value == "ProjGradNorm2")
        t = PANOCStopCrit::ProjGradNorm2;
    else if (s.value == "ProjGradUnitNorm")
        t = PANOCStopCrit::ProjGradUnitNorm;
    else if (s.value == "ProjGradUnitNorm2")
        t = PANOCStopCrit::ProjGradUnitNorm2;
    else if (s.value == "FPRNorm")
        t = PANOCStopCrit::FPRNorm;
    else if (s.value == "FPRNorm2")
        t = PANOCStopCrit::FPRNorm2;
    else if (s.value == "Ipopt")
        t = PANOCStopCrit::Ipopt;
    else if (s.value == "LBFGSBpp")
        t = PANOCStopCrit::LBFGSBpp;
    else
        throw std::invalid_argument("Invalid value '" + std::string(s.value) +
                                    "' for type 'PANOCStopCrit' in '" +
                                    std::string(s.full_key) + "'");
}

PARAMS_TABLE(LBFGSParams<config_t>,        //
             PARAMS_MEMBER(memory),        //
             PARAMS_MEMBER(min_div_fac),   //
             PARAMS_MEMBER(min_abs_s),     //
             PARAMS_MEMBER(cbfgs),         //
             PARAMS_MEMBER(force_pos_def), //
             PARAMS_MEMBER(stepsize),      //
);

PARAMS_TABLE(AndersonAccelParams<config_t>, //
             PARAMS_MEMBER(memory),         //
             PARAMS_MEMBER(min_div_fac),    //
);

PARAMS_TABLE(CBFGSParams<config_t>, //
             PARAMS_MEMBER(α),      //
             PARAMS_MEMBER(ϵ),      //
);

PARAMS_TABLE(LipschitzEstimateParams<config_t>, //
             PARAMS_MEMBER(L_0),                //
             PARAMS_MEMBER(δ),                  //
             PARAMS_MEMBER(ε),                  //
             PARAMS_MEMBER(Lγ_factor),          //
);

PARAMS_TABLE(PANTRParams<config_t>,                                         //
             PARAMS_MEMBER(Lipschitz),                                      //
             PARAMS_MEMBER(max_iter),                                       //
             PARAMS_MEMBER(max_time),                                       //
             PARAMS_MEMBER(L_min),                                          //
             PARAMS_MEMBER(L_max),                                          //
             PARAMS_MEMBER(stop_crit),                                      //
             PARAMS_MEMBER(max_no_progress),                                //
             PARAMS_MEMBER(print_interval),                                 //
             PARAMS_MEMBER(print_precision),                                //
             PARAMS_MEMBER(quadratic_upperbound_tolerance_factor),          //
             PARAMS_MEMBER(TR_tolerance_factor),                            //
             PARAMS_MEMBER(ratio_threshold_acceptable),                     //
             PARAMS_MEMBER(ratio_threshold_good),                           //
             PARAMS_MEMBER(radius_factor_rejected),                         //
             PARAMS_MEMBER(radius_factor_acceptable),                       //
             PARAMS_MEMBER(radius_factor_good),                             //
             PARAMS_MEMBER(initial_radius),                                 //
             PARAMS_MEMBER(min_radius),                                     //
             PARAMS_MEMBER(compute_ratio_using_new_stepsize),               //
             PARAMS_MEMBER(update_direction_on_prox_step),                  //
             PARAMS_MEMBER(recompute_last_prox_step_after_direction_reset), //
             PARAMS_MEMBER(disable_acceleration),                           //
             PARAMS_MEMBER(ratio_approx_fbe_quadratic_model),               //
);

PARAMS_TABLE(PANOCParams<config_t>,                                     //
             PARAMS_MEMBER(Lipschitz),                                  //
             PARAMS_MEMBER(max_iter),                                   //
             PARAMS_MEMBER(max_time),                                   //
             PARAMS_MEMBER(min_linesearch_coefficient),                 //
             PARAMS_MEMBER(force_linesearch),                           //
             PARAMS_MEMBER(linesearch_strictness_factor),               //
             PARAMS_MEMBER(L_min),                                      //
             PARAMS_MEMBER(L_max),                                      //
             PARAMS_MEMBER(stop_crit),                                  //
             PARAMS_MEMBER(max_no_progress),                            //
             PARAMS_MEMBER(print_interval),                             //
             PARAMS_MEMBER(print_precision),                            //
             PARAMS_MEMBER(quadratic_upperbound_tolerance_factor),      //
             PARAMS_MEMBER(linesearch_tolerance_factor),                //
             PARAMS_MEMBER(update_direction_in_candidate),              //
             PARAMS_MEMBER(recompute_last_prox_step_after_lbfgs_flush), //
);

PARAMS_TABLE(ZeroFPRParams<config_t>,                                   //
             PARAMS_MEMBER(Lipschitz),                                  //
             PARAMS_MEMBER(max_iter),                                   //
             PARAMS_MEMBER(max_time),                                   //
             PARAMS_MEMBER(min_linesearch_coefficient),                 //
             PARAMS_MEMBER(force_linesearch),                           //
             PARAMS_MEMBER(linesearch_strictness_factor),               //
             PARAMS_MEMBER(L_min),                                      //
             PARAMS_MEMBER(L_max),                                      //
             PARAMS_MEMBER(stop_crit),                                  //
             PARAMS_MEMBER(max_no_progress),                            //
             PARAMS_MEMBER(print_interval),                             //
             PARAMS_MEMBER(print_precision),                            //
             PARAMS_MEMBER(quadratic_upperbound_tolerance_factor),      //
             PARAMS_MEMBER(linesearch_tolerance_factor),                //
             PARAMS_MEMBER(update_direction_in_candidate),              //
             PARAMS_MEMBER(recompute_last_prox_step_after_lbfgs_flush), //
             PARAMS_MEMBER(update_direction_from_prox_step),            //
);

PARAMS_TABLE(LBFGSDirectionParams<config_t>,              //
             PARAMS_MEMBER(rescale_on_step_size_changes), //
);

PARAMS_TABLE(AndersonDirectionParams<config_t>,           //
             PARAMS_MEMBER(rescale_on_step_size_changes), //
);

PARAMS_TABLE(StructuredLBFGSDirectionParams<config_t>,      //
             PARAMS_MEMBER(hessian_vec),                    //
             PARAMS_MEMBER(hessian_vec_finite_differences), //
             PARAMS_MEMBER(full_augmented_hessian),         //
);

PARAMS_TABLE(NewtonTRDirectionParams<config_t>,           //
             PARAMS_MEMBER(rescale_on_step_size_changes), //
             PARAMS_MEMBER(hessian_vec_factor),           //
             PARAMS_MEMBER(finite_diff),                  //
             PARAMS_MEMBER(finite_diff_stepsize),         //
);

PARAMS_TABLE(SteihaugCGParams<config_t>,     //
             PARAMS_MEMBER(tol_scale),       //
             PARAMS_MEMBER(tol_scale_root),  //
             PARAMS_MEMBER(tol_max),         //
             PARAMS_MEMBER(max_iter_factor), //
);

PARAMS_TABLE(StructuredNewtonRegularizationParams<config_t>, //
             PARAMS_MEMBER(min_eig),                         //
             PARAMS_MEMBER(print_eig),                       //
);

PARAMS_TABLE(StructuredNewtonDirectionParams<config_t>, //
             PARAMS_MEMBER(hessian_vec),                //
);

PARAMS_TABLE(ALMParams<config_t>,                           //
             PARAMS_MEMBER(tolerance),                      //
             PARAMS_MEMBER(dual_tolerance),                 //
             PARAMS_MEMBER(penalty_update_factor),          //
             PARAMS_MEMBER(penalty_update_factor_lower),    //
             PARAMS_MEMBER(min_penalty_update_factor),      //
             PARAMS_MEMBER(initial_penalty),                //
             PARAMS_MEMBER(initial_penalty_factor),         //
             PARAMS_MEMBER(initial_penalty_lower),          //
             PARAMS_MEMBER(initial_tolerance),              //
             PARAMS_MEMBER(initial_tolerance_increase),     //
             PARAMS_MEMBER(tolerance_update_factor),        //
             PARAMS_MEMBER(ρ_increase),                     //
             PARAMS_MEMBER(ρ_max),                          //
             PARAMS_MEMBER(rel_penalty_increase_threshold), //
             PARAMS_MEMBER(max_multiplier),                 //
             PARAMS_MEMBER(max_penalty),                    //
             PARAMS_MEMBER(min_penalty),                    //
             PARAMS_MEMBER(max_iter),                       //
             PARAMS_MEMBER(max_time),                       //
             PARAMS_MEMBER(max_num_initial_retries),        //
             PARAMS_MEMBER(max_num_retries),                //
             PARAMS_MEMBER(max_total_num_retries),          //
             PARAMS_MEMBER(print_interval),                 //
             PARAMS_MEMBER(print_precision),                //
             PARAMS_MEMBER(single_penalty_factor),          //
);

#if ALPAQA_WITH_OCP
PARAMS_TABLE(PANOCOCPParams<config_t>,
             PARAMS_MEMBER(Lipschitz),                             //
             PARAMS_MEMBER(max_iter),                              //
             PARAMS_MEMBER(max_time),                              //
             PARAMS_MEMBER(min_linesearch_coefficient),            //
             PARAMS_MEMBER(linesearch_strictness_factor),          //
             PARAMS_MEMBER(L_min),                                 //
             PARAMS_MEMBER(L_max),                                 //
             PARAMS_MEMBER(L_max_inc),                             //
             PARAMS_MEMBER(stop_crit),                             //
             PARAMS_MEMBER(max_no_progress),                       //
             PARAMS_MEMBER(gn_interval),                           //
             PARAMS_MEMBER(gn_sticky),                             //
             PARAMS_MEMBER(reset_lbfgs_on_gn_step),                //
             PARAMS_MEMBER(lqr_factor_cholesky),                   //
             PARAMS_MEMBER(lbfgs_params),                          //
             PARAMS_MEMBER(print_interval),                        //
             PARAMS_MEMBER(print_precision),                       //
             PARAMS_MEMBER(quadratic_upperbound_tolerance_factor), //
             PARAMS_MEMBER(linesearch_tolerance_factor),           //
             PARAMS_MEMBER(disable_acceleration),                  //
);
#endif

namespace detail {

/// Check if @p A is equal to any of @p Bs.
template <class A, class... Bs>
constexpr bool any_is_same() {
    return (std::is_same_v<A, Bs> || ...);
}

/// Unused unique type tag for template specializations that were rejected
/// because some types were not distinct.
template <class...>
struct _dummy;

/// If @p NewAlias is not the same type as any of @p PossibleAliases, the result
/// is @p NewAlias. If @p NewAlias is not distinct from @p PossibleAliases, the
/// result is a dummy type, uniquely determined by @p NewAlias and
/// @p PossibleAliases.
template <class NewAlias, class... PossibleAliases>
using possible_alias_t =
    std::conditional_t<any_is_same<NewAlias, PossibleAliases...>(),
                       _dummy<NewAlias, PossibleAliases...>, NewAlias>;

} // namespace detail

template <class... Ts>
void set_param(detail::_dummy<Ts...> &, ParamString) {}

#define ALPAQA_SET_PARAM_INST(...)                                             \
    template void ALPAQA_EXPORT set_param(                                     \
        detail::possible_alias_t<__VA_ARGS__> &, ParamString)

ALPAQA_SET_PARAM_INST(float);
ALPAQA_SET_PARAM_INST(double, float);
ALPAQA_SET_PARAM_INST(long double, double, float);

ALPAQA_SET_PARAM_INST(int8_t);
ALPAQA_SET_PARAM_INST(uint8_t);
ALPAQA_SET_PARAM_INST(int16_t);
ALPAQA_SET_PARAM_INST(uint16_t);
ALPAQA_SET_PARAM_INST(int32_t);
ALPAQA_SET_PARAM_INST(int64_t);
ALPAQA_SET_PARAM_INST(uint32_t);
ALPAQA_SET_PARAM_INST(uint64_t);

// Here, we would like to instantiate alpaqa::params::set_param for all standard
// integer types, but the issue is that they might not be distinct types:
// For example, on some platforms, int32_t might be a weak alias to int, whereas
// on other platforms, it could be a distinct type.
// To resolve this issue, we use some metaprogramming to ensure distinct
// instantiations with unique dummy types.
#define ALPAQA_SET_PARAM_INST_INT(...)                                         \
    ALPAQA_SET_PARAM_INST(__VA_ARGS__, int8_t, uint8_t, int16_t, uint16_t,     \
                          int32_t, int64_t, uint32_t, uint64_t)

ALPAQA_SET_PARAM_INST_INT(short);
ALPAQA_SET_PARAM_INST_INT(int, short);
ALPAQA_SET_PARAM_INST_INT(long, int, short);
ALPAQA_SET_PARAM_INST_INT(long long, long, int, short);
ALPAQA_SET_PARAM_INST_INT(ptrdiff_t, long long, long, int, short);
ALPAQA_SET_PARAM_INST_INT(unsigned short);
ALPAQA_SET_PARAM_INST_INT(unsigned int, unsigned short);
ALPAQA_SET_PARAM_INST_INT(unsigned long, unsigned int, unsigned short);
ALPAQA_SET_PARAM_INST_INT(unsigned long long, unsigned long, unsigned int,
                          unsigned short);
ALPAQA_SET_PARAM_INST_INT(size_t, unsigned long long, unsigned long,
                          unsigned int, unsigned short);

ALPAQA_SET_PARAM_INST(std::chrono::nanoseconds);
ALPAQA_SET_PARAM_INST(std::chrono::microseconds);
ALPAQA_SET_PARAM_INST(std::chrono::milliseconds);
ALPAQA_SET_PARAM_INST(std::chrono::seconds);
ALPAQA_SET_PARAM_INST(std::chrono::minutes);
ALPAQA_SET_PARAM_INST(std::chrono::hours);

ALPAQA_SET_PARAM_INST(PANOCParams<config_t>);
ALPAQA_SET_PARAM_INST(ZeroFPRParams<config_t>);
ALPAQA_SET_PARAM_INST(PANTRParams<config_t>);
ALPAQA_SET_PARAM_INST(LBFGSParams<config_t>);
ALPAQA_SET_PARAM_INST(AndersonAccelParams<config_t>);
ALPAQA_SET_PARAM_INST(LBFGSDirectionParams<config_t>);
ALPAQA_SET_PARAM_INST(AndersonDirectionParams<config_t>);
ALPAQA_SET_PARAM_INST(StructuredLBFGSDirectionParams<config_t>);
ALPAQA_SET_PARAM_INST(NewtonTRDirectionParams<config_t>);
ALPAQA_SET_PARAM_INST(SteihaugCGParams<config_t>);
ALPAQA_SET_PARAM_INST(StructuredNewtonRegularizationParams<config_t>);
ALPAQA_SET_PARAM_INST(StructuredNewtonDirectionParams<config_t>);
ALPAQA_SET_PARAM_INST(ALMParams<config_t>);
#if ALPAQA_WITH_OCP
ALPAQA_SET_PARAM_INST(PANOCOCPParams<config_t>);
#endif

} // namespace alpaqa::params
