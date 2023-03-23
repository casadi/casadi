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
#if ALPAQA_WITH_OCP
#include <alpaqa/inner/panoc-ocp.hpp>
#endif

#include <alpaqa/implementation/params/params.tpp>

namespace alpaqa::params {

template <>
void set_param(bool &b, ParamString s) {
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
void set_param(std::string_view &v, ParamString s) {
    assert_key_empty<bool>(s);
    v = s.value;
}

template <>
void set_param(std::string &v, ParamString s) {
    assert_key_empty<bool>(s);
    v = s.value;
}

template <class T>
    requires((std::floating_point<T> || std::integral<T>) && !std::is_enum_v<T>)
void set_param(T &f, ParamString s) {
    assert_key_empty<T>(s);
    const auto *val_end = s.value.data() + s.value.size();
    auto [ptr, ec]      = std::from_chars(s.value.data(), val_end, f);
    if (ec != std::errc())
        throw std::invalid_argument(
            "Invalid value '" + std::string(s.value) + "' for type '" +
            demangled_typename(typeid(T)) + "' in '" + std::string(s.full_key) +
            "': " + std::make_error_code(ec).message());
    if (ptr != val_end)
        throw std::invalid_argument("Invalid suffix '" +
                                    std::string(ptr, val_end) + "' for type '" +
                                    demangled_typename(typeid(T)) + "' in '" +
                                    std::string(s.full_key) + "'");
}

template void ALPAQA_EXPORT set_param(real_t<config_t> &, ParamString);
template void ALPAQA_EXPORT set_param(index_t<config_t> &, ParamString);
template void ALPAQA_EXPORT set_param(int &, ParamString);
template void ALPAQA_EXPORT set_param(unsigned &, ParamString);

template <>
void set_param(alpaqa::vec<config_t> &v, ParamString s) {
    v.resize(std::count(s.value.begin(), s.value.end(), ',') + 1);
    std::string_view value, remainder = s.value;
    for (auto &e : v) {
        std::tie(value, remainder) = split_key(remainder, ',');
        set_param(e, {.full_key = s.full_key, .key = "", .value = value});
    }
}

template <class Rep, class Period>
void set_param(std::chrono::duration<Rep, Period> &t, ParamString s) {
    using Duration = std::remove_cvref_t<decltype(t)>;
    assert_key_empty<Duration>(s);
    const auto *val_end = s.value.data() + s.value.size();
    double value;
    auto [ptr, ec] = std::from_chars(s.value.data(), val_end, value);
    if (ec != std::errc())
        throw std::invalid_argument("Invalid value '" +
                                    std::string(ptr, val_end) + "' for type '" +
                                    demangled_typename(typeid(Duration)) +
                                    "' in '" + std::string(s.full_key) +
                                    "': " + std::make_error_code(ec).message());
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
void set_param(LBFGSStepSize &t, ParamString s) {
    using enum LBFGSStepSize;
    if (s.value == "BasedOnExternalStepSize")
        t = BasedOnExternalStepSize;
    else if (s.value == "BasedOnCurvature")
        t = BasedOnCurvature;
    else
        throw std::invalid_argument("Invalid value '" + std::string(s.value) +
                                    "' for type 'LBFGSStepSize' in '" +
                                    std::string(s.full_key) + "'");
}

template <>
void set_param(PANOCStopCrit &t, ParamString s) {
    using enum PANOCStopCrit;
    if (s.value == "ApproxKKT")
        t = ApproxKKT;
    else if (s.value == "ApproxKKT2")
        t = ApproxKKT2;
    else if (s.value == "ProjGradNorm")
        t = ProjGradNorm;
    else if (s.value == "ProjGradNorm2")
        t = ProjGradNorm2;
    else if (s.value == "ProjGradUnitNorm")
        t = ProjGradUnitNorm;
    else if (s.value == "ProjGradUnitNorm2")
        t = ProjGradUnitNorm2;
    else if (s.value == "FPRNorm")
        t = FPRNorm;
    else if (s.value == "FPRNorm2")
        t = FPRNorm2;
    else if (s.value == "Ipopt")
        t = Ipopt;
    else if (s.value == "LBFGSBpp")
        t = LBFGSBpp;
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
             PARAMS_MEMBER(μ1),                                             //
             PARAMS_MEMBER(μ2),                                             //
             PARAMS_MEMBER(c1),                                             //
             PARAMS_MEMBER(c2),                                             //
             PARAMS_MEMBER(c3),                                             //
             PARAMS_MEMBER(Δ_0),                                            //
             PARAMS_MEMBER(Δ_min),                                          //
             PARAMS_MEMBER(compute_ratio_using_new_γ),                      //
             PARAMS_MEMBER(update_direction_on_prox_step),                  //
             PARAMS_MEMBER(recompute_last_prox_step_after_direction_reset), //
             PARAMS_MEMBER(disable_acceleration),                           //
);

PARAMS_TABLE(PANOCParams<config_t>,                                     //
             PARAMS_MEMBER(Lipschitz),                                  //
             PARAMS_MEMBER(max_iter),                                   //
             PARAMS_MEMBER(max_time),                                   //
             PARAMS_MEMBER(τ_min),                                      //
             PARAMS_MEMBER(force_linesearch),                           //
             PARAMS_MEMBER(β),                                          //
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
             PARAMS_MEMBER(τ_min),                                      //
             PARAMS_MEMBER(force_linesearch),                           //
             PARAMS_MEMBER(β),                                          //
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

PARAMS_TABLE(LBFGSDirectionParams<config_t>,        //
             PARAMS_MEMBER(rescale_when_γ_changes), //
);

PARAMS_TABLE(StructuredLBFGSDirectionParams<config_t>,      //
             PARAMS_MEMBER(hessian_vec),                    //
             PARAMS_MEMBER(hessian_vec_finite_differences), //
             PARAMS_MEMBER(full_augmented_hessian),         //
);

PARAMS_TABLE(NewtonTRDirectionParams<config_t>,     //
             PARAMS_MEMBER(rescale_when_γ_changes), //
             PARAMS_MEMBER(hessian_vec_factor),     //
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

PARAMS_TABLE(ALMParams<config_t>,                    //
             PARAMS_MEMBER(ε),                       //
             PARAMS_MEMBER(δ),                       //
             PARAMS_MEMBER(Δ),                       //
             PARAMS_MEMBER(Δ_lower),                 //
             PARAMS_MEMBER(Δ_min),                   //
             PARAMS_MEMBER(Σ_0),                     //
             PARAMS_MEMBER(σ_0),                     //
             PARAMS_MEMBER(Σ_0_lower),               //
             PARAMS_MEMBER(ε_0),                     //
             PARAMS_MEMBER(ε_0_increase),            //
             PARAMS_MEMBER(ρ),                       //
             PARAMS_MEMBER(ρ_increase),              //
             PARAMS_MEMBER(ρ_max),                   //
             PARAMS_MEMBER(θ),                       //
             PARAMS_MEMBER(M),                       //
             PARAMS_MEMBER(penalty_alm_split),       //
             PARAMS_MEMBER(Σ_max),                   //
             PARAMS_MEMBER(Σ_min),                   //
             PARAMS_MEMBER(max_iter),                //
             PARAMS_MEMBER(max_time),                //
             PARAMS_MEMBER(max_num_initial_retries), //
             PARAMS_MEMBER(max_num_retries),         //
             PARAMS_MEMBER(max_total_num_retries),   //
             PARAMS_MEMBER(print_interval),          //
             PARAMS_MEMBER(print_precision),         //
             PARAMS_MEMBER(single_penalty_factor),   //
);

#if ALPAQA_WITH_OCP
PARAMS_TABLE(PANOCOCPParams<config_t>,
             PARAMS_MEMBER(Lipschitz),                             //
             PARAMS_MEMBER(max_iter),                              //
             PARAMS_MEMBER(max_time),                              //
             PARAMS_MEMBER(τ_min),                                 //
             PARAMS_MEMBER(β),                                     //
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

template void ALPAQA_EXPORT set_params(PANOCParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(ZeroFPRParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(PANTRParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(LBFGSParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(LBFGSDirectionParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT
set_params(StructuredLBFGSDirectionParams<config_t> &, std::string_view,
           std::span<const std::string_view>, std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(NewtonTRDirectionParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(SteihaugCGParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT
set_params(StructuredNewtonRegularizationParams<config_t> &, std::string_view,
           std::span<const std::string_view>, std::optional<std::span<bool>>);
template void ALPAQA_EXPORT
set_params(StructuredNewtonDirectionParams<config_t> &, std::string_view,
           std::span<const std::string_view>, std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(ALMParams<config_t> &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(length_t<config_t> &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(int &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(unsigned &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(std::string_view &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(std::string &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(bool &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
template void ALPAQA_EXPORT set_params(vec<config_t> &, std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
#if ALPAQA_WITH_OCP
template void ALPAQA_EXPORT set_params(PANOCOCPParams<config_t> &,
                                       std::string_view,
                                       std::span<const std::string_view>,
                                       std::optional<std::span<bool>>);
#endif

} // namespace alpaqa::params