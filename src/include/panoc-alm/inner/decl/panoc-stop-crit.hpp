#pragma once

#include <ostream>
#include <stdexcept>

namespace pa {

enum class PANOCStopCrit {
    /// Find an ε-approximate KKT point in the ∞-norm:
    /// @f[
    ///     \varepsilon = \left\| \gamma_k^{-1} (x^k - \hat x^k) +
    ///     \nabla \psi(\hat x^k) - \nabla \psi(x^k) \right\|_\infty
    /// @f]
    ApproxKKT = 0,
    /// Find an ε-approximate KKT point in the 2-norm:
    /// @f[
    ///     \varepsilon = \left\| \gamma_k^{-1} (x^k - \hat x^k) +
    ///     \nabla \psi(\hat x^k) - \nabla \psi(x^k) \right\|_2
    /// @f]
    ApproxKKT2,
    /// ∞-norm of the projected gradient with step size γ:
    /// @f[
    ///     \varepsilon = \left\| x^k -
    ///     \Pi_C\left(x^k - \gamma_k \nabla \psi(x^k)\right) \right\|_\infty
    /// @f]
    ProjGradNorm,
    /// 2-norm of the projected gradient with step size γ:
    /// @f[
    ///     \varepsilon = \left\| x^k -
    ///     \Pi_C\left(x^k - \gamma_k \nabla \psi(x^k)\right) \right\|_2
    /// @f]
    ProjGradNorm2,
    /// ∞-norm of the projected gradient with unit step size:
    /// @f[
    ///     \varepsilon = \left\| x^k -
    ///     \Pi_C\left(x^k - \nabla \psi(x^k)\right) \right\|_\infty
    /// @f]
    ProjGradUnitNorm,
    /// 2-norm of the projected gradient with unit step size:
    /// @f[
    ///     \varepsilon = \left\| x^k -
    ///     \Pi_C\left(x^k - \nabla \psi(x^k)\right) \right\|_2
    /// @f]
    ProjGradUnitNorm2,
    /// ∞-norm of fixed point residual:
    /// @f[
    ///     \varepsilon = \gamma_k^{-1} \left\| x^k -
    ///     \Pi_C\left(x^k - \gamma_k \nabla \psi(x^k)\right) \right\|_\infty
    /// @f]
    FPRNorm,
    /// 2-norm of fixed point residual:
    /// @f[
    ///     \varepsilon = \gamma_k^{-1} \left\| x^k -
    ///     \Pi_C\left(x^k - \gamma_k \nabla \psi(x^k)\right) \right\|_2
    /// @f]
    FPRNorm2,
};

inline const char *enum_name(PANOCStopCrit s) {
    switch (s) {
        case PANOCStopCrit::ApproxKKT: return "ApproxKKT";
        case PANOCStopCrit::ApproxKKT2: return "ApproxKKT2";
        case PANOCStopCrit::ProjGradNorm: return "ProjGradNorm";
        case PANOCStopCrit::ProjGradNorm2: return "ProjGradNorm2";
        case PANOCStopCrit::ProjGradUnitNorm: return "ProjGradUnitNorm";
        case PANOCStopCrit::ProjGradUnitNorm2: return "ProjGradUnitNorm2";
        case PANOCStopCrit::FPRNorm: return "FPRNorm";
        case PANOCStopCrit::FPRNorm2: return "FPRNorm2";
    }
    throw std::out_of_range("invalid value for pa::PANOCStopCrit");
}

inline std::ostream &operator<<(std::ostream &os, PANOCStopCrit s) {
    return os << enum_name(s);
}

} // namespace pa