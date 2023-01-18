#pragma once

#include <iosfwd>
#include <stdexcept>

namespace alpaqa {

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
    /// This is the same criterion as used by
    /// [OpEn](https://alphaville.github.io/optimization-engine/).
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
    /// The stopping criterion used by Ipopt, see
    /// https://link.springer.com/article/10.1007/s10107-004-0559-y equation (5).
    ///
    /// Given a candidate iterate @f$ \hat x^k @f$ and the corresponding
    /// candidate Lagrange multipliers @f$ \hat y^k @f$ for the general
    /// constraints @f$ g(x)\in D @f$,
    /// the multipliers @f$ w @f$ for the box constraints @f$ x\in C @f$
    /// (that are otherwise not computed explicitly) are given by
    /// @f[
    /// w^k = v^k - \Pi_C(v^k),
    /// @f]
    /// where
    /// @f[ \begin{aligned}
    /// v^k &\triangleq
    /// \hat x^k - \nabla f(\hat x^k) - \nabla g(\hat x^k)\, \hat y^k \\ &=
    /// \hat x^k - \nabla \psi(\hat x^k)
    /// \end{aligned} @f]
    /// The quantity that is compared to the (scaled) tolerance is then given by
    /// @f[ \begin{aligned}
    /// \varepsilon' &=
    /// \left\|
    ///     \nabla f(\hat x^k) + \nabla g(\hat x^k)\, \hat y^k + w^k
    /// \right\|_\infty \\ &=
    /// \left\|
    /// \hat x^k - \Pi_C\left(v^k\right)
    /// \right\|_\infty
    /// \end{aligned} @f]
    /// Finally, the quantity is scaled by the factor
    /// @f[
    /// s_d \triangleq \max\left\{
    /// s_\text{max},\;\frac{\|\hat y^k\|_1 + \|w^k\|_1}{2m + 2n}
    /// \right\} / s_\text{max},
    /// @f]
    /// i.e. @f$ \varepsilon = \varepsilon' / s_d @f$.
    Ipopt,
    /// The stopping criterion used by LBFGS++, see
    /// https://lbfgspp.statr.me/doc/classLBFGSpp_1_1LBFGSBParam.html#afb20e8fd6c6808c1f736218841ba6947
    ///
    /// @f[
    ///     \varepsilon = \frac{\left\| x^k -
    ///     \Pi_C\left(x^k - \nabla \psi(x^k)\right) \right\|_\infty}
    ///     {\max\left\{1, \|x\|_2 \right\}}
    /// @f]
    LBFGSBpp,
};

inline constexpr const char *enum_name(PANOCStopCrit s) {
    switch (s) {
        case PANOCStopCrit::ApproxKKT: return "ApproxKKT";
        case PANOCStopCrit::ApproxKKT2: return "ApproxKKT2";
        case PANOCStopCrit::ProjGradNorm: return "ProjGradNorm";
        case PANOCStopCrit::ProjGradNorm2: return "ProjGradNorm2";
        case PANOCStopCrit::ProjGradUnitNorm: return "ProjGradUnitNorm";
        case PANOCStopCrit::ProjGradUnitNorm2: return "ProjGradUnitNorm2";
        case PANOCStopCrit::FPRNorm: return "FPRNorm";
        case PANOCStopCrit::FPRNorm2: return "FPRNorm2";
        case PANOCStopCrit::Ipopt: return "Ipopt";
        case PANOCStopCrit::LBFGSBpp: return "LBFGSBpp";
        default:;
    }
    throw std::out_of_range("invalid value for alpaqa::PANOCStopCrit");
}

std::ostream &operator<<(std::ostream &os, PANOCStopCrit s);

} // namespace alpaqa