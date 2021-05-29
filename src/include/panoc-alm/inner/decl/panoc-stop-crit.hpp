#pragma once

#include <ostream>
#include <stdexcept>

namespace pa {

enum class PANOCStopCrit {
    ApproxKKT = 0,     ///< Find a ε-approximate KKT point
    ProjGradNorm,     ///< ∞-norm of the projected gradient with step size γ
    ProjGradUnitNorm, ///< ∞-norm of the projected gradient with unit step size
    FPRNorm,          ///< ∞-norm of fixed point residual
};

inline const char *enum_name(PANOCStopCrit s) {
    switch (s) {
        case PANOCStopCrit::ApproxKKT: return "ApproxKKT";
        case PANOCStopCrit::ProjGradNorm: return "ProjGradNorm";
        case PANOCStopCrit::ProjGradUnitNorm: return "ProjGradUnitNorm";
        case PANOCStopCrit::FPRNorm: return "FPRNorm";
    }
    throw std::out_of_range("invalid value for pa::PANOCStopCrit");
}

inline std::ostream &operator<<(std::ostream &os, PANOCStopCrit s) {
    return os << enum_name(s);
}

} // namespace pa