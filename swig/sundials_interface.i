%{
#include "interfaces/sundials/cvodes_integrator.hpp"
#include "interfaces/sundials/idas_integrator.hpp"
%}

namespace CasADi{
namespace Sundials{

/// Input arguments of the Jacobian in the nonlinear iteration: M = 1 - gamma*df/dy
enum MInput{M_T, M_Y, M_P, M_GAMMA, M_NUM_IN};

/// Output arguments of the Jacobian function
enum MOutput{M_M, M_NUM_OUT};

class CVodesIntegrator : public Integrator{
  public:
    /// Default (empty) constructor
    CVodesIntegrator();

    /// Create an integrator for explicit ODEs 
    explicit CVodesIntegrator(const FX& f, const FX& q=FX());
};

/// Input arguments of a jacobian function: J = df/dy + cj*df/dydot
enum JACInput{JAC_T, JAC_Y, JAC_YDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/// Output arguments of an DAE residual function 
enum JACOutput{JAC_J, JAC_NUM_OUT};

class IdasIntegrator : public Integrator{
public:
    /// Default (empty) constructor
    IdasIntegrator();
    
    /// Create an integrator for fully implicit DAEs
    explicit IdasIntegrator(const FX& f, const FX& q=FX());
};

} // namespace Sundials
} // namespace CasADi


