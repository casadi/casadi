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
    
    /// Create an integrator which integrates the ODE/DAE augmented with the forward sensitivity equations
    CVodesIntegrator jac(int iind=0, int oind=0);
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

    /// Create an integrator which integrates the ODE/DAE augmented with the forward sensitivity equations
    IdasIntegrator jac(int iind=0, int oind=0);

};

} // namespace Sundials
} // namespace CasADi


