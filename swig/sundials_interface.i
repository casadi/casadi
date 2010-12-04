%{
#include "sundials_interface/cvodes_integrator.hpp"
#include "sundials_interface/idas_integrator.hpp"
%}

namespace CasADi{
namespace Sundials{

/// Input arguments of a jacobian function: J = df/dy + cj*df/dydot
enum JACInput{JAC_T, JAC_Y, JAC_YDOT, JAC_P, JAC_CJ, JAC_NUM_IN};

/// Output arguments of an DAE residual function 
enum JACOutput{JAC_J, JAC_NUM_OUT};

class CVodesIntegrator : public Integrator{
  public:
    /// Default (empty) constructor
    CVodesIntegrator();

    /// Create an integrator for explicit ODEs 
    explicit CVodesIntegrator(const FX& f, const FX& q=FX());
};

class IdasIntegrator : public Integrator{
public:
    /// Default (empty) constructor
    IdasIntegrator();
    
    /// Create an integrator for fully implicit DAEs
    explicit IdasIntegrator(const FX& f, const FX& q=FX(), const FX& jacx=FX(), const FX& jacp=FX());
};

} // namespace Sundials
} // namespace CasADi


