%{
#include "sundials_interface/cvodes_integrator.hpp"
#include "sundials_interface/idas_integrator.hpp"
%}

namespace CasADi{
namespace Sundials{

class CVodesIntegrator : public Integrator{
  public:
    // Create an integrator for explicit ODEs 
    explicit CVodesIntegrator(const FX& f, const FX& q=FX());
};

class IdasIntegrator : public Integrator{
public:
  // Create an integrator for fully implicit DAEs
  explicit IdasIntegrator(const FX& f, const FX& q=FX());
};

} // namespace Sundials
} // namespace CasADi


