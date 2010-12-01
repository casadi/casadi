#include "integrator_c.hpp"

using namespace CasADi;
using namespace std;

// Get a reference
Integrator& get_integrator(fx_ref ref){
  if(ref==0) throw CasadiException("get_integrator failed: null pointer");
  Integrator *r = (Integrator*)(ref);
  return *r;
}

extern "C"
int casadi_integrator_integrate(fx_ref ref, double t_out){
  try{
    get_integrator(ref).integrate(t_out);
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }
}

extern "C"
int casadi_integrator_reset(fx_ref ref, int with_sens){
  try{
    get_integrator(ref).reset(bool(with_sens));
    return 0;
  } catch(exception &e){
    cerr << e.what();
    return 1;
  } catch(...){
    return 1;
  }  
}

