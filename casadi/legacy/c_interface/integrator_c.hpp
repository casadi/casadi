#ifndef INTEGRATOR_C_HPP
#define INTEGRATOR_C_HPP

#include "integrator_c.h"
#include "fx_c.hpp"
#include "../fx/integrator.hpp"

CasADi::Integrator& get_integrator(fx_ref ref);

#endif // INTEGRATOR_C_HPP 
