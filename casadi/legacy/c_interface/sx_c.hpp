#ifndef SX_C_HPP
#define SX_C_HPP

#include "sx_c.h"
#include "../sx/sx.hpp"
#include "c_interface.hpp"

namespace CasADi{

/// Convert a pointer to an SX to a reference to an SX
SX& sx_ref(sx_ptr ptr);

/// Convert a pointer to an vector<SX> to a reference to an vector<SX>
std::vector<SX>& sx_vec(sx_vec_ptr v);

}

#endif // SX_C_HPP
