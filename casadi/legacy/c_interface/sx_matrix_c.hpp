#ifndef SX_MATRIX_C_HPP
#define SX_MATRIX_C_HPP

#include "sx_matrix_c.h"
#include "../sx/sx_matrix.hpp"
#include "c_interface.hpp"
#include "sx_c.hpp"

namespace CasADi{

  SXMatrix& get_sx_matrix(sx_matrix_ref ref);

  std::vector<SXMatrix>& get_sx_matrix_vec(sx_matrix_vec v);

}
  
#endif // SX_MATRIX_C_HPP
