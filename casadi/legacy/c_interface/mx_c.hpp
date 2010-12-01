#ifndef MX_C_HPP
#define MX_C_HPP

#include "mx_c.h"
#include "../mx/mx.hpp"

CasADi::MX& get_mx(mx_ref ref);
std::vector<CasADi::MX>& get_mx_vec(mx_vec v);
std::string& get_stl_string(string_ptr str);

#endif // MX_C_HPP
