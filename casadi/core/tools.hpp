/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_TOOLS_HPP
#define CASADI_TOOLS_HPP

#include "casadi/core/function.hpp"

namespace casadi {

/** \brief Apply a transformation defined externally

\param name Name of the shared library
\param op Name of the operation
\param f Function to transform
\param opts Options

    \identifier{27i} */
CASADI_EXPORT Function external_transform(const std::string& name,
                    const std::string& op,
                    const Function& f,
                    const Dict& opts=Dict());


typedef void (*external_print_callback_t)(const char* s);
typedef const char* (*external_transform_t)(char api_version, const char* casadi_version,
    const char* in,
    external_print_callback_t cb_stdout, external_print_callback_t cb_stderr);

} // namespace casadi

#ifndef SWIG
extern "C" {
  CASADI_EXPORT const char* external_transform_test_success__f(char api_version,
    const char* casadi_version,
    const char* in,
    casadi::external_print_callback_t cb_stdout, casadi::external_print_callback_t cb_stderr);
  CASADI_EXPORT const char* external_transform_test_fail__f(char api_version,
    const char* casadi_version,
    const char* in,
    casadi::external_print_callback_t cb_stdout, casadi::external_print_callback_t cb_stderr);
}
#endif // SWIG

#endif // CASADI_TOOLS_HPP
