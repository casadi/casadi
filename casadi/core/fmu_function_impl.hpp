/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#ifndef CASADI_FMU_FUNCTION_IMPL_HPP
#define CASADI_FMU_FUNCTION_IMPL_HPP

#include "fmu_function.hpp"
#include "function_internal.hpp"

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 protected:
   // Global unique identifier
   const std::string guid_;

   // Path to the FMU resource directory
   std::string resource_loc_;

   // Value reference to the inputs and outputs
   std::vector<std::vector<casadi_int>> id_in_, id_out_;

 public:

  /** \brief Constructor */
  FmuFunction(const std::string& name, const std::string& guid,
      const std::string& resource_loc,
      const std::vector<std::vector<casadi_int>>& id_in,
      const std::vector<std::vector<casadi_int>>& id_out);

  /** \brief Destructor */
  ~FmuFunction() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunction";}

  /// Initialize
  void init(const Dict& opts) override;

  ///@{
  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return 2;}
  size_t get_n_out() override {return 2;}
  ///@}

  ///@{
  /** \brief Names of function input and outputs */
  std::string get_name_in(casadi_int i) override { return i == 0 ? "xd" : "xn";}
  std::string get_name_out(casadi_int i) override { return i == 0 ? "yd" : "yn";}
  /// @}

  /// @{
  /** \brief Retreive differentiability */
  bool get_diff_in(casadi_int i) override {return i == 0;}
  bool get_diff_out(casadi_int i) override {return i == 0;}
  /// @}

};

} // namespace casadi
/// \endcond

#endif // CASADI_FMU_FUNCTION_IMPL_HPP
