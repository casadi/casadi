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


#ifndef CASADI_ONNX_FUNCTION_HPP
#define CASADI_ONNX_FUNCTION_HPP

#include "function.hpp"

namespace casadi {

  /** \defgroup main_onnx Title
      \par

      Black-box numerical evaluation of an ONNX model.

      A model is evaluated by a runtime backend (e.g. ONNX Runtime) selected by name; the
      graph internals are opaque to CasADi. Construct one through the GraphBuilder facade:

      \verbatim
      Function f = GraphBuilder("model.onnx").create("net");
      std::vector<DM> res = f(DM(...));
      \endverbatim

      These free functions query and manage the runtime backend plugin registry.

      */

  ///@{
  /** \brief Check if a given ONNX runtime backend is available

      */
  CASADI_EXPORT bool has_onnx(const std::string& solver);

  /** \brief Load an ONNX runtime backend

      */
  CASADI_EXPORT void load_onnx(const std::string& solver);

  /** \brief List available ONNX runtime backends

      */
  CASADI_EXPORT std::vector<std::string> onnx_solvers();

  /** \brief Get documentation for an ONNX runtime backend

      */
  CASADI_EXPORT std::string onnx_doc(const std::string& solver);
  ///@}

} // namespace casadi

#endif // CASADI_ONNX_FUNCTION_HPP
