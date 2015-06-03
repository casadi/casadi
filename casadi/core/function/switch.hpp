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


#ifndef CASADI_SWITCH_HPP
#define CASADI_SWITCH_HPP

#include "function.hpp"

namespace casadi {

  /** \brief  Forward declaration of internal class */
  class SwitchInternal;

  /** Switch statement
      \author Joel Andersson
      \date 2015
  */
  class CASADI_EXPORT Switch : public Function {
  public:
    /** \brief Default constructor */
    Switch();

    /** \brief Constructor (generic switch) */
    Switch(const std::string& name, const std::vector<Function>& f,
           const Function& f_def, const Dict& opts=Dict());

    /** \brief Constructor (if-else) */
    Switch(const std::string& name, const Function& f_true,
           const Function& f_false, const Dict& opts=Dict());

    /** \brief  Access functions of the node */
    SwitchInternal* operator->();

    /** \brief  Const access functions of the node */
    const SwitchInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);
  };

} // namespace casadi

#endif // CASADI_SWITCH_HPP
