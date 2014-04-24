/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SYMBOLIC_NLP_HPP
#define SYMBOLIC_NLP_HPP

#include <casadi/core/casadi.hpp>

#include <casadi/nonlinear_programming/casadi_nonlinear_programming_export.h>


namespace casadi {

// Forward declaration
class SymbolicNLPInternal;

/** \brief A symbolic NLP representation
  \date 2012
  \author Joel Andersson
*/
class CASADI_NONLINEAR_PROGRAMMING_EXPORT SymbolicNLP : public PrintableObject {
  public:

    /** @name Symbolic representation of the NLP
    *  Data members
    */
    ///@{

      /// Variables
      SX x;

      /// Objective functions
      SX f;

      /// Constraint functions
      SX g;

      /// Bounds on x
      DMatrix x_lb, x_ub;

      /// Bounds on g
      DMatrix g_lb, g_ub;

      /// Primal initial guess
      DMatrix x_init;

      /// Dual initial guess
      DMatrix lambda_init;

    ///@}

    /// Parse an AMPL och PyOmo NL-file
    void parseNL(const std::string& filename, const Dictionary& options = Dictionary());

#ifndef SWIG
    /// Print a description of the object
    virtual void print(std::ostream &stream=std::cout) const;

    /// Print a representation of the object
    virtual void repr(std::ostream &stream=std::cout) const;

  protected:

    /// Read an expression from an NL-file (Polish infix format)
    static SXElement readExpressionNL(std::istream &stream, const std::vector<SXElement>& v);

#endif // SWIG




};

} // namespace casadi

#endif //SYMBOLIC_NLP_HPP
