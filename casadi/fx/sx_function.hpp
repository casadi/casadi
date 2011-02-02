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

#ifndef SX_FUNCTION_HPP
#define SX_FUNCTION_HPP

#include "fx.hpp"
#include "../sx/sx.hpp"
#include <vector>

namespace CasADi{

/** \brief  An elemenent of the algorithm, namely a binary operation */
struct SXAlgEl{
  unsigned short op; // operator
  int ind; // index of the binary operaton to be evaluated
  int ch[2]; // indices of the arguments
};
  
/// Forward declaration of internal class
class SXFunctionInternal;

/**   \brief Dynamically created function that can be expanded into a series of scalar operations.
\author Joel Andersson 
\date 2010
*/

class SXFunction : public FX{

public:
  /// Default constructor
  SXFunction();

  /// Multiple (matrix valued) input, multiple (matrix valued) output 
  SXFunction(const std::vector< SXMatrix>& arg, const std::vector<SXMatrix>& res);

#ifndef SWIG

  /// Multiple (vector valued) input, multiple (vector valued) output 
  SXFunction(const std::vector< std::vector<SX> >& arg, const std::vector< std::vector<SX> >& res);

  /// Single (scalar/matrix/vector valued) input, single (scalar/matrix/vector valued) output  
  SXFunction(const SXMatrix& arg, const SXMatrix& res);

  /// Multiple (vector valued) input, single (scalar/vector/matrix valued) output 
  SXFunction(const std::vector< std::vector<SX> >& arg, const SXMatrix& res);

  /// Multiple (matrix valued) input, single (scalar/vector/matrix valued) output 
  SXFunction(const std::vector< SXMatrix>& arg, const SXMatrix& res);
#endif // SWIG

  /// Access functions of the node 
  SXFunctionInternal* operator->();

  /// Const access functions of the node 
  const SXFunctionInternal* operator->() const;
  
  /// evaluate symbolically 
  std::vector<SXMatrix> eval(const std::vector<SXMatrix>& arg);

  /// evaluate symbolically (pass and get non-zero entries) 
  std::vector< std::vector<SX> > eval(const std::vector< std::vector<SX> >& arg);

#ifndef SWIG
  /// evaluate symbolically, single input, single output 
  SXMatrix eval(const SXMatrix& arg);

  /// evaluate symbolically, single input, single output (pass and get non-zero entries) 
  std::vector<SX> eval(const std::vector<SX>& arg);
#endif // SWIG
  
  /// Jacobian of output oind with respect to input iind 
  SXFunction jacobian(int iind=0, int oind=0);
  
  /// Hessian of output oind with respect to input iind 
  SXFunction hessian(int iind=0, int oind=0);

  /// Jacobian via source code transformation
  SXMatrix jac(int iind=0, int oind=0);

  /// Gradient via source code transformation
  SXMatrix grad(int iind=0, int oind=0);
  
  /// Hessian (forward over adjoint) via source code transformation
  SXMatrix hess(int iind=0, int oind=0);
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /** \brief get function input */
  const SXMatrix& inputSX(int ind=0) const;
  
  /** \brief get function output */
  const SXMatrix& outputSX(int ind=0) const;
  
  /** \brief Generate C code for the function */
  void generateCode(const std::string& filename) const;
  
  /** \brief Access the algorithm */
  const std::vector<SXAlgEl>& algorithm() const;
};

#ifdef SWIG
%extend SXFunction{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif

} // namespace CasADi

#endif // SX_FUNCTION_HPP
