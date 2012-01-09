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

#include "x_function.hpp"

namespace CasADi{

/** \brief  An elemenent of the algorithm, namely a binary operation */
struct SXAlgEl{
  /// operator
  unsigned char op; 
  /// index of the result
  int ind; 
  /// indices of the arguments
  int ch[2]; 
};

} // namespace CasADi

#ifdef SWIG
%template(SXAlgElVector) std::vector<CasADi::SXAlgEl>;
#endif // SWIG

namespace CasADi{

/// Forward declaration of internal class
class SXFunctionInternal;

/// Forward declaration of MXFunction
class MXFunction;

/**   \brief Dynamically created function that can be expanded into a series of scalar operations.
\author Joel Andersson 
\date 2010
*/

class SXFunction : public XFunction{

public:
  /// Default constructor
  SXFunction();
  
  /// Expand an MXFunction
  explicit SXFunction(const MXFunction &f);

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

  /// Single (scalar/vector/matrix valued) input, multiple (vector valued) output 
  SXFunction(const SXMatrix& arg, const std::vector< std::vector<SX> >& res);

  /// Single (scalar/vector/matrix valued) input, multiple (matrix valued) output 
  SXFunction(const SXMatrix& arg, const std::vector< SXMatrix>& res);
#endif // SWIG

  /// Access functions of the node 
  SXFunctionInternal* operator->();

  /// Const access functions of the node 
  const SXFunctionInternal* operator->() const;
    
  /** \brief Calculate the jacobian of output oind with respect to input iind 
  *
  * This is just the result of CasADi::SXFunction::jac,
  * wrapped in an SXFunction.
  *
  * \see CasADi::Jacobian for an AD approach
  */
  SXFunction jacobian(int iind=0, int oind=0);
  
  /// Hessian of output oind with respect to input iind 
  SXFunction hessian(int iind=0, int oind=0);

  /** \brief Jacobian via source code transformation
  *
  * \see CasADi::Jacobian for an AD approach
  */
  SXMatrix jac(int iind=0, int oind=0, bool compact=false);

  /// Gradient via source code transformation
  SXMatrix grad(int iind=0, int oind=0);
  
  /// Hessian (forward over adjoint) via source code transformation
  SXMatrix hess(int iind=0, int oind=0);
  
  /** \brief Calculate the expression for the jacobian of a number of function outputs with respect to a number of function inputs, optionally include the function outputs */
  std::vector<Matrix<SX> > jac(const std::vector<std::pair<int,int> >& jblocks, bool compact=false);
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /** \brief get function input */
  const SXMatrix& inputSX(int ind=0) const;
  
  /** \brief get function output */
  const SXMatrix& outputSX(int ind=0) const;
  
  /** \brief get function inputs */
  const std::vector<SXMatrix>& inputsSX() const;
  
  /** \brief get function outputs */
  const std::vector<SXMatrix> & outputsSX() const;
  
  /** \brief Generate C code for the function */
  void generateCode(const std::string& filename);
  
  /** \brief Access the algorithm */
  const std::vector<SXAlgEl>& algorithm() const;
  
  /** \brief Number of nodes in the algorithm */
  int countNodes() const;
  
  /** \brief Clear the function from its symbolic representation, to free up memory, no symbolic evaluations are possible after this */
  void clearSymbolic();
  
#ifndef SWIG 
  /// Construct a function that has only the k'th output
  SXFunction operator[](int k) const;
#endif //SWIG 

SXFunction indexed_one_based(int k) const{ return operator[](k-1);}
SXFunction indexed_zero_based(int k) const{ return operator[](k);}

};

} // namespace CasADi

#endif // SX_FUNCTION_HPP
