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

#ifndef MX_FUNCTION_HPP
#define MX_FUNCTION_HPP

#include <set>
#include <iostream>

#include "../mx/mx.hpp"
#include "x_function.hpp"

namespace CasADi{

/** \brief  An elemenent of the algorithm, namely an MX node */
struct MXAlgEl{
  // Function to be evaluated
  MX mx;
  
  // Numerical value of the node
  FunctionIO val;

  // Pointers to be passed to evaluate
  VDptr input;
  Dptr output;
  VVDptr fwdSeed;
  VDptr  fwdSens;
  VDptr adjSeed;
  VVDptr adjSens;
  
  // Indices of the children nodes
  std::vector<int> ch;
};

} // namespace CasADi

#ifdef SWIG
// Template instantiation
%template(MXAlgElVector) std::vector<CasADi::MXAlgEl>;
#endif // SWIG

// Lifting function to be passed to the evalutator in order to lift the evaluations
typedef void (*LiftingFunction)(double *v, int n, void *user_data);

namespace CasADi{

/** \brief  Forward declaration of internal class */
class MXFunctionInternal;

  /** \brief  General function mapping from/to MX
  \author Joel Andersson 
  \date 2010
*/
class MXFunction : public XFunction{
public:

  /** \brief  Ddefault constructor */
  MXFunction();

#ifndef SWIG  
  /** \brief  Single input, single output */
  MXFunction(const MX& input, const MX& output);

  /** \brief  Single input, multiple output */
  MXFunction(const MX& input, const std::vector<MX>& output);

  /** \brief  Multiple input, single output */
  MXFunction(const std::vector<MX>& input, const MX& output);
#endif // SWIG  

  /** \brief  Multiple input, multiple output*/
  MXFunction(const std::vector<MX>& input, const std::vector<MX>& output);

  /** \brief  Access functions of the node */
  MXFunctionInternal* operator->();

  /** \brief  Const access functions of the node */
  const MXFunctionInternal* operator->() const;

  /// get function input argument 
  const MX inputMX(int iind=0) const;
  
  /// get function output argument
  const MX outputMX(int oind=0) const;
  
  /** \brief Access the algorithm */
  const std::vector<MXAlgEl>& algorithm() const;
  
  /** \brief Set the lifting function */
  void setLiftingFunction(LiftingFunction liftfun, void* user_data=0);
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /** \brief Jacobian via source code transformation (work in progress! */
  std::vector<MX> jac(int iind=0);
  
};

#ifdef SWIG
%extend MXFunction{
  // Print (why is this not inherited?)
  std::string __repr__()  { return $self->getRepresentation(); }
}
#endif


  
} // namespace CasADi


#endif // MX_FUNCTION_HPP

