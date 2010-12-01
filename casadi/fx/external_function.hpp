/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef EXTERNAL_FUNCTION_HPP
#define EXTERNAL_FUNCTION_HPP

#include "c_function.hpp"
#include <string>

namespace CasADi{

  
/** \brief  Forward declaration of internal class */
class ExternalFunctionNode;

/** \brief  Interface for a function that is not implemented by CasADi symbolics
  \author Joel Andersson 
  \date 2010
	*/
class ExternalFunction : public CFunction{

public:

/** \brief  CONSTRUCTORS: */
/** \brief  default constructor */
  ExternalFunction();

/** \brief  Create an empty function */
  explicit ExternalFunction(const std::string& bin_name);

/** \brief  Access functions of the node */
  ExternalFunctionNode* operator->();
  const ExternalFunctionNode* operator->() const;
    
};
  
  
#if 0  
class ExternalFunctionNode : public FXNode{
  friend class ExternalFunction;
  public:
/** \brief  no public constructors */
  ~ExternalFunctionNode();

/** \brief  Set the input ind of derivative order ord */
  virtual void setArgument(const double *x, int ind=0, int ord=0);

/** \brief  Get the input */
  virtual void getArgument(double *x, int ind=0, int ord=0) const;

/** \brief  Set the output */
  virtual void setResult(const double *x, int ind=0, int ord=0);

/** \brief  Get the output */
  virtual void getResult(double *x, int ind=0, int ord=0) const;

/** \brief  initialize */
  virtual void init();

/** \brief  Clear the memory */
  virtual void clear(int ord=0);

/** \brief  Evaluate */
  virtual void evaluate(int tape_order=0);

/** \brief  Evaluate forward derivatives */
  virtual void evaluateFwd(bool use_tape=false);

/** \brief  Evaluate adjoint derivatives */
  virtual void evaluateAdj();

  protected:
/** \brief  constructor */
  explicit ExternalFunctionNode(const std::string& bin_name);

//@{
/** \brief  Function pointer types */
  typedef int (*evaluaterPtr)();
  typedef int (*clearerPtr)(int);
  typedef int (*setterPtr)(const double*, int, int);
  typedef int (*getterPtr)(double*, int, int);
//@}
  
//@{
/** \brief  Function pointers */
  setterPtr setArgument_ptr;
  getterPtr getArgument_ptr;
  setterPtr setResult_ptr;
  getterPtr getResult_ptr;
  clearerPtr clear_ptr;
  evaluaterPtr evaluate_ptr;
  evaluaterPtr evaluateFwd_ptr;
  evaluaterPtr evaluateAdj_ptr;
//@}
  
/** \brief  handle to the dll */
  void* handle;
};


#endif

} // namespace CasADi


#endif // EXTERNAL_FUNCTION_HPP
