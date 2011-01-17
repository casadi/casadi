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

#ifndef MX_FUNCTION_INTERNAL_HPP
#define MX_FUNCTION_INTERNAL_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "../mx/mx.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/** \brief  Internal node class for MX
  \author Joel Andersson 
  \date 2010
A regular user should never work with any Node class. Use MX directly.
*/
class MXFunctionInternal : public FXInternal{
  friend class MXFunction;

  public:
/** \brief  no public constructors */

/** \brief  Make a deep copy */
  virtual MXFunctionInternal* clone() const;

/** \brief  Destructor */
  ~MXFunctionInternal();

/** \brief  Order all nodes of a matrix syntax tree in the order of calculation */
  static void makeAlgorithm(MXNode* root, std::vector<MXNode*> &nodes, std::map<const MXNode*,int>  &nodemap);

/** \brief  Find a runtime element corresponding to a matrix expression */
  int findEl(const MX& mx) const;

/** \brief  Set the value of an node in the algorithm */
//  virtual void setElement(int ind, const double* x, int ord=0);

/** \brief  Get the value of an node in the algorithm */
//  virtual void getElement(int ind, double *x, int ord=0) const;

/** \brief  Evaluate the algorithm */
  virtual void evaluate(int fsens_order, int asens_order);

  /** \brief  Print representation */
  virtual void repr(std::ostream &stream) const;

  /** \brief  Print description */
  virtual void print(std::ostream &stream) const;

/** \brief  Initialize */
  virtual void init();
  
  
  protected:
/** \brief  Multiple input, multiple output constructor, only to be accessed from MXFunction, therefore protected */
  MXFunctionInternal(const std::vector<MX>& input, const std::vector<MX>& output);

/** \brief  All the runtime elements in the order of evaluation */
  std::vector<MX> alg;
  
/** \brief  Maps for quickly finding the place in the algorithm of a runtime element */
  std::map<const MXNode*,int>  nodemap;

/** \brief  Matrix expressions that are to be evaluated */
  std::vector<MX> outputv;

/** \brief  Dependent expressions */
  std::vector<MX> inputv;

  
/** \brief  Does an element exist in the algorithm */  
  bool hasEl(const MX& mx) const;

  
};

} // namespace CasADi


#endif // MX_FUNCTION_INTERNAL_HPP

