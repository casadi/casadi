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

#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP

#include <set>
#include <map>
#include <vector>
#include <iostream>

#include "fx.hpp"

namespace CasADi{


  
  
/** \brief  Forward declaration of internal class */
class JacobianNode;

/** \brief Jacobian class	
  Universal Jacobian class, calculates the Jacobian of a function based on compression techniques 
  \author Joel Andersson 
  \date 2010
    joel.andersson@esat.kuleuven.be
*/ 
class Jacobian : public FX{
public:

/** \brief  Default constructor */
  Jacobian();

/** \brief  Copy constructor */
  Jacobian(const Jacobian& ref);

/** \brief  Create a Jacobian */
  explicit Jacobian(const FX& fcn, int iind=0, int oind=0);

/** \brief  Access functions of the node */
  JacobianNode* operator->();
  const JacobianNode* operator->() const;
};


/** \brief  Internal node class for Jacobian
  \author Joel Andersson 
  \date 2010
*/
class JacobianNode : public FXNode{
  friend class Jacobian;
  
  protected:
/** \brief  Multiple input, multiple output constructor, only to be accessed from Jacobian, therefore protected */
  JacobianNode(const FX& fcn, int iind, int oind);

  public:
  
/** \brief  Destructor */
  virtual ~JacobianNode();

/** \brief  Evaluate the algorithm */
  virtual void evaluate(int fsens_order, int asens_order);

/** \brief  Initialize */
  virtual void init();
  
  bool sparse_jac_, use_fd_, use_ad_fwd_;
  
  std::vector<int> elind_;
  std::vector<int> rr_, cc_, rind_;
  std::vector<double> epsilon_; // perturbations
  
  protected:
    
  // Dimensions
  int n_,m_;
  
/** \brief  Function to be differentiated */
  FX fcn_;
  
/** \brief  Output and input indices */
  int oind_, iind_;

};



} // namespace CasADi


#endif // JACOBIAN_HPP

