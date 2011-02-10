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

#include "nlp_solver_internal.hpp"

using namespace std;
namespace CasADi{

NLPSolverInternal::NLPSolverInternal(const FX& F, const FX& G, const FX& H, const FX& J, const FX& GF) : F_(F), G_(G), H_(H), J_(J), GF_(GF){
  // set default options
  setOption("name",            "unnamed NLP solver"); // name of the function
    
  n_ = F_.input(0).numel();
  m_ = G_.isNull() ? 0 : G_.output(0).numel();

  input_.resize(6);
  input(NLP_X_INIT)      = DMatrix(n_,1,0);
  input(NLP_LBX)         = DMatrix(n_,1,0);
  input(NLP_UBX)         = DMatrix(n_,1,0);
  input(NLP_LBG)         = DMatrix(m_,1,0);
  input(NLP_UBG)         = DMatrix(m_,1,0);
  input(NLP_LAMBDA_INIT) = DMatrix(m_,1,0);
  
  // Allocate space for outputs
  output_.resize(5);
  output(NLP_X_OPT)      = DMatrix(n_,1,0);
  output(NLP_COST)       = DMatrix(1,1,0);
  output(NLP_LAMBDA_OPT) = DMatrix(m_,1,0);
  output(NLP_LAMBDA_LBX) = DMatrix(n_,1,0);
  output(NLP_LAMBDA_UBX) = DMatrix(n_,1,0);

  // Create a Jacobian if it does not already exists
  if(!G_.isNull() && J_.isNull()){
    J_ = G_.jacobian();
  }
}


NLPSolverInternal::~NLPSolverInternal(){
  
}

void NLPSolverInternal::init(){
  // Call the initialization method of the base class
  FXInternal::init();

  // Initialize the functions
  F_.init();
  if(!G_.isNull()) G_.init();
  if(!J_.isNull()) J_.init();
  if(!H_.isNull()) H_.init();
  if(!GF_.isNull()) GF_.init();
}


} // namespace CasADi
