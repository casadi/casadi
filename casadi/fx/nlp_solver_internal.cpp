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

NLPSolverInternal::NLPSolverInternal(){
  // set default options
  setOption("name",            "unnamed NLP solver"); // name of the function
  
  n_ = 0;
  m_ = 0;
  pn_ = 0;
  pm_ = 0;
}

NLPSolverInternal::~NLPSolverInternal(){
}

void NLPSolverInternal::init(){
  input_.resize(NLP_NUM_IN);
  input(NLP_X_INIT)      = DMatrix(n_,1,0);
  input(NLP_LBX)         = DMatrix(n_,1,0);
  input(NLP_UBX)         = DMatrix(n_,1,0);
  input(NLP_LBG)         = DMatrix(m_,1,0);
  input(NLP_UBG)         = DMatrix(m_,1,0);
  input(NLP_LAMBDA_INIT) = DMatrix(m_,1,0);
  input(NLP_P)           = DMatrix(pn_,pm_,0);
  
  // Allocate space for outputs
  output_.resize(NLP_NUM_OUT);
  output(NLP_X_OPT)      = DMatrix(n_,1,0);
  output(NLP_COST)       = DMatrix(1,1,0);
  output(NLP_LAMBDA_OPT) = DMatrix(m_,1,0);
  output(NLP_LAMBDA_LBX) = DMatrix(n_,1,0);
  output(NLP_LAMBDA_UBX) = DMatrix(n_,1,0);

  // Call the initialization method of the base class
  FXInternal::init();
}


} // namespace CasADi
