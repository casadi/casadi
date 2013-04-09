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

#include "liftopt_internal.hpp"
#include "casadi_lifter.hpp"
#include "../../symbolic/stl_vector_tools.hpp"
#include "../../symbolic/matrix/matrix_tools.hpp"
using namespace std;

namespace CasADi{

LiftoptInternal::LiftoptInternal(const MXFunction& fcn) : fcn_(fcn){
  casadi_warning("LiftoptInternal: the LIFTOPT interface is still experimental, more tests are needed");
  problem_ = 0;
  opt_ = 0;
  qp_solver_ = 0;
  
  addOption("optimizer", OT_STRING, "sqp");
  addOption("lifted", OT_BOOLEAN, true);
}

LiftoptInternal::~LiftoptInternal(){
  if(opt_!=0)     delete opt_;
  if(problem_!=0) delete problem_;
}

void LiftoptInternal::init(){
  // Initialize the NLP function
  fcn_.init();
  
  m_nCtrls_ = fcn_.input(LO_U).size();
  m_nObjRes_ = fcn_.output(LO_OBJRES).size();
  m_nEq_     = fcn_.output(LO_EQ).size();
  m_nIneq_   = fcn_.output(LO_INEQ).size();

  // Size of the NLP
  n_ = m_nCtrls_;
  m_ = m_nEq_ + m_nIneq_;

  // Call the init function for the base class
  NLPSolverInternal::init();
  
  uInit_ = liftopt::DVec(getPtr(input(NLP_SOLVER_X0)),n_,1);
  loCtrlBounds_ = liftopt::DVec(getPtr(input(NLP_SOLVER_LBX)),n_,1);
  upCtrlBounds_ = liftopt::DVec(getPtr(input(NLP_SOLVER_UBX)),n_,1);
  lambdaInit_ = liftopt::DVec(getPtr(input(NLP_SOLVER_LAM_G0)),m_,1);
  nodeInit_ = liftopt::DVec(getPtr(nodeInit),nodeInit.size(),1);
  problem_ = new CasadiLifter(this);
  
  if(getOption("optimizer")=="newton"){
    opt_ = liftopt::IOptimizer::create( liftopt::IOptimizer::Impl_Newton );
  } else if(getOption("optimizer")=="gaussnewton"){
    opt_ = liftopt::IOptimizer::create( liftopt::IOptimizer::Impl_GaussNewton );
  } else if(getOption("optimizer")=="gaussnewton_lr"){
    opt_ = liftopt::IOptimizer::create( liftopt::IOptimizer::Impl_GaussNewton_LR );
  } else if(getOption("optimizer")=="sqp"){
    opt_ = liftopt::IOptimizer::create( liftopt::IOptimizer::Impl_SQP );
  } else {
    throw CasadiException("Unknown optimizer: " + getOption("optimizer").toString());
  }

  // Allocate a QP solver
  qp_solver_ = liftopt::IQPSolver::create(liftopt::IQPSolver::Impl_qpOASES);
  opt_->setQPSolver(qp_solver_);
  opt_->setProblem(problem_);

  // Should the NLP be lifted?
  if(getOption("lifted").toInt()){
    opt_->setLiftMode( liftopt::LIFTED );
  } else {
    opt_->setLiftMode( liftopt::UNLIFTED );
  }

  // Pass the lifting function to the MXFunction
  fcn_.setLiftingFunction(CasadiLifter::liftfun, problem_);
}

void LiftoptInternal::evaluate(int nfdir, int nadir){
  // Evaluate
  opt_->setStartControls ( uInit_ );
  opt_->setStartEqMult( lambdaInit_ );
  opt_->setLowerCtrlBounds ( loCtrlBounds_ );
  opt_->setUpperCtrlBounds ( upCtrlBounds_ );
  opt_->setStartNodeValues( nodeInit_ );
  opt_->solveProblem( );
}

} // namespace CasADi

