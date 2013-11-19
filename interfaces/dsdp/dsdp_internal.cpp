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

#include "dsdp_internal.hpp"

#include "../../symbolic/stl_vector_tools.hpp"
#include "../../symbolic/matrix/matrix_tools.hpp"
#include "../../symbolic/mx/mx_tools.hpp"
#include "../../symbolic/fx/mx_function.hpp"
/**
Some implementation details

"Multiple cones can be created for the same solver, but it is usually more efficient to group all blocks into the same conic structure." user manual

*/
using namespace std;
namespace CasADi {

DSDPInternal* DSDPInternal::clone() const{
  // Return a deep copy
  DSDPInternal* node = new DSDPInternal(st_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
DSDPInternal::DSDPInternal(const std::vector<CRSSparsity> &st) : SDPSolverInternal(st){
 
  casadi_assert_message(double(m_)*(double(m_)+1)/2 < std::numeric_limits<int>::max(),"Your problem size m is too large to be handled by DSDP.");

  addOption("gapTol",OT_REAL,1e-8,"Convergence criterion based on distance between primal and dual objective");
  addOption("maxIter",OT_INTEGER,500,"Maximum number of iterations");
  addOption("dualTol",OT_REAL,1e-4,"Tolerance for dual infeasibility (translates to primal infeasibility in dsdp terms)");
  addOption("primalTol",OT_REAL,1e-4,"Tolerance for primal infeasibility (translates to dual infeasibility in dsdp terms)");
  addOption("stepTol",OT_REAL,5e-2,"Terminate the solver if the step length in the primal is below this tolerance. ");
  addOption("infinity",OT_REAL,1e30,"Treat numbers higher than this as infinity");
  addOption("_use_penalty", OT_BOOLEAN, true, "Modifies the algorithm to use a penality gamma on r.");
  addOption("_penalty", OT_REAL, 1e5, "Penality parameter lambda. Must exceed the trace of Y. This parameter heavily influences the ability of DSDP to treat linear equalities. The DSDP standard default (1e8) will make a problem with linear equality return unusable solutions.");
  addOption("_rho", OT_REAL,4.0,"Potential parameter. Must be >=1");
  addOption("_zbar", OT_REAL,1e10,"Initial upper bound on the objective of the dual problem.");
  addOption("_reuse",OT_INTEGER,4,"Maximum on the number of times the Schur complement matrix is reused");
  addOption("_printlevel",OT_INTEGER,1,"A printlevel of zero will disable all output. Another number indicates how often a line is printed.");
  addOption("_loglevel",OT_INTEGER,0,"An integer that specifies how much logging is done on stdout.");
  
  // Set DSDP memory blocks to null
  dsdp_ = 0;
  sdpcone_ = 0;
}

DSDPInternal::~DSDPInternal(){ 
  if(dsdp_!=0){
    DSDPDestroy(dsdp_);
    dsdp_ = 0;
  }
}

void DSDPInternal::init(){
  // Initialize the base classes
  SDPSolverInternal::init();
  log("DSDPInternal::init","Enter");

  terminationReason_[DSDP_CONVERGED]="DSDP_CONVERGED";
  terminationReason_[DSDP_MAX_IT]="DSDP_MAX_IT";
  terminationReason_[DSDP_INFEASIBLE_START]="DSDP_INFEASIBLE_START";
  terminationReason_[DSDP_INDEFINITE_SCHUR_MATRIX]="DSDP_INDEFINITE SCHUR";
  terminationReason_[DSDP_SMALL_STEPS]="DSDP_SMALL_STEPS";
  terminationReason_[DSDP_NUMERICAL_ERROR]="DSDP_NUMERICAL_ERROR";
  terminationReason_[DSDP_UPPERBOUND]="DSDP_UPPERBOUND";
  terminationReason_[DSDP_USER_TERMINATION]="DSDP_USER_TERMINATION";
  terminationReason_[CONTINUE_ITERATING]="CONTINUE_ITERATING";
  
  solutionType_[DSDP_PDFEASIBLE] = "DSDP_PDFEASIBLE";
  solutionType_[DSDP_UNBOUNDED] = "DSDP_UNBOUNDED";
  solutionType_[DSDP_INFEASIBLE] = "DSDP_INFEASIBLE";
  solutionType_[DSDP_PDUNKNOWN] = "DSDP_PDUNKNOWN";

  // Fill the data structures that hold DSDP-style sparse symmetric matrix
  pattern_.resize(n_+1);
  values_.resize(n_+1);
  
  for (int i=0;i<n_+1;++i) {
    pattern_[i].resize(nb_);
    values_[i].resize(nb_);
    for (int j=0;j<nb_;++j) {
      CRSSparsity CAij = mapping_.output(i*nb_+j).sparsity();
      pattern_[i][j].resize(CAij.sizeL());
      values_[i][j].resize(pattern_[i][j].size());
      int nz=0;
      vector<int> rowind,col;
      CAij.getSparsityCRS(rowind,col);
      for(int r=0; r<rowind.size()-1; ++r) {
        for(int el=rowind[r]; el<rowind[r+1]; ++el){
         if(r>=col[el]){
           pattern_[i][j][nz++] = r*(r + 1)/2 + col[el];
         }
        }
      }
      mapping_.output(i*nb_+j).get(values_[i][j],SPARSESYM);
    }
  }
  
  if (nc_>0) {
    // Fill in the linear program structure
    MX A = msym("A",input(SDP_SOLVER_A).sparsity());
    MX LBA = msym("LBA",input(SDP_SOLVER_LBA).sparsity());
    MX UBA = msym("UBA",input(SDP_SOLVER_UBA).sparsity());
    
    std::vector< MX >  syms;
    syms.push_back(A);
    syms.push_back(LBA);
    syms.push_back(UBA);
    
    // DSDP has no infinities -- replace by a big number
    // There is a way to deal with this properly, but requires modifying the dsdp source
    // We already did this for the variable bounds
    MX lba = trans(horzcat(fmax(fmin(LBA,getOption("infinity")),-(double)getOption("infinity")),A));
    MX uba = trans(horzcat(fmax(fmin(UBA,getOption("infinity")),-(double)getOption("infinity")),A));
    
    mappingA_ = MXFunction(syms,horzcat(-lba,uba));
    mappingA_.init();
  
  }
  
  if (calc_dual_) {
    store_X_.resize(nb_);
    for (int j=0;j<nb_;++j) {
      store_X_[j].resize(block_sizes_[j]*(block_sizes_[j]+1)/2);
    }
  }
  if (calc_p_) {
    store_P_.resize(nb_);
    for (int j=0;j<nb_;++j) {
      store_P_[j].resize(block_sizes_[j]*(block_sizes_[j]+1)/2);
    }
  }
  

}

void DSDPInternal::evaluate() {

  if (inputs_check_) checkInputs();
  if (print_problem_) printProblem();
  
  int info;
  
  // Seems unavoidable to do DSDPCreate here
  
  // Destroy existing DSDP instance if already allocated
  if(dsdp_!=0){
    DSDPDestroy(dsdp_);
    dsdp_ = 0;
  }
  
  // Allocate DSDP solver memory
  info = DSDPCreate(n_, &dsdp_);
  DSDPSetGapTolerance(dsdp_, getOption("gapTol"));
  DSDPSetMaxIts(dsdp_, getOption("maxIter"));
  DSDPSetPTolerance(dsdp_,getOption("dualTol"));
  DSDPSetRTolerance(dsdp_,getOption("primalTol"));
  DSDPSetStepTolerance(dsdp_,getOption("stepTol"));
  DSDPSetStandardMonitor(dsdp_,getOption("_printlevel"));
  
  info = DSDPCreateSDPCone(dsdp_,nb_,&sdpcone_);
  for (int j=0;j<nb_;++j) {
    log("DSDPInternal::init","Setting");
    info = SDPConeSetBlockSize(sdpcone_, j, block_sizes_[j]);
    info = SDPConeSetSparsity(sdpcone_, j, block_sizes_[j]);
  }
  if (nc_>0) {
    info = DSDPCreateLPCone( dsdp_, &lpcone_);
    info = LPConeSetData(lpcone_, nc_*2, &mappingA_.output(0).rowind()[0], &mappingA_.output(0).col()[0], &mappingA_.output(0).data()[0]);
  }
  
  info = DSDPCreateBCone( dsdp_, &bcone_);
  info = BConeAllocateBounds( bcone_, n_);
  
  info = DSDPUsePenalty(dsdp_,getOption("_use_penalty") ? 1: 0);
  info = DSDPSetPenaltyParameter(dsdp_, getOption("_penalty") );
  info = DSDPSetPotentialParameter(dsdp_, getOption("_rho"));
  info = DSDPSetZBar(dsdp_, getOption("_zbar"));
  info = DSDPReuseMatrix(dsdp_,getOption("_reuse") ? 1: 0);
  DSDPLogInfoAllow(getOption("_loglevel"),0);
  
  // Copy bounds
  for (int i=0;i<n_;++i) {
    if(input(SDP_SOLVER_LBX).at(i)==-std::numeric_limits< double >::infinity()) {
      info = BConeSetUnboundedLower(bcone_,i+1);   
    } else {
      info = BConeSetLowerBound(bcone_,i+1,input(SDP_SOLVER_LBX).at(i));
    }
    if(input(SDP_SOLVER_UBX).at(i)==std::numeric_limits< double >::infinity()) {
      info = BConeSetUnboundedUpper(bcone_,i+1);   
    } else {
      info = BConeSetUpperBound(bcone_,i+1,input(SDP_SOLVER_UBX).at(i));
    }
  }
  
  if (nc_>0) {
    // Copy linear constraints
    mappingA_.setInput(input(SDP_SOLVER_A),0);
    mappingA_.setInput(input(SDP_SOLVER_LBA),1);
    mappingA_.setInput(input(SDP_SOLVER_UBA),2);
    mappingA_.evaluate();
    
    // TODO: this can be made non-allocating bu hacking into DSDP source code
    info = LPConeSetDataC(lpcone_, nc_*2, &mappingA_.output(0).data()[0]);

  }

  // Copy b vector
  for (int i=0;i<n_;++i) {
    info = DSDPSetDualObjective(dsdp_, i+1, -input(SDP_SOLVER_C).at(i));
  }

  
  if (nb_>0) {
    for (int i=0;i<n_+1;++i) {
      for (int j=0;j<nb_;++j) {
        info = SDPConeSetASparseVecMat(sdpcone_, j, i, block_sizes_[j], 1, 0, &pattern_[i][j][0], &values_[i][j][0], pattern_[i][j].size() );
      }
    }
  }
  
  
  if (nb_>0) {
    // Get Ai from supplied A
    mapping_.setInput(input(SDP_SOLVER_G),0);
    mapping_.setInput(input(SDP_SOLVER_F),1);
    mapping_.evaluate();

    for (int i=0;i<n_+1;++i) {
      for (int j=0;j<nb_;++j) {
        mapping_.output(i*nb_+j).get(values_[i][j],SPARSESYM);
      }
    }
  }
  
  // Manual: Do not set data into the cones after calling this routine.
  // Manual: Should be called only once for each DSDP solver created
  info = DSDPSetup(dsdp_);
  info = DSDPSolve(dsdp_);

  casadi_assert_message(info==0,"DSDPSolver failed");
  
  
  DSDPTerminationReason reason;
  DSDPStopReason(dsdp_, &reason);
  stats_["termination_reason"] = (*terminationReason_.find(reason)).second;

  DSDPSolutionType pdfeasible;
  DSDPGetSolutionType(dsdp_,&pdfeasible);
  stats_["solution_type"] =  (*solutionType_.find(pdfeasible)).second;
  
  info = DSDPGetY(dsdp_,&output(SDP_SOLVER_X).at(0),n_);
  
  double temp;
  DSDPGetDDObjective(dsdp_, &temp);
  output(SDP_SOLVER_COST).set(-temp);
  DSDPGetPPObjective(dsdp_, &temp);
  output(SDP_SOLVER_DUAL_COST).set(-temp);
  
  if (calc_dual_) {
    for (int j=0;j<nb_;++j) {
      info = SDPConeComputeX(sdpcone_, j, block_sizes_[j], &store_X_[j][0], store_X_[j].size());
      Pmapper_.input(j).set(store_X_[j],SPARSESYM);
    }
    Pmapper_.evaluate();
    std::copy(Pmapper_.output().data().begin(),Pmapper_.output().data().end(),output(SDP_SOLVER_DUAL).data().begin());
  }
  
  if (calc_p_) {
    for (int j=0;j<nb_;++j) {
      info = SDPConeComputeS(sdpcone_, j, 1.0,  &output(SDP_SOLVER_X).at(0), n_, 0, block_sizes_[j] , &store_P_[j][0], store_P_[j].size());
      Pmapper_.input(j).set(store_P_[j],SPARSESYM);
    }
    Pmapper_.evaluate();
    std::copy(Pmapper_.output().data().begin(),Pmapper_.output().data().end(),output(SDP_SOLVER_P).data().begin());
  }
      
  DSDPComputeX(dsdp_); 
  
  info = BConeCopyXSingle( bcone_, &output(SDP_SOLVER_LAM_X).at(0), n_);
  
  if (nc_>0) {
    int dummy;
    double *lam;
   
    info = LPConeGetXArray(lpcone_, &lam, &dummy);
    std::transform (lam + nc_, lam + 2*nc_, lam, &output(SDP_SOLVER_LAM_A).at(0), std::minus<double>());
  }
  
}

} // namespace CasADi
