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

#include "nlp_solver.hpp"
#include "nlp_solver_internal.hpp"
#include "../sx/sx_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../fx/fx_tools.hpp"

using namespace std;
namespace CasADi{

  NLPSolver::NLPSolver(){
  }

  NLPSolverInternal* NLPSolver::operator->(){
    return static_cast<NLPSolverInternal*>(FX::operator->());
  }

  const NLPSolverInternal* NLPSolver::operator->() const{
    return static_cast<const NLPSolverInternal*>(FX::operator->());
  }
    
  bool NLPSolver::checkNode() const{
    return dynamic_cast<const NLPSolverInternal*>(get())!=0;
  }

  void NLPSolver::reportConstraints(std::ostream &stream) { 
    (*this)->reportConstraints();
  }

  void NLPSolver::setQPOptions() {
    (*this)->setQPOptions();
  }

  FX NLPSolver::nlp(){ 
    return (*this)->nlp_;
  }

  FX NLPSolver::gradF(){ 
    return (*this)->gradF();
  }
  
  FX NLPSolver::jacG(){ 
    return (*this)->jacG();
  }

  FX NLPSolver::hessLag(){  
    return (*this)->hessLag();
  }

  FX NLPSolver::joinFG(FX F, FX G){    
    if(G.isNull()){
      // unconstrained
      if(is_a<SXFunction>(F)){
        SXFunction F_sx = shared_cast<SXFunction>(F);
        vector<SXMatrix> nlp_in = F_sx.inputExpr();
        nlp_in.resize(NLP_NUM_IN);
        vector<SXMatrix> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_F] = F_sx.outputExpr(0);        
        return SXFunction(nlp_in,nlp_out);
      } else if(is_a<MXFunction>(F)){
        MXFunction F_mx = shared_cast<MXFunction>(F);
        vector<MX> nlp_in = F_mx.inputExpr();
        nlp_in.resize(NLP_NUM_IN);
        vector<MX> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_F] = F_mx.outputExpr(0);        
        return MXFunction(nlp_in,nlp_out);
      } else {
        vector<MX> F_in = F.symbolicInput();
        vector<MX> nlp_in(NLP_NUM_IN);
        nlp_in[NLP_X] = F_in.at(0);
        if(F_in.size()>1) nlp_in[NLP_P] = F_in.at(1);        
        vector<MX> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_F] = F.call(F_in).front();
        return MXFunction(nlp_in,nlp_out);
      }
    } else if(F.isNull()){
      // feasibility problem
      if(is_a<SXFunction>(G)){
        SXFunction G_sx = shared_cast<SXFunction>(G);
        vector<SXMatrix> nlp_in = G_sx.inputExpr();
        nlp_in.resize(NLP_NUM_IN);
        vector<SXMatrix> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_G] = G_sx.outputExpr(0);        
        return SXFunction(nlp_in,nlp_out);
      } else if(is_a<MXFunction>(G)){
        MXFunction G_mx = shared_cast<MXFunction>(F);
        vector<MX> nlp_in = G_mx.inputExpr();
        nlp_in.resize(NLP_NUM_IN);        
        vector<MX> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_G] = G_mx.outputExpr(0);        
        nlp_out.resize(NLP_NUM_OUT);
        return MXFunction(nlp_in,nlp_out);
      } else {
        vector<MX> G_in = G.symbolicInput();
        vector<MX> nlp_in(NLP_NUM_IN);
        nlp_in[NLP_X] = G_in.at(0);
        if(G_in.size()>1) nlp_in[NLP_P] = G_in.at(1);        
        vector<MX> nlp_out(NLP_NUM_OUT);
        nlp_out[NLP_G] = G.call(G_in).at(0);
        return MXFunction(nlp_in,nlp_out);
      }
    } else {
      // Standard (constrained) NLP
      
      // SXFunction if both functions are SXFunction
      if(is_a<SXFunction>(F) && is_a<SXFunction>(G)){
        vector<SXMatrix> nlp_in(NLP_NUM_IN), nlp_out(NLP_NUM_OUT);
        SXFunction F_sx = shared_cast<SXFunction>(F);
        SXFunction G_sx = shared_cast<SXFunction>(G);        
        nlp_in[NLP_X] = G_sx.inputExpr(0);
        if(G_sx.getNumInputs()>1){
          nlp_in[NLP_P] = G_sx.inputExpr(1);
        } else {
          nlp_in[NLP_P] = ssym("p",0,1);
        }

        // Expression for f and g
        nlp_out[NLP_G] = G_sx.outputExpr(0);
        nlp_out[NLP_F] = substitute(F_sx.outputExpr(),F_sx.inputExpr(),G_sx.inputExpr()).front();

        return SXFunction(nlp_in,nlp_out);
      } else { // MXFunction otherwise
        vector<MX> nlp_in(NLP_NUM_IN), nlp_out(NLP_NUM_OUT);

        // Try to cast into MXFunction
        MXFunction F_mx = shared_cast<MXFunction>(F);
        MXFunction G_mx = shared_cast<MXFunction>(G);

        // Convert to MX if cast failed and make sure that they use the same expressions if cast was successful
        if(!G_mx.isNull()){
          nlp_in[NLP_X] = G_mx.inputExpr(0);
          if(G_mx.getNumInputs()>1){
            nlp_in[NLP_P] = G_mx.inputExpr(1);
          } else {
            nlp_in[NLP_P] = msym("p",0,1);
          }
          nlp_out[NLP_G] = G_mx.outputExpr(0);
          if(!F_mx.isNull()){ // Both are MXFunction, make sure they use the same variables
            nlp_out[NLP_F] = substitute(F_mx.outputExpr(),F_mx.inputExpr(),G_mx.inputExpr()).front();
          } else { // G_ but not F_ MXFunction
            nlp_out[NLP_F] = F.call(G_mx.inputExpr()).front();
          }
        } else {
          if(!F_mx.isNull()){ // F but not G MXFunction
            nlp_in[NLP_X] = F_mx.inputExpr(0);
            if(F_mx.getNumInputs()>1){
              nlp_in[NLP_P] = F_mx.inputExpr(1);
            } else {
              nlp_in[NLP_P] = msym("p",0,1);              
            }
            nlp_out[NLP_F] = F_mx.outputExpr(0);
            nlp_out[NLP_G] = G.call(F_mx.inputExpr()).front();
          } else { // None of them MXFunction
            vector<MX> FG_in = G.symbolicInput();
            nlp_in[NLP_X] = FG_in.at(0);
            if(FG_in.size()>1) nlp_in[NLP_P] = FG_in.at(1);
            nlp_out[NLP_G] = G.call(FG_in).front();
            nlp_out[NLP_F] = F.call(FG_in).front();
          }
        }
        return MXFunction(nlp_in,nlp_out);
      } // SXFunction/MXFunction
    } // constrained/unconstrained
  }


} // namespace CasADi

