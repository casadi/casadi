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

#include "symbolic_qr_internal.hpp"
#include "../sx/sx_tools.hpp"
#include "sx_function.hpp"

#ifdef WITH_DL 
#include <cstdlib>
#endif // WITH_DL 

using namespace std;
namespace CasADi{

  SymbolicQRInternal::SymbolicQRInternal(const CRSSparsity& sparsity, int nrhs) : LinearSolverInternal(sparsity,nrhs){
    addOption("codegen",           OT_BOOLEAN,  false,               "C-code generation");
    addOption("compiler",          OT_STRING,    "gcc -fPIC -O2",    "Compiler command to be used for compiling generated code");
  }

  SymbolicQRInternal::~SymbolicQRInternal(){
  }
  
  void SymbolicQRInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    LinearSolverInternal::deepCopyMembers(already_copied);
    fact_fcn_ = deepcopy(fact_fcn_,already_copied);
    solv_fcn_N_ = deepcopy(solv_fcn_N_,already_copied);
    solv_fcn_T_ = deepcopy(solv_fcn_T_,already_copied);
  }

  void SymbolicQRInternal::init(){
    // Call the base class initializer
    LinearSolverInternal::init();

    // Read options
    bool codegen = getOption("codegen");
    string compiler = getOption("compiler");

    // Make sure that command processor is available
    if(codegen){
#ifdef WITH_DL 
      int flag = system(static_cast<const char*>(0));
      casadi_assert_message(flag!=0, "No command procesor available");
#else // WITH_DL 
      casadi_error("Codegen requires CasADi to be compiled with option \"WITH_DL\" enabled");
#endif // WITH_DL 
    }

    // Symbolic expression for A
    SXMatrix A = ssym("A",input(0).sparsity());

    // Make a BLT transformation of A
    std::vector<int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
    A.sparsity().dulmageMendelsohn(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);

    // Get the inverted row permutation
    std::vector<int> inv_rowperm(rowperm.size());
    for(int k=0; k<rowperm.size(); ++k)
      inv_rowperm[rowperm[k]] = k;

    // Get the inverted column permutation
    std::vector<int> inv_colperm(colperm.size());
    for(int k=0; k<colperm.size(); ++k)
      inv_colperm[colperm[k]] = k;
    
    // Permute the linear system
    SXMatrix Aperm(0,A.size2());
    for(int i=0; i<A.size1(); ++i){
      Aperm.resize(i+1,A.size2());
      for(int el=A.rowind(rowperm[i]); el<A.rowind(rowperm[i]+1); ++el){
        Aperm(i,inv_colperm[A.col(el)]) = A[el];
      }
    }

    // Generate the QR factorization function
    vector<SXMatrix> QR(2);
    qr(Aperm,QR[0],QR[1]);
    SXFunction fact_fcn(A,QR);

    // Optionally generate c code and load as DLL
    if(codegen){
      stringstream ss;
      ss << "symbolic_qr_fact_fcn_" << this;
      fact_fcn_ = dynamicCompilation(fact_fcn,ss.str(),"Symbolic QR factorization function",compiler);
    } else {
      fact_fcn_ = fact_fcn;
    }

    // Initialize factorization function
    fact_fcn_.setOption("name","QR_fact");
    fact_fcn_.init();

    // Symbolic expressions for solve function
    SXMatrix Q = ssym("Q",QR[0].sparsity());
    SXMatrix R = ssym("R",QR[1].sparsity());
    SXMatrix b = ssym("b",input(1).sparsity());
    
    // Solve non-transposed
    // We have inv(A) = inv(Px) * inv(R) * Q' * Pb

    // Permute the right hand side
    SXMatrix bperm(0,b.size2());
    for(int i=0; i<b.size1(); ++i){
      bperm.resize(i+1,b.size2());
      for(int el=b.rowind(rowperm[i]); el<b.rowind(rowperm[i]+1); ++el){
        bperm(i,b.col(el)) = b[el];
      }
    }

    // Solve the factorized system
    SXMatrix xperm = CasADi::solve(R,mul(trans(Q),bperm));

    // Permute back the solution
    SXMatrix x(0,xperm.size2());
    for(int i=0; i<xperm.size1(); ++i){
      x.resize(i+1,xperm.size2());
      for(int el=xperm.rowind(inv_colperm[i]); el<xperm.rowind(inv_colperm[i]+1); ++el){
        x(i,xperm.col(el)) = xperm[el];
      }
    }

    // Generate the QR solve function
    vector<SXMatrix> solv_in(3);
    solv_in[0] = Q;
    solv_in[1] = R;
    solv_in[2] = b;
    SXFunction solv_fcn(solv_in,x);

    // Optionally generate c code and load as DLL
    if(codegen){
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_T_" << this;
      solv_fcn_T_ = dynamicCompilation(solv_fcn,ss.str(),"Symbolic QR solve function",compiler);
    } else {
      solv_fcn_T_ = solv_fcn;
    }

    // Initialize solve function
    solv_fcn_T_.setOption("name","QR_solv_T");
    solv_fcn_T_.init();

    // Solve transposed
    // We have inv(A)' = inv(Pb) * Q *inv(R') * Px

    // Permute the right hand side
    bperm = SXMatrix(0,b.size2());
    for(int i=0; i<b.size1(); ++i){
      bperm.resize(i+1,b.size2());
      for(int el=b.rowind(colperm[i]); el<b.rowind(colperm[i]+1); ++el){
        bperm(i,b.col(el)) = b[el];
      }
    }

    // Solve the factorized system
    xperm = mul(Q,CasADi::solve(trans(R),bperm));

    // Permute back the solution
    x = SXMatrix(0,xperm.size2());
    for(int i=0; i<xperm.size1(); ++i){
      x.resize(i+1,xperm.size2());
      for(int el=xperm.rowind(inv_rowperm[i]); el<xperm.rowind(inv_rowperm[i]+1); ++el){
        x(i,xperm.col(el)) = xperm[el];
      }
    }

    // Mofify the QR solve function
    solv_fcn = SXFunction(solv_in,x);

    // Optionally generate c code and load as DLL
    if(codegen){
      stringstream ss;
      ss << "symbolic_qr_solv_fcn_N_" << this;
      solv_fcn_N_ = dynamicCompilation(solv_fcn,ss.str(),"QR_solv_N",compiler);
    } else {
      solv_fcn_N_ = solv_fcn;
    }

    // Initialize solve function
    solv_fcn_N_.setOption("name","QR_solvT");
    solv_fcn_N_.init();

    // Allocate storage for QR factorization
    Q_ = DMatrix::zeros(Q.sparsity());
    R_ = DMatrix::zeros(R.sparsity());      
  }

  void SymbolicQRInternal::prepare(){
    // Factorize
    fact_fcn_.setInput(input(0));
    fact_fcn_.evaluate();
    fact_fcn_.getOutput(Q_,0);
    fact_fcn_.getOutput(R_,1);
  }

  void SymbolicQRInternal::solve(double* x, int nrhs, bool transpose){
    // Select solve function
    FX& solv = transpose ? solv_fcn_N_ : solv_fcn_T_;

    // Pass QR factorization
    solv.setInput(Q_,0);
    solv.setInput(R_,1);

    // Solve for all right hand sides
    for(int i=0; i<nrhs; ++i){
      solv.setInput(x,2);
      solv.evaluate();
      solv.getOutput(x);
      x += solv.output().size();
    }
  }
  
  void SymbolicQRInternal::generateDeclarations(std::ostream &stream, const std::string& type, CodeGenerator& gen) const{

    // Generate code for the embedded functions
    gen.addDependency(fact_fcn_);
    gen.addDependency(solv_fcn_N_);
    gen.addDependency(solv_fcn_T_);
  }

  void SymbolicQRInternal::generateBody(std::ostream &stream, const std::string& type, CodeGenerator& gen) const{
    
    // Data structures to hold A, Q and R
    stream << "  static int prepared = 0;" << endl;
    stream << "  static d A[" << input(LINSOL_A).size() << "];" << endl;
    stream << "  static d Q[" << Q_.size() << "];" << endl;
    stream << "  static d R[" << R_.size() << "];" << endl;

    // Store matrix to be factorized and check if up-to-date
    stream << "  int i;" << endl;
    stream << "  for(i=0; i<" << input(LINSOL_A).size() << "; ++i){" << endl;
    stream << "    prepared = prepared && A[i] == x0[i];" << endl;
    stream << "    A[i] = x0[i];" << endl;
    stream << "  }" << endl;

    // Factorize if needed
    int fact_ind = gen.getDependency(fact_fcn_);
    stream << "  if(!prepared){" << endl;
    stream << "    f" << fact_ind << "(A,Q,R);" << endl;
    stream << "    prepared = 1;" << endl;
    stream << "  }" << endl;

    // Solve
    int solv_ind_N = gen.getDependency(solv_fcn_N_);
    int solv_ind_T = gen.getDependency(solv_fcn_T_);
    stream << "  if(*x2==0){" << endl;
    stream << "    f" << solv_ind_N << "(Q,R,x1,r0);" << endl;    
    stream << "  } else {" << endl;
    stream << "    f" << solv_ind_T << "(Q,R,x1,r0);" << endl;    
    stream << "  }" << endl;
  }

} // namespace CasADi

  


