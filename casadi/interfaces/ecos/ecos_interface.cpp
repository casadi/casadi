/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#include "ecos_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SOCPSOLVER_ECOS_EXPORT
  casadi_register_socpsolver_ecos(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = EcosInterface::creator;
    plugin->name = "ecos";
    plugin->doc = EcosInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_ECOS_EXPORT casadi_load_socpsolver_ecos() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_ecos);
  }

  EcosInterface* EcosInterface::clone() const {
    // Return a deep copy
    EcosInterface* node = new EcosInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  EcosInterface::EcosInterface(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    // Define options
  }

  EcosInterface::~EcosInterface() {
    if (is_init_) {
      // Destructor
    }
  }

  void EcosInterface::init() {
    // Initialize the base classes
    SocpSolverInternal::init();
    log("EcosInterface::init", "Enter");

    // Initialize c
    ecos_c_vec_.resize(m_+N_+2*nc_+2*n_);

    // Initialize q
    for (int i=0;i<m_;++i) ecos_q_vec_.push_back(ni_[i]+1);

    // Initalize h
    for (int i=0;i<m_+N_;++i) ecos_h_vec_.push_back(0.0);

    // Initialize G
    DMatrix G = DMatrix(m_+N_,m_+N_+2*nc_+2*n_);
    int sum_ni = 0;
    for (int i=0;i<m_;++i) {
      G(sum_ni+i,i) = -1.0;
      for (int ii=0;ii<ni_[i];++ii) {
        G(sum_ni+i+ii+1,m_+sum_ni+i+ii) = -1.0;  
      }
      sum_ni += ni_[i];
    }
    std::cout << G << std::endl;
    const Sparsity& Gsparse = G.sparsity();
    ecos_Gpr_vec_ = std::vector<pfloat>(G.data());
    
    ecos_Gjc_vec_ = std::vector<idxint>(Gsparse.getColind().begin(), Gsparse.getColind().end());
    ecos_Gir_vec_ = std::vector<idxint>(Gsparse.getRow().begin(), Gsparse.getRow().end());
   
  }

  void EcosInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // ECOS expects the SOCP in conic form, obtain this form by formulating dual SOCP
    convertToDualSocp();

    /** Assemble input for ECOS solver
             n (idxint)    : number of variables
             m (idxint)    : number of inequality constraints
             p (idxint)    : number of equality constraints
             l (idxint)    : dimension of positive orthant
        ncones (idxint)    : number of second order cones
             q (idxint[])  : dimension of each cone
           Gpr (pfloat[])  : matrix G in CCS format, data
           Gjc (idxint[])  : matrix G in CCS format, column indices
           Gir (idxint[])  : matrix G in CCS format, row index arrays
           Apr (pfloat[])  : matrix A in CCS format, data
           Ajc (idxint[])  : matrix A in CCS format, column indices
           Air (idxint[])  : matrix A in CCS format, row index arrays
             c (pfloat[])  : vector c (dense)
             h (pfloat[])  : vector h (dense)
             b (pfloat[])  : vector b (dense)
       
        For problem:
          min c'*x
         s.t. A*x = b
              G*x <=_K h
       
        where the last inequality is generalized, i.e. h - G*x belongs to the cone K. ECOS supports the
        positive orthant R_+ and second-order cones Q_n defined as
             Q_n = { (t,x) | t >= || x ||_2 }
       
        In the definition above, t is a scalar and x is in R_{n-1}. The cone K is therefore a direct product of
        the positive orthant and second-order cones:
             K = R_+ x Q_n1 x ... x Q_nN
    */

    std::vector<idxint> Ajc_vec = std::vector<idxint>(dual_A_colind_.begin(),dual_A_colind_.end());
    std::vector<idxint> Air_vec = std::vector<idxint>(dual_A_row_.begin(),dual_A_row_.end());
    std::transform(dual_c_.begin(), dual_c_.end(), ecos_c_vec_.begin(), std::negate<pfloat>());

    idxint n        = dual_n_;
    idxint m        = m_+N_;
    idxint p        = dual_nc_;
    idxint l        = 0;
    idxint ncones   = m_;
    idxint* q       = &ecos_q_vec_[0];
    pfloat* Gpr     = &ecos_Gpr_vec_[0];
    idxint* Gjc     = &ecos_Gjc_vec_[0];
    idxint* Gir     = &ecos_Gir_vec_[0];
    pfloat* Apr     = &dual_A_data_[0];
    idxint* Ajc     = &Ajc_vec[0];
    idxint* Air     = &Air_vec[0];
    pfloat* c       = &ecos_c_vec_[0];
    pfloat* h       = &ecos_h_vec_[0];
    pfloat* b       = &dual_b_[0];


    // Setup ECOS for new problem based on problem definition. We should be able to place 
    // this in init(), it requires implementing an extra feature to cleanup ECOS properly.
    pwork* w =  ECOS_setup(n, m, p, l, ncones, q, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);

    // Solve problem with ECOS
    idxint ret = ECOS_solve(w);

    std::cout << "ECOS primal cost:" << w->best_info->pcost << std::endl;
    std::cout << "ECOS dual cost:" << w->best_info->dcost << std::endl;


    // Cleanup memory
    ECOS_cleanup(w,0);

  }


} // namespace casadi
