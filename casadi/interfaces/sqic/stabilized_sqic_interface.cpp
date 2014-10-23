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


#include "stabilized_sqic_interface.hpp"

#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

#include "casadi/core/std_vector_tools.hpp"

#include "wsqic.hpp"
#include "casadi/interfaces/sqic/resource_sqic.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_STABILIZEDQPSOLVER_SQIC_EXPORT
  casadi_register_stabilizedqpsolver_sqic(StabilizedQpSolverInternal::Plugin* plugin) {
    plugin->creator = StabilizedSqicInterface::creator;
    plugin->name = "sqic";
    plugin->doc = StabilizedSqicInterface::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_STABILIZEDQPSOLVER_SQIC_EXPORT casadi_load_stabilizedqpsolver_sqic() {
    StabilizedQpSolverInternal::registerPlugin(casadi_register_stabilizedqpsolver_sqic);
  }

  StabilizedSqicInterface* StabilizedSqicInterface::clone() const {
    // Return a deep copy
    StabilizedSqicInterface* node = new StabilizedSqicInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  StabilizedSqicInterface::StabilizedSqicInterface(const std::vector<Sparsity>& st)
    : StabilizedQpSolverInternal(st) {
    is_init_ = false;
  }

  StabilizedSqicInterface::~StabilizedSqicInterface() {
    sqicDestroy();
  }

  void StabilizedSqicInterface::evaluate() {
    if (inputs_check_) checkInputs();

    std::copy(input(STABILIZED_QP_SOLVER_X0).begin(),
              input(STABILIZED_QP_SOLVER_X0).end(), x_.begin());
    std::fill(x_.begin()+n_, x_.end(), 0);

    std::transform(input(STABILIZED_QP_SOLVER_LAM_X0).begin(),
                   input(STABILIZED_QP_SOLVER_LAM_X0).end(), rc_.begin(), negate<double>());
    std::fill(rc_.begin()+n_, rc_.end(), 0);

    std::copy(input(STABILIZED_QP_SOLVER_LBX).begin(), input(STABILIZED_QP_SOLVER_LBX).end(),
              bl_.begin());
    std::copy(input(STABILIZED_QP_SOLVER_UBX).begin(), input(STABILIZED_QP_SOLVER_UBX).end(),
              bu_.begin());

    std::copy(input(STABILIZED_QP_SOLVER_LBA).begin(), input(STABILIZED_QP_SOLVER_LBA).end(),
              bl_.begin()+n_);
    std::copy(input(STABILIZED_QP_SOLVER_UBA).begin(), input(STABILIZED_QP_SOLVER_UBA).end(),
              bu_.begin()+n_);

    std::copy(input(STABILIZED_QP_SOLVER_MUE).begin(), input(STABILIZED_QP_SOLVER_MUE).end(),
              piE_.begin());


    for (int i=0;i<n_+nc_+1;++i) {
      if (bl_[i]==-std::numeric_limits<double>::infinity()) bl_[i]=-inf_;
      if (bu_[i]==std::numeric_limits<double>::infinity()) bu_[i]=inf_;
    }

    formatA_.setInput(input(STABILIZED_QP_SOLVER_A), 0);
    formatA_.setInput(input(STABILIZED_QP_SOLVER_G), 1);
    formatA_.evaluate();

    int m = nc_+1;

    sqicSolveStabilized(&output(QP_SOLVER_COST).data()[0],
                        &input(STABILIZED_QP_SOLVER_MU).data()[0],
                        &m, &piE_[0]);

    std::copy(x_.begin(), x_.begin()+n_, output(QP_SOLVER_X).begin());
    std::transform(rc_.begin(), rc_.begin()+n_, output(QP_SOLVER_LAM_X).begin(),
                   negate<double>());
    std::transform(rc_.begin()+n_, rc_.begin()+n_+nc_, output(QP_SOLVER_LAM_A).begin(),
                   negate<double>());

    output(QP_SOLVER_COST)[0]+= x_[n_+nc_];
  }

  void StabilizedSqicInterface::init() {
    // Call the init method of the base class
    StabilizedQpSolverInternal::init();

    if (is_init_) sqicDestroy();

    inf_ = 1.0e+20;

    // Allocate data structures for SQIC
    bl_.resize(n_+nc_+1, 0);
    bu_.resize(n_+nc_+1, 0);
    x_.resize(n_+nc_+1, 0);
    hs_.resize(n_+nc_+1, 0);
    hEtype_.resize(n_+nc_+1, 0);
    pi_.resize(nc_+1, 0);
    piE_.resize(nc_+1, 0);
    rc_.resize(n_+nc_+1, 0);

    locH_ = st_[QP_STRUCT_H].colind();
    indH_ = st_[QP_STRUCT_H].row();

    // Fortran indices are one-based
    for (int i=0;i<indH_.size();++i) indH_[i]+=1;
    for (int i=0;i<locH_.size();++i) locH_[i]+=1;

    // Sparsity of augmented linear constraint matrix
    Sparsity A_ = vertcat(st_[QP_STRUCT_A], Sparsity::dense(1, n_));
    locA_ = A_.colind();
    indA_ = A_.row();

    // Fortran indices are one-based
    for (int i=0;i<indA_.size();++i) indA_[i]+=1;
    for (int i=0;i<locA_.size();++i) locA_[i]+=1;

    // helper functions for augmented linear constraint matrix
    MX a = MX::sym("A", st_[QP_STRUCT_A]);
    MX g = MX::sym("g", n_);
    std::vector<MX> ins;
    ins.push_back(a);
    ins.push_back(g);
    formatA_ = MXFunction(ins, vertcat(a, g.T()));
    formatA_.init();

    // Set objective row of augmented linear constraints
    bu_[n_+nc_] = inf_;
    bl_[n_+nc_] = -inf_;

    is_init_ = true;

    int n = n_;
    int m = nc_+1;

    int nnzA=formatA_.output().size();
    int nnzH=input(STABILIZED_QP_SOLVER_H).size();

    std::fill(hEtype_.begin()+n_, hEtype_.end(), 3);

    sqic(&m , &n, &nnzA, &indA_[0], &locA_[0], &formatA_.output().data()[0], &bl_[0], &bu_[0],
         &hEtype_[0], &hs_[0], &x_[0], &pi_[0], &rc_[0], &nnzH, &indH_[0], &locH_[0],
         &input(STABILIZED_QP_SOLVER_H).data()[0]);
  }

  map<int, string> StabilizedSqicInterface::calc_flagmap() {
    map<int, string> f;

    return f;
  }

  map<int, string> StabilizedSqicInterface::flagmap = StabilizedSqicInterface::calc_flagmap();

  void StabilizedSqicInterface::sqic_error(const string& module, int flag) {
    // Find the error
    map<int, string>::const_iterator it = flagmap.find(flag);

    stringstream ss;
    if (it == flagmap.end()) {
      ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
    } else {
      ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
    }
    ss << " Consult SQIC documentation.";
    casadi_error(ss.str());
  }

  void StabilizedSqicInterface::generateNativeCode(std::ostream& file) const {

    // Dump the contents of resource_sqic, but filter out the C bind stuff
    std::string resource_sqic_input(resource_sqic);
    std::istringstream stream(resource_sqic_input);
    std::string line;
    while (std::getline(stream, line)) {
      size_t b_i = line.find("bind ( C, ");
      if (b_i!=std::string::npos) {
        file << line.substr(0, b_i) << std::endl;
      } else {
        file << line << std::endl;
      }
    }

    file.precision(std::numeric_limits<double>::digits10+2);
    file << std::scientific; // This is really only to force a decimal dot,
    // would be better if it can be avoided

    file << "program exported" << std::endl;
    file << "  use SQICModule" << std::endl;
    file << "  implicit none" << std::endl;
    file << "  integer(ip)               :: m, n, nInf, nnH, nnzH, nnzA, nS, lenpi" << std::endl;


    file << "  real(rp)                  :: Obj, mu" << std::endl;

    file << "  real(rp), allocatable:: bl(:), bu(:), x(:), valA(:), valH(:) , pi(:), piE(:), rc(:)"
         << std::endl;
    file << "  integer(ip), allocatable:: indA(:), locA(:), indH(:), locH(:), hEtype(:), hs(:)"
         << std::endl;

    int n = n_;
    int m = nc_+1;
    int nnzA=formatA_.output().size();
    int nnzH=input(STABILIZED_QP_SOLVER_H).size();

    file << "  n = " << n << std::endl;
    file << "  m = " << m << std::endl;
    file << "  nnzA = " << nnzA << std::endl;
    file << "  nnzH = " << nnzH << std::endl;

    file << "  allocate ( bl(n+m), bu(n+m) )" << std::endl;
    file << "  allocate ( hEtype(n+m) )" << std::endl;
    file << "  allocate ( locA(n+1), valA(nnzA), indA(nnzA) )" << std::endl;
    file << "  allocate ( pi(m), piE(m), rc(n+m), x(n+m) )" << std::endl;
    file << "  allocate ( hs(n+m) )" << std::endl;
    file << "  allocate ( valH(nnzH), locH(n+1), indH(nnzH) )" << std::endl;

    for (int i=0;i<indA_.size();++i) {
      file << "  indA(" << i +1 << ") = " << indA_[i] << std::endl;
    }
    for (int i=0;i<locA_.size();++i) {
      file << "  locA(" << i +1 << ") = " << locA_[i] << std::endl;
    }
    for (int i=0;i<formatA_.output().size();++i) {
      file << "  valA(" << i +1 << ") = " << formatA_.output().at(i) << std::endl;
    }
    for (int i=0;i<bl_.size();++i) {
      file << "  bl(" << i +1 << ") = " << bl_[i] << std::endl;
      file << "  bu(" << i +1 << ") = " << bu_[i] << std::endl;
    }
    for (int i=0;i<hEtype_.size();++i) {
      file << "  hEtype(" << i +1 << ") = " << hEtype_[i] << std::endl;
    }
    for (int i=0;i<hs_.size();++i) {
      file << "  hs(" << i +1 << ") = " << hs_[i] << std::endl;
    }
    for (int i=0;i<indH_.size();++i) {
      file << "  indH(" << i +1 << ") = " << indH_[i] << std::endl;
    }
    for (int i=0;i<locH_.size();++i) {
      file << "  locH(" << i +1 << ") = " << locH_[i] << std::endl;
    }
    for (int i=0;i<input(STABILIZED_QP_SOLVER_H).size();++i) {
      file << "  valH(" << i +1 << ") = " << input(STABILIZED_QP_SOLVER_H).at(i) << std::endl;
    }
    for (int i=0;i<input(QP_SOLVER_X0).size();++i) {
      file << "  x(" << i +1 << ") = " << input(QP_SOLVER_X0).at(i) << std::endl;
    }
    for (int i=0;i<pi_.size();++i) {
      file << "  pi(" << i +1 << ") = " <<  0 << std::endl; //pi_[i] << std::endl;
    }
    for (int i=0;i<rc_.size();++i) {
      file << "  rc(" << i +1 << ") = "
           << ((i<input(QP_SOLVER_LAM_X0).size()) ? -input(QP_SOLVER_LAM_X0).at(i) : 0.0)
           << std::endl;
    }
    file << "  lenpi = " << m << std::endl;
    file << "  mu = " << input(STABILIZED_QP_SOLVER_MUR).at(0) << std::endl;
    for (int i=0;i<piE_.size();++i) {
      file << "  piE(" << i +1 << ") = " << piE_[i] << std::endl;
    }

    file << "  call wsqic (m, n, nnzA, indA, locA, valA, bl, bu, hEtype, hs, x, "
         << "pi, rc, nnzH, indH, locH, valH)" << std::endl;
    /**for (int i=0;i<input(QP_SOLVER_X0).size();++i) {
       file << "  x(" << i +1 << ") = " << input(QP_SOLVER_X0).at(i) << std::endl;
       }
       for (int i=0;i<pi_.size();++i) {
       file << "  pi(" << i +1 << ") = " << pi_[i] << std::endl;
       }
       for (int i=0;i<rc_.size();++i) {
       file << "  rc(" << i +1 << ") = "
       << ((i<input(QP_SOLVER_LAM_X0).size()) ? -input(QP_SOLVER_LAM_X0).at(i) : 0.0)
       << std::endl;
       }*/
    file << "  call sqicSolveStabilized (Obj, mu, lenpi, piE)" << std::endl;
    /**for (int i=0;i<input(QP_SOLVER_X0).size();++i) {
       file << "  x(" << i +1 << ") = " << input(QP_SOLVER_X0).at(i) << std::endl;
       }
       for (int i=0;i<pi_.size();++i) {
       file << "  pi(" << i +1 << ") = " << pi_[i] << std::endl;
       }
       for (int i=0;i<rc_.size();++i) {
       file << "  rc(" << i +1 << ") = "
       << ((i<input(QP_SOLVER_LAM_X0).size()) ? -input(QP_SOLVER_LAM_X0).at(i) : 0.0)
       << std::endl;
       }
       file << "  call sqicSolveStabilized (Obj, mu, lenpi, piE)" << std::endl;**/
    file << "  deallocate ( bl, bu )" << std::endl;
    file << "  deallocate ( hEtype )" << std::endl;
    file << "  deallocate ( locA, valA, indA )" << std::endl;
    file << "  deallocate ( pi, piE, rc, x )" << std::endl;
    file << "  deallocate ( valH, locH, indH )" << std::endl;
    file << "  call sqicDestroy()" << std::endl;
    file << "end program exported" << std::endl;
  }

} // namespace casadi
