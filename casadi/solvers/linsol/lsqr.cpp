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


#include "lsqr.hpp"

#ifdef WITH_DL
#include <cstdlib>
#endif // WITH_DL

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LINSOL_LSQR_EXPORT
  casadi_register_linsol_lsqr(LinsolInternal::Plugin* plugin) {
    plugin->creator = Lsqr::creator;
    plugin->name = "lsqr";
    plugin->doc = Lsqr::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Lsqr::options_;
    return 0;
  }

  extern "C"
  void CASADI_LINSOL_LSQR_EXPORT casadi_load_linsol_lsqr() {
    LinsolInternal::registerPlugin(casadi_register_linsol_lsqr);
  }

  Lsqr::Lsqr(const std::string& name) :
    LinsolInternal(name) {
  }

  Lsqr::~Lsqr() {
    clear_mem();
  }

  Options Lsqr::options_
  = {{&FunctionInternal::options_},
    {{"codegen",
      {OT_BOOL,
       "C-code generation"}},
      {"compiler",
      {OT_STRING,
       "Compiler command to be used for compiling generated code"}}
   }
  };

  void Lsqr::init(const Dict& opts) {
    // Call the base class initializer
    LinsolInternal::init(opts);

    // Default options
    bool codegen = false;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="codegen") {
        codegen = op.second;
      } else if (op.first=="compiler") {
        casadi_error("Option \"compiler\" has been removed");
      }
    }

    // Codegen options
    if (codegen) {
      fopts_["compiler"] = compilerplugin_;
      fopts_["jit_options"] = jit_options_;
    }
  }

  int Lsqr::init_mem(void* mem) const {
    return LinsolInternal::init_mem(mem);
  }

  void Lsqr::reset(void* mem, const int* sp) const {
    LinsolInternal::reset(mem, sp);
    auto m = static_cast<LsqrMemory*>(mem);

    Sparsity spA = Sparsity::compressed(m->sparsity);

    int m_ = sp[0];
    int n_ = sp[1];

    // Temporary storage
    m->w.resize(m_+4*n_);
    m->A.resize(spA.nnz());

  }

  void Lsqr::factorize(void* mem, const double* A) const {
    auto m = static_cast<LsqrMemory*>(mem);

    std::copy(A, A+m->A.size(), get_ptr(m->A));
  }


  void sym_ortho(double a, double b, double& cs, double&sn, double& rho) {
    if (b == 0) {
      cs = sign(a);
      sn = 0;
      rho = fabs(a);
    } else if (a==0) {
      cs = 0;
      sn = sign(b);
      rho = fabs(b);
    } else if (fabs(b)>fabs(a)) {
      double tau = a/b;
      sn = sign(b)/sqrt(1+tau*tau);
      cs = sn*tau;
      rho = b/sn;
    } else {
      double tau = b/a;
      cs = sign(a)/sqrt(1+tau*tau);
      sn = cs*tau;
      rho = a/cs;
    }
  }

  void solve_(void* mem, double* x, bool tr) {
    auto m = static_cast<LsqrMemory*>(mem);

    const double*A = get_ptr(m->A);

    int m_ = m->sparsity[0];
    int n_ = m->sparsity[1];

    double damp = 0;
    double atol=1e-15;
    double btol=1e-15;
    double conlim=1e8;
    int iter_lim = 10000;


    double* w = get_ptr(m->w);

    int itn = 0;
    int istop = 0;
    //int nstop = 0;
    double ctol = 0;
    if (conlim > 0) ctol = 1/conlim;
    double anorm = 0;
    double acond = 0;
    double dampsq = damp*damp;
    double ddnorm = 0;
    double res2 = 0;
    double xnorm = 0;
    double xxnorm = 0;
    double z = 0;
    double cs2 = -1;
    double sn2 = 0;

    double *u = w;  w+= m_; std::copy(x, x+m_, u);
    double *v = w;  w+= n_; fill_n(v, n_, 0.0);
    double *xx = w; w+= n_; fill_n(xx, n_, 0.0);
    double *ww = w; w+= n_; fill_n(v, n_, 0.0);
    double *dk = w; w+= n_;

    double alpha = 0;
    double beta = casadi_norm_2(m_, u);

    if (beta>0) {
      for (int i=0;i<m_;++i) u[i]*=1/beta;
      casadi_mv(A, get_ptr(m->sparsity), u, v, !tr);
      alpha = casadi_norm_2(n_, v);
    }

    if (alpha>0) {
      for (int i=0;i<n_;++i) v[i]*=1/alpha;
      std::copy(v, v+n_, ww);
    }

    double rhobar = alpha;
    double phibar = beta;
    double bnorm = beta;
    double rnorm = beta;
    //double r1norm = rnorm;
    //double r2norm = rnorm;

    double arnorm = alpha * beta;


  //  userOut() << "   Itn      x[0]       r1norm     r2norm "
  //               "Compatible    LS      Norm A   Cond A" << std::endl;

  //  double test1 = 1;
  //  double test2 = alpha / beta;


  //  userOut() << itn << ":" << xx[0] << ":" << r1norm << ":" <<
  //r2norm << ":" << test1 << ":" << test2 << std::endl;

    while (itn<iter_lim) {
      itn++;
      for (int i=0;i<m_;++i) u[i]*=-alpha;
      casadi_mv(A, get_ptr(m->sparsity), v, u, tr);
      beta = casadi_norm_2(m_, u);

      if (beta>0) {
        for (int i=0;i<m_;++i) u[i]*=1/beta;
        anorm = sqrt(anorm*anorm + alpha*alpha+beta*beta+damp*damp);
        for (int i=0;i<n_;++i) v[i]*=-beta;
        casadi_mv(A, get_ptr(m->sparsity), u, v, !tr);
        alpha = casadi_norm_2(n_, v);
        if (alpha>0) for (int i=0;i<n_;++i) v[i]*=1/alpha;
      }

      double rhobar1 = sqrt(rhobar*rhobar+damp*damp);

      double cs1 = rhobar / rhobar1;
      double sn1 = damp / rhobar1;
      double psi = sn1 * phibar;
      phibar *= cs1;

      double cs, sn, rho;
      sym_ortho(rhobar1, beta, cs, sn, rho);

      double theta = sn * alpha;
      rhobar = -cs * alpha;
      double phi = cs * phibar;
      phibar *= sn;
      double tau = sn * phi;

      double t1 = phi / rho;
      double t2 = -theta / rho;

      for (int i=0;i<n_;++i) dk[i]=ww[i]/rho;

      for (int i=0; i<n_; ++i) xx[i] += t1*ww[i];
      for (int i=0; i<n_; ++i) ww[i] = v[i] + t2*ww[i];

      double n2dk = casadi_norm_2(n_, dk);
      ddnorm += n2dk*n2dk;

      double delta = sn2 * rho;
      double gambar = -cs2 * rho;
      double rhs = phi - delta * z;
      double zbar = rhs / gambar;
      xnorm = sqrt(xxnorm + zbar*zbar);
      double gamma = sqrt(gambar*gambar + theta*theta);
      cs2 = gambar / gamma;
      sn2 = theta / gamma;
      z = rhs / gamma;
      xxnorm += z*z;

      acond = anorm * sqrt(ddnorm);
      double res1 = phibar*phibar;
      res2 += psi*psi;
      rnorm = sqrt(res1+res2);
      arnorm = alpha*fabs(tau);

      double r1sq = rnorm*rnorm - dampsq * xxnorm;
      double r1norm = sqrt(fabs(r1sq));
      if (r1sq < 0) r1norm = -r1norm;
      //double r2norm = rnorm;

      double test1 = rnorm / bnorm;
      double test2 = arnorm / (anorm * rnorm);
      double test3 = 1 / acond;
      t1 = test1 / (1 + anorm * xnorm / bnorm);
      double rtol = btol + atol * anorm * xnorm / bnorm;

      if (itn >= iter_lim) istop = 7;
      if (1 + test3 <= 1) istop = 6;
      if (1 + test2 <= 1) istop = 5;
      if (1 + t1 <= 1) istop = 4;

      if (test3 <= ctol) istop = 3;
      if (test2 <= atol) istop = 2;
      if (test1 <= rtol) istop = 1;

      if (istop != 0) break;

    }
    std::copy(xx, xx+m_, x);

  }

  void Lsqr::solve(void* mem, double* x, int nrhs, bool tr) const {
    auto m = static_cast<LsqrMemory*>(mem);

    int n_ = m->sparsity[1];

    for (int i=0;i<nrhs;++i) {
      solve_(mem, x+i*n_, tr);
    }

  }


} // namespace casadi
