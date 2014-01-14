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

#include "irk_integrator_internal.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

  IRKIntegratorInternal::IRKIntegratorInternal(const FX& f, const FX& g) : RKBaseInternal(f,g){
    addOption("interpolation_order",           OT_INTEGER,  3,  "Order of the interpolating polynomials");
    addOption("collocation_scheme",            OT_STRING,  "radau",  "Collocation scheme","radau|legendre");
    addOption("implicit_solver",               OT_IMPLICITFUNCTION,  GenericType(), "An implicit function solver");
    addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
    setOption("name","unnamed_irk_integrator");
  }

  void IRKIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    RKBaseInternal::deepCopyMembers(already_copied);
    implicit_solver_ = deepcopy(implicit_solver_,already_copied);
    explicit_fcn_ = deepcopy(explicit_fcn_,already_copied);
  }

  IRKIntegratorInternal::~IRKIntegratorInternal(){
  }

  void IRKIntegratorInternal::init(){
  
    // Call the base class init
    RKBaseInternal::init();
  
    // Interpolation order
    int deg = getOption("interpolation_order");

    // All collocation time points
    std::vector<long double> tau_root = collocationPointsL(deg,getOption("collocation_scheme"));

    // Coefficients of the collocation equation
    vector<vector<double> > C(deg+1,vector<double>(deg+1,0));
      
    // Coefficients of the continuity equation
    vector<double> D(deg+1,0);
      
    // Coefficients of the quadratures
    vector<double> B(deg+1,0);

    // For all collocation points
    for(int j=0; j<deg+1; ++j){

      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      Polynomial p = 1;
      for(int r=0; r<deg+1; ++r){
        if(r!=j){
          p *= Polynomial(-tau_root[r],1)/(tau_root[j]-tau_root[r]);
        }
      }
    
      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      D[j] = zeroIfSmall(p(1.0L));
    
      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      Polynomial dp = p.derivative();
      for(int r=0; r<deg+1; ++r){
        C[j][r] = zeroIfSmall(dp(tau_root[r]));
      }
        
      // Integrate polynomial to get the coefficients of the quadratures
      Polynomial ip = p.anti_derivative();
      B[j] = zeroIfSmall(ip(1.0L));
    }

    // Symbolic inputs
    MX x0 = msym("x0",f_.input(DAE_X).sparsity());
    MX p = msym("p",f_.input(DAE_P).sparsity());
    MX t = msym("t",f_.input(DAE_T).sparsity());

    // Implicitly defined variables (z and x)
    MX v = msym("v",deg*(nx_+nz_));
    vector<int> v_offset(1,0);
    for(int d=0; d<deg; ++d){
      v_offset.push_back(v_offset.back()+nx_);
      v_offset.push_back(v_offset.back()+nz_);
    }
    vector<MX> vv = vertsplit(v,v_offset);
    vector<MX>::const_iterator vv_it = vv.begin();

    // Collocated states
    vector<MX> x(deg+1), z(deg+1);
    x[0] = x0;
    for(int d=1; d<=deg; ++d){
      x[d] = *vv_it++;
      z[d] = *vv_it++;
    }

    // Collocation time points
    vector<MX> tt(deg+1);
    for(int d=0; d<=deg; ++d){
      tt[d] = t + h_*tau_root[d];
    }

    // Equations that implicitly define v
    vector<MX> eq;

    // Quadratures
    MX qf = MX::zeros(f_.output(DAE_QUAD).sparsity());

    // End state
    MX xf = D[0]*x0;

    // For all collocation points
    for(int j=1; j<deg+1; ++j){

      // Get an expression for the state derivative at the collocation point
      MX xp_j = 0;
      for(int r=0; r<deg+1; ++r){
        xp_j += (C[r][j]/h_) * x[r];
      }

      // Evaluate the DAE
      vector<MX> f_in(DAE_NUM_IN);
      f_in[DAE_T] = tt[j];
      f_in[DAE_P] = p;
      f_in[DAE_X] = x[j];
      f_in[DAE_Z] = z[j];
      vector<MX> f_out = f_.call(f_in);
      
      // Add collocation equation
      eq.push_back(f_out[DAE_ODE] - xp_j);
        
      // Add the algebraic conditions
      eq.push_back(f_out[DAE_ALG]);

      // Add contribution to the final state
      xf += D[j]*x[j];
        
      // Add contribution to quadratures
      qf += (B[j]*h_)*f_out[DAE_QUAD];
    }

    // Form forward discrete time dynamics
    vector<MX> F_in(DAE_NUM_IN);
    F_in[DAE_T] = t;
    F_in[DAE_X] = x0;
    F_in[DAE_P] = p;
    F_in[DAE_Z] = v;
    vector<MX> F_out(DAE_NUM_OUT);
    F_out[DAE_ODE] = xf;
    F_out[DAE_ALG] = vertcat(eq);
    F_out[DAE_QUAD] = qf;
    F_ = MXFunction(F_in,F_out);
    F_.init();
    
    

    {

    // MX version of the same
    MX h_mx = h_;
    
    // Coefficients of the collocation equation
    vector<vector<MX> > C(deg+1,vector<MX>(deg+1));
  
    // Coefficients of the collocation equation as DMatrix
    DMatrix C_num = DMatrix::zeros(deg+1,deg+1);

    // Coefficients of the continuity equation
    vector<MX> D(deg+1);
  
    // Coefficients of the collocation equation as DMatrix
    DMatrix D_num = DMatrix(deg+1,1,0);

    // Collocation point
    SXMatrix tau = ssym("tau");

    // For all collocation points
    for(int j=0; j<deg+1; ++j){
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      SXMatrix L = 1;
      for(int j2=0; j2<deg+1; ++j2){
        if(j2 != j){
          L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2]);
        }
      }
    
      SXFunction lfcn(tau,L);
      lfcn.init();
  
      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      lfcn.setInput(1.0);
      lfcn.evaluate();
      D[j] = lfcn.output();
      D_num(j) = lfcn.output();

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      FX tfcn = lfcn.tangent();
      tfcn.init();
      for(int j2=0; j2<deg+1; ++j2){
        tfcn.setInput(tau_root[j2]);        
        tfcn.evaluate();
        C[j][j2] = tfcn.output();
        C_num(j,j2) = tfcn.output();
      }
    }


    C_num(std::vector<int>(1,0),ALL) = 0;
    C_num(0,0)   = 1;
  
    // Coefficients of the quadrature
    DMatrix Q = solve(C_num,D_num);
    
    Q = B;
    //   cout << "Q = " << Q << endl;
  

    casadi_assert_message(fabs(sumAll(Q)-1).at(0)<1e-9,"Check on quadrature coefficients");
    casadi_assert_message(fabs(sumAll(D_num)-1).at(0)<1e-9,"Check on collocation coefficients");
  
    // Initial state
    MX X0("X0",nx_);
  
    // Parameters
    MX P("P",np_);
  
    // Backward state
    MX RX0("RX0",nrx_);
  
    // Backward parameters
    MX RP("RP",nrp_);
  
    // Collocated differential states and algebraic variables
    int nX = (nk_*(deg+1)+1)*(nx_+nrx_);
    int nZ = nk_*deg*(nz_+nrz_);
  
    // Unknowns
    MX V("V",nX+nZ);
    int offset = 0;
  
    // Get collocated states, algebraic variables and times
    vector<vector<MX> > X(nk_+1);
    vector<vector<MX> > RX(nk_+1);
    vector<vector<MX> > Z(nk_);
    vector<vector<MX> > RZ(nk_);
    coll_time_.resize(nk_+1);
    for(int k=0; k<nk_+1; ++k){
      // Number of time points
      int nj = k==nk_ ? 1 : deg+1;
    
      // Allocate differential states expressions at the time points
      X[k].resize(nj);
      RX[k].resize(nj);
      coll_time_[k].resize(nj);

      // Allocate algebraic variable expressions at the collocation points
      if(k!=nk_){
        Z[k].resize(nj-1);
        RZ[k].resize(nj-1);
      }

      // For all time points
      for(int j=0; j<nj; ++j){
        // Get expressions for the differential state
        X[k][j] = V[range(offset,offset+nx_)];
        offset += nx_;
        RX[k][j] = V[range(offset,offset+nrx_)];
        offset += nrx_;
      
        // Get the local time
        coll_time_[k][j] = t0_ + h_*(k + tau_root[j]);
      
        // Get expressions for the algebraic variables
        if(j>0){
          Z[k][j-1] = V[range(offset,offset+nz_)];
          offset += nz_;
          RZ[k][j-1] = V[range(offset,offset+nrz_)];
          offset += nrz_;
        }
      }
    }
  
    // Check offset for consistency
    casadi_assert(offset==V.size());

    // Constraints
    vector<MX> g;
    g.reserve(2*(nk_+1));
  
    // Quadrature expressions
    MX QF = MX::zeros(nq_);
    MX RQF = MX::zeros(nrq_);
  
    // Counter
    int jk = 0;
  
    // Add initial condition
    g.push_back(X[0][0]-X0);
  
    // For all finite elements
    for(int k=0; k<nk_; ++k, ++jk){
  
      // For all collocation points
      for(int j=1; j<deg+1; ++j, ++jk){
        // Get the time
        MX tkj = coll_time_[k][j];
      
        // Get an expression for the state derivative at the collocation point
        MX xp_jk = 0;
        for(int j2=0; j2<deg+1; ++j2){
          xp_jk += C[j2][j]*X[k][j2];
        }
      
        // Add collocation equations to the NLP
        vector<MX> f_in(DAE_NUM_IN);
        f_in[DAE_T] = tkj;
        f_in[DAE_P] = P;
        f_in[DAE_X] = X[k][j];
        f_in[DAE_Z] = Z[k][j-1];
      
        vector<MX> f_out;
        f_out = f_.call(f_in);
        g.push_back(h_mx*f_out[DAE_ODE] - xp_jk);
      
        // Add the algebraic conditions
        if(nz_>0){
          g.push_back(f_out[DAE_ALG]);
        }
      
        // Add the quadrature
        if(nq_>0){
          QF += Q[j]*h_mx*f_out[DAE_QUAD];
        }
      
        // Now for the backward problem
        if(nrx_>0){
        
          // Get an expression for the state derivative at the collocation point
          MX rxp_jk = 0;
          for(int j2=0; j2<deg+1; ++j2){
            rxp_jk += C[j2][j]*RX[k][j2];
          }
        
          // Add collocation equations to the NLP
          vector<MX> g_in(RDAE_NUM_IN);
          g_in[RDAE_T] = tkj;
          g_in[RDAE_X] = X[k][j];
          g_in[RDAE_Z] = Z[k][j-1];
          g_in[RDAE_P] = P;
          g_in[RDAE_RP] = RP;
          g_in[RDAE_RX] = RX[k][j];
          g_in[RDAE_RZ] = RZ[k][j-1];
        
          vector<MX> g_out;
          g_out = g_.call(g_in);
          g.push_back(h_mx*g_out[RDAE_ODE] + rxp_jk);
        
          // Add the algebraic conditions
          if(nrz_>0){
            g.push_back(g_out[RDAE_ALG]);
          }
        
          // Add the backward quadrature
          if(nrq_>0){
            RQF += Q[j]*h_mx*g_out[RDAE_QUAD];
          }
        }
      }
    
      // Get an expression for the state at the end of the finite element
      MX xf_k = 0;
      for(int j=0; j<deg+1; ++j){
        xf_k += D[j]*X[k][j];
      }

      // Add continuity equation to NLP
      g.push_back(X[k+1][0] - xf_k);
    
      if(nrx_>0){
        // Get an expression for the state at the end of the finite element
        MX rxf_k = 0;
        for(int j=0; j<deg+1; ++j){
          rxf_k += D[j]*RX[k][j];
        }

        // Add continuity equation to NLP
        g.push_back(RX[k+1][0] - rxf_k);
      }
    }
  
    // Add initial condition for the backward integration
    if(nrx_>0){
      g.push_back(RX[nk_][0]-RX0);
    }
  
    // Constraint expression
    MX gv = vertcat(g);
    
    // Make sure that the dimension is consistent with the number of unknowns
    casadi_assert_message(gv.size()==V.size(),"Implicit function unknowns and equations do not match");

    // Implicit function
    vector<MX> ifcn_in(1+INTEGRATOR_NUM_IN);
    ifcn_in[0] = V;
    ifcn_in[1+INTEGRATOR_X0] = X0;
    ifcn_in[1+INTEGRATOR_P] = P;
    ifcn_in[1+INTEGRATOR_RX0] = RX0;
    ifcn_in[1+INTEGRATOR_RP] = RP;
    FX ifcn = MXFunction(ifcn_in,gv);
    std::stringstream ss_ifcn;
    ss_ifcn << "collocation_implicit_residual_" << getOption("name");
    ifcn.setOption("name",ss_ifcn.str());
    ifcn.init(); 
  
    // Auxiliary output function
    vector<MX> afcn_out(1+INTEGRATOR_NUM_OUT);
    afcn_out[0] = V;
    afcn_out[1+INTEGRATOR_XF] = X[nk_][0];
    afcn_out[1+INTEGRATOR_QF] = QF;
    afcn_out[1+INTEGRATOR_RXF] = RX[0][0];
    afcn_out[1+INTEGRATOR_RQF] = RQF;
    FX afcn = MXFunction(ifcn_in,afcn_out);
    std::stringstream ss_afcn;
    ss_afcn << "collocation_output_" << getOption("name");
    afcn.setOption("name",ss_afcn.str());
    afcn.init();
  
    // Get the NLP creator function
    implicitFunctionCreator implicit_function_creator = getOption("implicit_solver");
  
    // Allocate an NLP solver
    implicit_solver_ = implicit_function_creator(ifcn,FX(),LinearSolver());
    std::stringstream ss_implicit_solver;
    ss_implicit_solver << "collocation_implicitsolver_" << getOption("name");
    implicit_solver_.setOption("name",ss_implicit_solver.str());
    
    // Pass options
    if(hasSetOption("implicit_solver_options")){
      const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
      implicit_solver_.setOption(implicit_solver_options);
    }
  
    // Initialize the solver
    implicit_solver_.init();
  
    // Nonlinear constraint function input
    vector<MX> gfcn_in(INTEGRATOR_NUM_IN);
    gfcn_in[INTEGRATOR_X0] = X0;
    gfcn_in[INTEGRATOR_P] = P;
    gfcn_in[INTEGRATOR_RX0] = RX0;
    gfcn_in[INTEGRATOR_RP] = RP;
    ifcn_in[0] = implicit_solver_.call(gfcn_in).front();
    explicit_fcn_ = MXFunction(gfcn_in,afcn.call(ifcn_in));
    std::stringstream ss_explicit_fcn;
    ss_explicit_fcn << "collocation_explicit_" << getOption("name");
    explicit_fcn_.setOption("name",ss_explicit_fcn.str());
    explicit_fcn_.init();

    }
  }
  
  void IRKIntegratorInternal::reset(){
    
    // Call the base class method
    RKBaseInternal::reset();
      
    // Pass the inputs
    for(int iind=0; iind<INTEGRATOR_NUM_IN; ++iind){
      explicit_fcn_.input(iind).set(input(iind));
    }
  
    vector<double>& v = implicit_solver_.output().data();
            
    // Integrate, stopping at all time points
    int offs=0;
    for(int k=0; k<coll_time_.size(); ++k){
      for(int j=0; j<coll_time_[k].size(); ++j){
        
        // Save the differential states
        const DMatrix& x = input(INTEGRATOR_X0);
        for(int i=0; i<nx_; ++i){
          v.at(offs++) = x.at(i);
        }

        // Skip algebraic variables (for now) // FIXME
        if(j>0){
          offs += nz_;
        }
        
        // Skip backward states // FIXME
        const DMatrix& rx = input(INTEGRATOR_RX0);
        for(int i=0; i<nrx_; ++i){
          v.at(offs++) = rx.at(i);
        }
        
        // Skip backward algebraic variables // FIXME
        if(j>0){
          offs += nrz_;
        }
      }
    }
    
    // Solve the system of equations
    explicit_fcn_.evaluate();
  }

  void IRKIntegratorInternal::integrate(double t_out){
    for(int oind=0; oind<INTEGRATOR_NUM_OUT; ++oind){
      output(oind).set(explicit_fcn_.output(1+oind));
    }
  }

  void IRKIntegratorInternal::integrateB(double t_out){
  }

  double IRKIntegratorInternal::zeroIfSmall(double x){
    return fabs(x) < numeric_limits<double>::epsilon() ? 0 : x;
  }

} // namespace CasADi
