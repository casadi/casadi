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

#include "collocation_integrator_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

CollocationIntegratorInternal::CollocationIntegratorInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
  addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  addOption("interpolation_order",           OT_INTEGER,  3,  "Order of the interpolating polynomials");
  addOption("collocation_scheme",            OT_STRING,  "radau",  "Collocation scheme","radau|legendre");
  addOption("implicit_solver",               OT_IMPLICITFUNCTION,  GenericType(), "An implicit function solver");
  addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  addOption("expand_f",                      OT_BOOLEAN,  false, "Expand the ODE/DAE residual function in an SX graph");
  addOption("expand_q",                      OT_BOOLEAN,  false, "Expand the quadrature function in an SX graph");
  addOption("hotstart",                      OT_BOOLEAN,  true, "Initialize the trajectory at the previous solution");
  addOption("quadrature_solver",             OT_LINEARSOLVER,  GenericType(), "An linear solver to solver the quadrature equations");
  addOption("quadrature_solver_options",     OT_DICTIONARY, GenericType(), "Options to be passed to the quadrature solver");
  addOption("startup_integrator",            OT_INTEGRATOR,  GenericType(), "An ODE/DAE integrator that can be used to generate a startup trajectory");
  addOption("startup_integrator_options",    OT_DICTIONARY, GenericType(), "Options to be passed to the startup integrator");
  setOption("name","unnamed_collocation_integrator");
}

void CollocationIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
}

CollocationIntegratorInternal::~CollocationIntegratorInternal(){
}

void CollocationIntegratorInternal::init(){
  
  // Call the base class init
  IntegratorInternal::init();
  
  // Legendre collocation points
  double legendre_points[][6] = {
    {0},
    {0,0.500000},
    {0,0.211325,0.788675},
    {0,0.112702,0.500000,0.887298},
    {0,0.069432,0.330009,0.669991,0.930568},
    {0,0.046910,0.230765,0.500000,0.769235,0.953090}};
    
  // Radau collocation points
  double radau_points[][6] = {
    {0},
    {0,1.000000},
    {0,0.333333,1.000000},
    {0,0.155051,0.644949,1.000000},
    {0,0.088588,0.409467,0.787659,1.000000},
    {0,0.057104,0.276843,0.583590,0.860240,1.000000}};

  // Read options
  bool use_radau;
  if(getOption("collocation_scheme")=="radau"){
    use_radau = true;
  } else if(getOption("collocation_scheme")=="legendre"){
    use_radau = false;
  }
  
  // Hotstart?
  hotstart_ = getOption("hotstart");
  
  // Number of finite elements
  int nk = getOption("number_of_finite_elements");
  
  // Interpolation order
  deg_ = getOption("interpolation_order");

  // All collocation time points
  double* tau_root = use_radau ? radau_points[deg_] : legendre_points[deg_];

  // Size of the finite elements
  double h = (tf_-t0_)/nk;
  
  // MX version of the same
  MX h_mx = h;
    
  // Coefficients of the collocation equation
  vector<vector<MX> > C(deg_+1,vector<MX>(deg_+1));
  
  // Coefficients of the collocation equation as DMatrix
  DMatrix C_num = DMatrix(deg_+1,deg_+1,0);

  // Coefficients of the continuity equation
  vector<MX> D(deg_+1);
  
  // Coefficients of the collocation equation as DMatrix
  DMatrix D_num = DMatrix(deg_+1,1,0);

  // Collocation point
  SXMatrix tau = ssym("tau");
  
  // Collection of polynomial functions
  std::vector< SXFunction > lfcnv;

  // For all collocation points
  for(int j=0; j<deg_+1; ++j){
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    SXMatrix L = 1;
    for(int j2=0; j2<deg_+1; ++j2){
      if(j2 != j){
        L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2]);
      }
    }
    
    SXFunction lfcn(tau,L);
    lfcn.init();
    lfcnv.push_back(lfcn);
  
    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn.setInput(1.0);
    lfcn.evaluate();
    D[j] = lfcn.output();
    D_num(j) = lfcn.output();

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    for(int j2=0; j2<deg_+1; ++j2){
      lfcn.setInput(tau_root[j2]);
      lfcn.setFwdSeed(1.0);
      lfcn.evaluate(1,0);
      C[j][j2] = lfcn.fwdSens();
      C_num(j,j2) = lfcn.fwdSens();
    }
  }

  C_num(std::vector<int>(1,0),ALL) = 0;
  C_num(0,0)   = 1;
  
  // Coefficients of the quadrature
  DMatrix Q = solve(C_num,D_num);
  
  casadi_assert_message(fabs(sumAll(Q)-1)<1e-9,"Check on quadrature coefficients");
  casadi_assert_message(fabs(sumAll(D_num)-1)<1e-9,"Check on collocation coefficients");
  
  // Initial state
  MX X0("X0",nx_);
  
  // Parameters
  MX P("P",np_);
  
  // Backward state
  MX RX0("RX0",nrx_);
  
  // Backward parameters
  MX RP("RP",nrp_);
  
  // Collocated differential states and algebraic variables
  int nX = (nk*(deg_+1)+1)*(nx_+nrx_);
  int nZ = nk*deg_*(nz_+nrz_);
  
  // Unknowns
  MX V("V",nX+nZ);
  int offset = 0;
  
  // Get collocated states, algebraic variables and times
  vector<vector<MX> > X(nk+1);
  vector<vector<MX> > RX(nk+1);
  vector<vector<MX> > Z(nk);
  vector<vector<MX> > RZ(nk);
  coll_time_.resize(nk+1);
  for(int k=0; k<nk+1; ++k){
    // Number of time points
    int nj = k==nk ? 1 : deg_+1;
    
    // Allocate differential states expressions at the time points
    X[k].resize(nj);
    RX[k].resize(nj);
    coll_time_[k].resize(nj);

    // Allocate algebraic variable expressions at the collocation points
    if(k!=nk){
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
      coll_time_[k][j] = t0_ + h*(k + tau_root[j]);
      
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
  
  MX Tau("tau");
  MX rTau("rtau");
  
  std::vector< MX > Lfcn_inputs;
  Lfcn_inputs.push_back(Tau);
  Lfcn_inputs.push_back(rTau);
  int kmul = (deg_+1)*(nx_+nrx_)+(deg_)*(nz_+nrz_);
  
  MX Vk("Vk",kmul);
  MX Vrk("Vrk",kmul);
  
  MX QFk("QFk",nq_*(deg_+1));
  MX RQFk("RQFk",nrq_*(deg_+1));
  
  Lfcn_inputs.push_back(Vk); // Collocation interval
  Lfcn_inputs.push_back(Vrk);
  Lfcn_inputs.push_back(QFk);
  Lfcn_inputs.push_back(RQFk);
  
  std::vector< MX > Lfcn_outputs(INTEGRATOR_NUM_OUT);
  Lfcn_outputs[INTEGRATOR_XF]  = MX::zeros(nx_,1);
  Lfcn_outputs[INTEGRATOR_RXF] = MX::zeros(nrx_,1);
  Lfcn_outputs[INTEGRATOR_QF]  = MX::zeros(nq_,1);
  Lfcn_outputs[INTEGRATOR_RQF] = MX::zeros(nrq_,1);
  
  offset = 0;
  // For all collocation points
  for(int j=0; j<deg_+1; ++j) {
    MX f = lfcnv[j].call(Tau)[0];
    MX rf = lfcnv[j].call(rTau)[0];
    Lfcn_outputs[INTEGRATOR_XF]+=Vk[range(offset,offset+nx_)]*f;
    offset += nx_;
    Lfcn_outputs[INTEGRATOR_RXF]+=Vrk[range(offset,offset+nrx_)]*rf;
    offset += nrx_;
    if (j>0) {
      offset += nz_;
      offset += nrz_;
    }
    if(nq_>0){
      Lfcn_outputs[INTEGRATOR_QF]+=QFk[range(j*nq_,(j+1)*nq_)]*f;
    }
    if(nrq_>0){
      Lfcn_outputs[INTEGRATOR_RQF]+=RQFk[range(j*nrq_,(j+1)*nrq_)]*f;
    }
  }
  
  // Check offset for consistency
  casadi_assert(offset==Vk.size());
  casadi_assert(offset==Vrk.size());
  
  Lfcn_ = MXFunction(Lfcn_inputs,Lfcn_outputs);
  Lfcn_.init();
  
  Lfcn_ = SXFunction(Lfcn_);
  Lfcn_.init();

  // Constraints
  vector<MX> g;
  g.reserve(2*(nk+1));
  
  // Quadrature expressions
  MX QF = MX::zeros(nq_);
  MX RQF = MX::zeros(nrq_);
  
  // Quadratures states values at collocation intervals
  std::vector<MX> Qk(nk*(deg_+1));
  std::vector<MX> RQk(nk*(deg_+1));
  
  // Counter
  int jk = 0;
  
  // Add initial condition
  g.push_back(X[0][0]-X0);
  
  // For all finite elements
  for(int k=0; k<nk; ++k, ++jk){

    if(nrq_>0) RQk.at(k*(deg_+1)) = RQF;
    if(nq_>0)   Qk.at(k*(deg_+1)) =  QF;
        
    // For all collocation points
    for(int j=1; j<deg_+1; ++j, ++jk){
      // Get the time
      MX tkj = coll_time_[k][j];
      
      // Get an expression for the state derivative at the collocation point
      MX xp_jk = 0;
      for(int j2=0; j2<deg_+1; ++j2){
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
        Qk.at(k*(deg_+1)+j) =  QF;
      }
      
      // Now for the backward problem
      if(nrx_>0){
        
        // Get an expression for the state derivative at the collocation point
        MX rxp_jk = 0;
        for(int j2=0; j2<deg_+1; ++j2){
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
          RQk.at(k*(deg_+1)+j) = RQF;
        }

      }
      
    }
    
    // Get an expression for the state at the end of the finite element
    MX xf_k = 0;
    for(int j=0; j<deg_+1; ++j){
      xf_k += D[j]*X[k][j];
    }

    // Add continuity equation to NLP
    g.push_back(X[k+1][0] - xf_k);
    
    if(nrx_>0){
      // Get an expression for the state at the end of the finite element
      MX rxf_k = 0;
      for(int j=0; j<deg_+1; ++j){
        rxf_k += D[j]*RX[k][j];
      }

      // Add continuity equation to NLP
      g.push_back(RX[k+1][0] - rxf_k);
    }
  }
  
  // Add initial condition for the backward integration
  if(nrx_>0){
    g.push_back(RX[nk][0]-RX0);
  }
  
  // Constraint expression
  MX gv = vertcat(g);
    
  // Make sure that the dimension is consistent with the number of unknowns
  casadi_assert_message(gv.size()==V.size(),"Implicit function unknowns and equations do not match");

  // Nonlinear constraint function input
  vector<MX> gfcn_in(1+INTEGRATOR_NUM_IN);
  gfcn_in[0] = V;
  gfcn_in[1+INTEGRATOR_X0] = X0;
  gfcn_in[1+INTEGRATOR_P] = P;
  gfcn_in[1+INTEGRATOR_RX0] = RX0;
  gfcn_in[1+INTEGRATOR_RP] = RP;

  vector<MX> gfcn_out(1+INTEGRATOR_NUM_OUT+2);
  gfcn_out[0] = gv;
  gfcn_out[1+INTEGRATOR_XF] = X[nk][0];
  gfcn_out[1+INTEGRATOR_QF] = QF;
  gfcn_out[1+INTEGRATOR_RXF] = RX[0][0];
  gfcn_out[1+INTEGRATOR_RQF] = RQF;
  gfcn_out[1+INTEGRATOR_NUM_IN+0] = vertcat(Qk);
  gfcn_out[1+INTEGRATOR_NUM_IN+1] = vertcat(RQk);
  
  // Nonlinear constraint function
  FX gfcn = MXFunction(gfcn_in,gfcn_out);
  
  // Expand f?
  bool expand_f = getOption("expand_f");
  if(expand_f){
    gfcn.init();
    gfcn = SXFunction(shared_cast<MXFunction>(gfcn));
  }
  
  // Get the NLP creator function
  implicitFunctionCreator implicit_function_creator = getOption("implicit_solver");
  
  // Allocate an NLP solver
  implicit_solver_ = implicit_function_creator(gfcn);
  
  // Pass options
  if(hasSetOption("implicit_solver_options")){
    const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
    implicit_solver_.setOption(implicit_solver_options);
  }
  
  // Initialize the solver
  implicit_solver_.init();
  
  if(hasSetOption("startup_integrator")){
    
    // Create the linear solver
    integratorCreator startup_integrator_creator = getOption("startup_integrator");
    
    // Allocate an NLP solver
    startup_integrator_ = startup_integrator_creator(f_,g_);
    
    // Pass options
    startup_integrator_.setOption("number_of_fwd_dir",0); // not needed
    startup_integrator_.setOption("number_of_adj_dir",0); // not needed
    startup_integrator_.setOption("t0",coll_time_.front().front());
    startup_integrator_.setOption("tf",coll_time_.back().back());
    if(hasSetOption("startup_integrator_options")){
      const Dictionary& startup_integrator_options = getOption("startup_integrator_options");
      startup_integrator_.setOption(startup_integrator_options);
    }
    
    // Initialize the startup integrator
    startup_integrator_.init();
  }

  // Mark the system not yet integrated
  integrated_once_ = false;
}
  
void CollocationIntegratorInternal::initAdj(){
}

void CollocationIntegratorInternal::reset(int nsens, int nsensB, int nsensB_store){

  // Call the base class method
  IntegratorInternal::reset(nsens,nsensB,nsensB_store);
  
  // Pass the inputs
  for(int iind=0; iind<INTEGRATOR_NUM_IN; ++iind){
    implicit_solver_.input(iind).set(input(iind));
  }
  
  // Pass the forward seeds
  for(int dir=0; dir<nsens; ++dir){
    for(int iind=0; iind<INTEGRATOR_NUM_IN; ++iind){
      implicit_solver_.fwdSeed(iind,dir).set(fwdSeed(iind,dir));
    }
  }

  // Pass solution guess (if this is the first integration or if hotstart is disabled)
  if(hotstart_==false || integrated_once_==false){
    vector<double>& v = implicit_solver_.output().data();
      
    // Check if an integrator for the startup trajectory has been supplied
    bool has_startup_integrator = !startup_integrator_.isNull();
    
    // Use supplied integrator, if any
    if(has_startup_integrator){
      for(int iind=0; iind<INTEGRATOR_NUM_IN; ++iind){
        startup_integrator_.input(iind).set(input(iind));
      }
      
      // Reset the integrator
      startup_integrator_.reset();
    }
      
    // Integrate, stopping at all time points
    int offs=0;
    for(int k=0; k<coll_time_.size(); ++k){
      for(int j=0; j<coll_time_[k].size(); ++j){
        
        if(has_startup_integrator){
          // Integrate to the time point
          startup_integrator_.integrate(coll_time_[k][j]);
        }
          
        // Save the differential states
        const DMatrix& x = has_startup_integrator ? startup_integrator_.output(INTEGRATOR_XF) : input(INTEGRATOR_X0);
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
      
    // Print
    if(has_startup_integrator && verbose()){
      cout << "startup trajectory generated, statistics:" << endl;
      startup_integrator_.printStats();
    }
  }
    
  // Solve the system of equations
  implicit_solver_.evaluate(nsens);
  
  // Mark the system integrated at least once
  integrated_once_ = true;
}

void CollocationIntegratorInternal::resetB(){
}

void CollocationIntegratorInternal::integrate(double t_out){
  // We will interpolate in the existing grid
  casadi_assert_message(t_out>=t0_ && t_out<=tf_,"CollocationIntegratorInternal::integrate(" << t_out << ") called, but supplied time is not within [" << t0_ << "," << tf_ << "].");
  
  double rt_out = t0_ + (tf_-t_out); 

  // Find out which collocation interval to use
  int k; for (k=0; t_out>coll_time_[k].back();++k) { };
  
  // Find out which collocation interval to use for revsersed system
  int rk; for (rk=0; rt_out>coll_time_[rk].back();++rk) { };
  
  // Local time tau [0,1] on collocation interval k
  double tau = (t_out-coll_time_[k].front()) / (coll_time_[k].back() - coll_time_[k].front());
  
  // Local time tau [0,1] on collocation interval rk for revsersed system
  double rtau = (rt_out-coll_time_[rk].front()) / (coll_time_[rk].back() - coll_time_[rk].front());

  casadi_assert(tau>=0 &&   tau<=1);
  casadi_assert(rtau>=0 &&   rtau<=1);
  
  // Set interpolation function tau input
  Lfcn_.setInput(tau,0);
  Lfcn_.setInput(rtau,1);
  
  int kmul = (nx_+nrx_)*(deg_+1)+(nz_+nrz_)*deg_;
  
  // Pass the saved states of collocation interval k to the interpolating function
  std::copy(implicit_solver_.output(0).data().begin()+k*kmul,implicit_solver_.output(0).data().begin()+(k+1)*kmul,Lfcn_.input(2).data().begin());
  std::copy(implicit_solver_.output(0).data().begin()+rk*kmul,implicit_solver_.output(0).data().begin()+(rk+1)*kmul,Lfcn_.input(3).data().begin());
  std::copy(implicit_solver_.output(1+INTEGRATOR_NUM_IN+0).data().begin()+k*nq_*(deg_+1),implicit_solver_.output(1+INTEGRATOR_NUM_IN+0).data().begin()+(k+1)*nq_*(deg_+1),Lfcn_.input(4).data().begin());
  std::copy(implicit_solver_.output(1+INTEGRATOR_NUM_IN+1).data().begin()+k*nrq_*(deg_+1),implicit_solver_.output(1+INTEGRATOR_NUM_IN+1).data().begin()+(k+1)*nrq_*(deg_+1),Lfcn_.input(5).data().begin()); 
  
  // Pass the forward sensitivities of collocation interval k to the  interpolating function
  for(int dir=0; dir<nsens_; ++dir){
    std::copy(implicit_solver_.fwdSens(0,dir).data().begin()+k*kmul,implicit_solver_.fwdSens(0,dir).data().begin()+(k+1)*kmul,Lfcn_.fwdSeed(2).data().begin());
    std::copy(implicit_solver_.fwdSens(0,dir).data().begin()+rk*kmul,implicit_solver_.fwdSens(0,dir).data().begin()+(rk+1)*kmul,Lfcn_.fwdSeed(3).data().begin());
    std::copy(implicit_solver_.fwdSens(1+INTEGRATOR_NUM_IN+0,dir).data().begin()+k*nq_*(deg_+1),implicit_solver_.fwdSens(1+INTEGRATOR_NUM_IN+0,dir).data().begin()+(k+1)*nq_*(deg_+1),Lfcn_.fwdSeed(4,dir).data().begin());
    std::copy(implicit_solver_.fwdSens(1+INTEGRATOR_NUM_IN+1,dir).data().begin()+k*nrq_*(deg_+1),implicit_solver_.fwdSens(1+INTEGRATOR_NUM_IN+1,dir).data().begin()+(k+1)*nrq_*(deg_+1),Lfcn_.fwdSeed(5,dir).data().begin()); 
  }

  Lfcn_.evaluate(nsens_);
  
  for(int oind=0; oind<INTEGRATOR_NUM_OUT; ++oind){
    output(oind).set(Lfcn_.output(oind));
    for(int dir=0; dir<nsens_; ++dir){
      fwdSens(oind,dir).set(Lfcn_.fwdSens(oind,dir));
    }
  }

}

void CollocationIntegratorInternal::integrateB(double t_out){
  casadi_assert_message(t_out>=t0_ && t_out<=tf_,"CollocationIntegratorInternal::integrateB(" << t_out << ") called, but supplied time is not within [" << t0_ << "," << tf_ <<"].");
}

} // namespace CasADi
