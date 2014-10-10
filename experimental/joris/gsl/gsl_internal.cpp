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


#include "gsl_internal.hpp"
#include "core/function/sx_function_internal.hpp"
#include "core/std_vector_tools.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/function/linear_solver_internal.hpp"
#include "core/function/mx_function.hpp"
#include "core/function/jacobian.hpp"

using namespace std;
namespace casadi {
namespace GSL {

GslInternal* GslInternal::clone() const{
  // Return a deep copy
  Function f = deepcopy(f_);
  Function q = deepcopy(q_);
  GslInternal* node = new GslInternal(f,q);
  node->setOption(dictionary());
  node->jac_f_ = deepcopy(jac_f_);
  node->dt_f_ = deepcopy(dt_f_);
  if(!node->is_init)
    node->init();
  return node;
}
  
GslInternal::GslInternal(const Function& f, const Function& q) : IntegratorInternal(f,q){
  is_init = false;
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "reset", true);
  std::cout << "Warning: GslIntegrator is highly experimental" << std::endl;
}

GslInternal::~GslInternal(){ 
  gsl_odeiv_evolve_free (evolve_ptr);
  gsl_odeiv_control_free (control_ptr);
  gsl_odeiv_step_free (step_ptr);
}


void GslInternal::resetAdj() { return;}

void GslInternal::integrateAdj(double t_out) {return;}

Function GslInternal::getJacobian() {return Function();}
  
LinearSolver GslInternal::getLinearSolver() { return LinearSolver();}

void GslInternal::setLinearSolver(const LinearSolver& linsol, const Function& jac) {
  return;
}
  
void GslInternal::init(){
  // Init ODE rhs function and quadrature functions
  f_.init();
  casadi_assert(f_.getNumInputs()==DAE_NUM_IN);
  casadi_assert(f_.getNumOutputs()==DAE_NUM_OUT);
  if(!q_.isNull()){
    q_.init();
    casadi_assert(q_.getNumInputs()==DAE_NUM_IN);
    casadi_assert(q_.getNumOutputs()==DAE_NUM_OUT);
  }

  // Number of states
  int nx = f_.output(INTEGRATOR_XF).numel();

  // Add quadratures, if any
  if(!q_.isNull()) nx += q_.output().numel();

  // Number of parameters
  int np = f_.input(DAE_P).numel();

  setDimensions(nx,np);

  // If time was not specified, initialise it.
  if (f_.input(DAE_T).numel()==0) {
    std::vector<MX> in1(DAE_NUM_IN);
    in1[DAE_T] = MX("T");
    in1[DAE_Y] = MX("Y",f_.input(DAE_Y).size1(),f_.input(DAE_Y).size2());
    in1[DAE_YDOT] = MX("YDOT",f_.input(DAE_YDOT).size1(),f_.input(DAE_YDOT).size2());
    in1[DAE_P] = MX("P",f_.input(DAE_P).size1(),f_.input(DAE_P).size2());
    std::vector<MX> in2(in1);
    in2[DAE_T] = MX();
    f_ = MXFunction(in1,f_.call(in2));
    f_.init();
  }
  
  // We only allow for 0-D time
  casadi_assert_message(f_.input(DAE_T).numel()==1, "IntegratorInternal: time must be zero-dimensional, not (" <<  f_.input(DAE_T).size1() << 'x' << f_.input(DAE_T).size2() << ")");
  
  // ODE right hand side must be a dense matrix
  casadi_assert_message(f_.output(DAE_RES).dense(),"ODE right hand side must be dense: reformulate the problem");
  
  // States and RHS should match 
  casadi_assert_message(f_.output(DAE_RES).size()==f_.input(DAE_Y).size(),
    "IntegratorInternal: rhs of ODE is (" <<  f_.output(DAE_RES).size1() << 'x' << f_.output(DAE_RES).size2() << ") - " << f_.output(DAE_RES).size() << " non-zeros" << std::endl <<
    "              ODE state matrix is (" <<  f_.input(DAE_Y).size1() << 'x' << f_.input(DAE_Y).size2() << ") - " << f_.input(DAE_Y).size() << " non-zeros" << std::endl <<
    "Mismatch between number of non-zeros"
  );

  IntegratorInternal::init();
  
  jac_f_ = Jacobian(f_,DAE_Y,DAE_RES);
  dt_f_ = Jacobian(f_,DAE_T,DAE_RES);
  
  jac_f_.init();
  dt_f_.init();
  
  // define the type of routine for making steps: 
  type_ptr = gsl_odeiv_step_rkf45;
  // some other possibilities (see GSL manual):          
  //   = gsl_odeiv_step_rk4;
  //   = gsl_odeiv_step_rkck;
  //   = gsl_odeiv_step_rk8pd;
  //   = gsl_odeiv_step_rk4imp;
  //   = gsl_odeiv_step_bsimp;  
  //   = gsl_odeiv_step_gear1;
  //   = gsl_odeiv_step_gear2;
  
  step_ptr = gsl_odeiv_step_alloc (type_ptr, nx);
  control_ptr = gsl_odeiv_control_y_new (abstol_, reltol_);
  evolve_ptr = gsl_odeiv_evolve_alloc (nx);
  
  my_system.function = rhs_wrapper;	// the right-hand-side functions dy[i]/dt 
  my_system.jacobian = jac_wrapper;	// the Jacobian df[i]/dy[j] 
  my_system.dimension = nx;	// number of diffeq's 
  my_system.params = this;	// parameters to pass to rhs and jacobian
  
  is_init = true;
  
}
  
void GslInternal::reset(int fsens_order, int asens_order){
  if(monitored("reset")){
    cout << "initial state: " << endl;
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
  }

  // Reset timers
  t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;
  
  // Get the time horizon
  t_ = t0_;
  
  int flag = gsl_odeiv_evolve_reset(evolve_ptr);
	if(flag!=GSL_SUCCESS) gsl_error("Reset",flag);
	
  output(INTEGRATOR_XF).set(input(INTEGRATOR_X0));

}

void GslInternal::integrate(double t_out){
  double h = 1e-6;		// starting step size for ode solver 
  std::cout << "I wanna integrate from " << t_ << " to " << t_out << std::endl;
	int flag = gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t_, t_out, &h, &output(INTEGRATOR_XF).data()[0]);
	if(flag!=GSL_SUCCESS) gsl_error("Integrate",flag);
	std::cout << "GSL returned succes: " << flag << std::endl;
}

void GslInternal::printStats(std::ostream &stream) const{

    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    
    int flag;

    stream << "number of steps taken by CVODES: " << nsteps << std::endl;
    stream << "number of calls to the user's f function: " << nfevals << std::endl;
    stream << "number of calls made to the linear solver setup function: " << nlinsetups << std::endl;
    stream << "number of error test failures: " << netfails << std::endl;
    stream << "method order used on the last internal step: " << qlast << std::endl;
    stream << "method order to be used on the next internal step: " << qcur << std::endl;
    stream << "actual value of initial step size: " << hinused << std::endl;
    stream << "step size taken on the last internal step: " << hlast << std::endl;
    stream << "step size to be attempted on the next internal step: " << hcur << std::endl;
    stream << "current internal time reached: " << tcur << std::endl;
    stream << std::endl;

    stream << "Time spent in the ODE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: " << t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: " << t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: " << t_lsetup_fac << " s." << endl;
    stream << std::endl;

}
  
map<int,string> GslInternal::calc_flagmap(){
  map<int,string> f;
  f[GSL_SUCCESS] = "GSL_SUCCESS";
  f[GSL_FAILURE] = "GSL_FAILURE";
  f[GSL_CONTINUE] = "GSL_CONTINUE";
  f[GSL_EDOM] = "GSL_EDOM";
  f[GSL_ERANGE] = "GSL_ERANGE";
  f[GSL_EFAULT] = "GSL_EFAULT";
  f[GSL_EINVAL] = "GSL_EINVAL";
  f[GSL_EFAILED] = "GSL_EFAILED";
  f[GSL_EFACTOR] = "GSL_EFACTOR";
  f[GSL_ESANITY] = "GSL_ESANITY";
  f[GSL_ENOMEM] = "GSL_ENOMEM";
  f[GSL_EBADFUNC] = "GSL_EBADFUNC";
  f[GSL_ERUNAWAY] = "GSL_ERUNAWAY";
  f[GSL_EMAXITER] = "GSL_EMAXITER";
  f[GSL_EZERODIV] = "GSL_EZERODIV";
  f[GSL_EBADTOL] = "GSL_EBADTOL";
  f[GSL_ETOL] = "GSL_ETOL";
  f[GSL_EUNDRFLW] = "GSL_EUNDRFLW";
  f[GSL_EOVRFLW] = "GSL_EOVRFLW";
  f[GSL_ELOSS] = "GSL_ELOSS";
  f[GSL_EROUND] = "GSL_EROUND";
  f[GSL_EBADLEN] = "GSL_EBADLEN";
  f[GSL_ENOTSQR] = "GSL_ENOTSQR";
  f[GSL_ESING] = "GSL_ESING";
  f[GSL_EDIVERGE] = "GSL_EDIVERGE";
  f[GSL_EUNSUP] = "GSL_EUNSUP";
  f[GSL_EUNIMPL] = "GSL_EUNIMPL";
  f[GSL_ECACHE] = "GSL_ECACHE";
  f[GSL_ETABLE] = "GSL_ETABLE";
  f[GSL_ENOPROG] = "GSL_ENOPROG";
  f[GSL_ENOPROGJ] = "GSL_ENOPROGJ";
  f[GSL_ETOLF] = "GSL_ETOLF";
  f[GSL_ETOLX] = "GSL_ETOLX";
  f[GSL_ETOLG] = "GSL_ETOLG";
  f[GSL_EOF] = "GSL_EOF";
  return f;
}
  
map<int,string> GslInternal::flagmap = GslInternal::calc_flagmap();

void GslInternal::gsl_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult GSL documentation.";
  casadi_error(ss.str());
}
  
  
int GslInternal::rhs_wrapper(double t, const double y[], double f[], void *userdata) {
  try{
    casadi_assert(userdata);
    GslInternal *this_ = (GslInternal*)userdata;
    this_->rhs(t,y,f);
    return 0;
  } catch(exception& e){
    cerr << "fun failed: " << e.what() << endl;;
    return 1;
  }
}

void GslInternal::rhs(double t, const double y[], double f[]) {

  // Get time
  time1 = clock();
  
  f_.setInput(t,"t");
  f_.setInput(y,"y");
  f_.setInput(f_.input(DAE_P),"p");
  
  f_.evaluate();
  
  std::cout << "rhs @ " << t << ": " << f_.input(DAE_Y) << " -> " << f_.output() << std::endl;
  
  std::copy(f_.output().data().begin(),f_.output().data().end(),f);
  
  // Log time duration
  time2 = clock();
  t_res += double(time2-time1)/CLOCKS_PER_SEC;
}
  
int GslInternal::jac_wrapper(double t, const double y[], double *dfdy, double dfdt[], void *userdata) {
  try{
    casadi_assert(userdata);
    GslInternal *this_ = (GslInternal*)userdata;
    this_->jac(t,y,dfdy,dfdt);
    return 0;
  } catch(exception& e){
    cerr << "jac failed: " << e.what() << endl;;
    return 1;
  }
}

void GslInternal::jac(double t, const double y[], double *dfdy, double dfdt[]){
  // Get time
  time1 = clock();

  // Pass inputs to the jacobian function
  jac_f_.setInput(t,"t");
  jac_f_.setInput(y,"y");
  jac_f_.setInput(f_.input(DAE_P),"p");

  // Pass inputs to the jacobian function
  dt_f_.setInput(t,"t");
  dt_f_.setInput(y,"y");
  dt_f_.setInput(f_.input(DAE_P),"p");
  
  // Evaluate
  jac_f_.evaluate();
  
  // Evaluate
  dt_f_.evaluate();
  
  std::copy(jac_f_.output().data().begin(),jac_f_.output().data().end(),dfdy);
  std::copy(dt_f_.output().data().begin(),dt_f_.output().data().end(),dfdt);
  
  std::cout << "jac @ " << t << ": " << jac_f_.output() << " , " << dt_f_.output() << std::endl;
    
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

void GslInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
}


} // namespace GSL
} // namespace casadi

