/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "cvodes_integrator.hpp"

#include <cvodes/cvodes.h>            /* prototypes for CVODE fcts. and consts. */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h> 
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>

#include <limits>
#include <iterator>
#include <vector>
#include <cassert>

#include "../../casadi/stl_vector_tools.hpp"

#define ASENS 0

using namespace std;

namespace OPTICON{

CvodesIntegratorNode::CvodesIntegratorNode(OCP_old &ocp_) : SundialsIntegrator(ocp_){
  addOption("linear_multistep_method",   OT_STRING,  "bdf"); // "bdf" or "adams"
  addOption("nonlinear_solver_iteration",OT_STRING,  "newton"); // "newton" or "functional"
  addOption("interpolation_type",        OT_STRING,  "polynomial"); // "polynomial" or "hermite"
  addOption("max_num_steps",             OT_INTEGER, 10000); // "maximum number of steps"
  xx = 0;
  quad = 0;
}

CvodesIntegratorNode* CvodesIntegratorNode::create(OCP_old& ocp) const{
  return new CvodesIntegratorNode(ocp);
}

IntegratorNode_old* CvodesIntegratorNode::creator(OCP_old& ocp){
  return new CvodesIntegratorNode(ocp);
}

CvodesIntegratorNode::~CvodesIntegratorNode(){
  for(int i=0; i<cvode_mem.size(); ++i) CVodeFree(&cvode_mem[i]);

  for(int i=0; i<fsens.size(); ++i) N_VDestroy(fsens[i]);

#if ASENS
  for(int i=0; i<lambda.size(); ++i) N_VDestroy(lambda[i]);
  for(int i=0; i<bquad.size(); ++i) N_VDestroy(bquad[i]);
#endif
  if(quad) N_VDestroy(quad);
  if(xx) N_VDestroy(xx);
}

void CvodesIntegratorNode::init(){
  SundialsIntegrator::init();

  // Collect all the quadratures
  int n_quad = 0;
  for(int ord=0; ord<3; ++ord){
    for(int i=0; i<ocp.numOut(ord); ++i){
      int oind = output_ind[ord][i];
      quad_ind[ord].push_back(n_quad);
      if(ocp.output(oind).has_q) n_quad += ocp.function(OUT_G,ocp.output(oind).g_ind).fcn.size();
      }
      quad_ind[ord].push_back(n_quad);
    }

  // Collect all the forward sensitivities
  int n_fsens = 0;
  for(int ord=1; ord<3; ++ord){ // from order 1!    
    // indices of the jacobians and the initial conditions
    for(int i=0; i<ocp.numOut(ord); ++i){
      int oind = output_ind[ord][i];
      fsens_ind[ord].push_back(n_fsens);
      n_fsens += ocp.function(OUT_D,ocp.output(oind).d_ind).fcn.size2();
    }
  }

  // Request a jacobian of the de right hand side
  ocp.function(DE_F).addJac(ocp.getVars(VAR_X));
  ocp.function(DE_F).addJac(ocp.getVars(VAR_DEP));


  jacf_ind.clear();
  jacic_ind.clear();
  jacoutg_ind.clear();
  for(int ord=1; ord<3; ++ord){ // from order 1!    
    // indices of the jacobians and the initial conditions
    for(int i=0; i<output_ind[ord].size(); ++i){
      int oind = output_ind[ord][i];
      if(ocp.function(OUT_D,ocp.output(oind).d_ind).fcn.size2() == 0) throw "error, no dir";

      // Get all directions
      const Matrix& dirs = ocp.function(OUT_D,ocp.output(oind).d_ind).fcn;

      // Get parameters
      const Matrix& p = ocp.output(oind).p;

      ocp.function(OUT_G,oind).addJac(ocp.getVars(VAR_X));

      // loop over directions
      for(int j=0; j<dirs.size2(); ++j){
        Matrix d;
        getColumn(d,dirs,j);
  
        // Request jacobians
        jacf_ind.push_back(ocp.function(DE_F).addJac(p,d));
        jacic_ind.push_back(ocp.function(IC_EXPFUN).addJac(p,d));
        jacoutg_ind.push_back(ocp.function(OUT_G,output_ind[ord][i]).addJac(p,d));
        }
      }
      fsens_ind[ord].push_back(n_fsens);
    }

  // Initialize the ocp
  ocp.init();
  
  // Allocate memory for the state vector
  xx = N_VNew_Serial(ocp.numVar(VAR_X));

  // Allocate N-Vectors quadrature functions
  if(quad) N_VDestroy(quad);
  if(n_quad>0) quad = N_VNew_Serial(n_quad);

  // Allocate N_Vectors for the forward sensitivities
  fsens.resize(n_fsens);
  for(int i=0; i<n_fsens; ++i)
    fsens[i] =  N_VNew_Serial(NV_LENGTH_S(xx));

#if ASENS
  // Allocate N_Vectors for the adjoints (better: split up quad/no quad)
  lambda.resize(ocp.numFcn(ASENS_FUN) + ocp.numFcn(QSENS_FUN));
  for(int i=0; i<lambda.size(); ++i)
    lambda[i] = N_VNew_Serial(NV_LENGTH_S(xx));

  // Allocate N_vectors for adjoint sensitivities
  asens.resize(ocp.numFcn(ASENS_FUN) + ocp.numFcn(QSENS_FUN)); // split up!!
  for(int i=0; i<asens.size(); ++i){
    if(i<ocp.numFcn(ASENS_FUN)) asens[i] = N_VNew_Serial(ocp.numAdj(i));
    else               asens[i] = N_VNew_Serial(ocp.numAdjQuad(i));
  }

  // Allocate N_Vectors for the backward quadrature
  bquad.resize(lambda.size());
  for(int i=0; i<bquad.size(); ++i)
    bquad[i] = N_VNew_Serial(ocp.numAdjQuad(i));
#endif

  // Create and initialize cvode memory block(s)
  if(!output_ind[1].empty() || !output_ind[1].empty()){
    // Initialize a cvodes memory block for each substage
    cvode_mem.resize(ocp.nsub());
    adjointID.resize(ocp.nsub());
    adjointUD.resize(ocp.nsub());
    for(int i=0; i<ocp.nsub(); ++i)
      initCvodes(i);
  } else {
    // Just initialize one memory block
    cvode_mem.resize(1);
    adjointID.resize(1);
    adjointUD.resize(1);
    initCvodes(0);
  }
  
  // Set the outputs
  output.resize(ocp.numOut());
  for(int i=0; i<output.size(); ++i){
     // Reference to the OCPOutput
     const OCPOutput& oo = ocp.output(i);

     // Get the size
     MatrixSize sz(oo.nf, oo.nt);

     // Get the derivative order
     int nder = oo.ord;

     output[i] = FunctionIO(sz,nder);
  }
}

void CvodesIntegratorNode::initCvodes(int k){

  int lmm; // linear multistep method
  if(checkOption("linear_multistep_method","adams"))  lmm = CV_ADAMS;
  else if(checkOption("linear_multistep_method","bdf")) lmm = CV_BDF;
  else throw "Unknown linear multistep method";

  int iter; // nonlinear solver iteration
  if(checkOption("nonlinear_solver_iteration","newton")) iter = CV_NEWTON;
  else if(checkOption("nonlinear_solver_iteration","functional")) iter = CV_FUNCTIONAL;
  else throw "Unknown nonlinear solver iteration";

  // allocate a cvodes solver
  cvode_mem[k] = CVodeCreate(lmm, iter); 
  if(cvode_mem[k] == 0) throw "Error in CVodeCreate";

  // Set the error handling function
  sundialsAssert(CVodeSetErrHandlerFn(cvode_mem[k], ehfun_static, this));

  // initialize cvodes at t=k
  sundialsAssert(CVodeInit(cvode_mem[k], &f_static, k, xx));

  // set tolerances 
  sundialsAssert(CVodeSStolerances(cvode_mem[k], reltol, abstol));
  
  // attach linear solver
  if(checkOption("linear_solver","dense")){
    // Dense jacobian
    sundialsAssert(CVDense(cvode_mem[k], NV_LENGTH_S(xx)));
    if(exact_jacobian) sundialsAssert(CVDlsSetDenseJacFn(cvode_mem[k], djac_static));
  } else if(checkOption("linear_solver","band")) {
    // Banded jacobian
    sundialsAssert(CVBand(cvode_mem[k], NV_LENGTH_S(xx), int(getOption("mupper")), int(getOption("mlower"))));
    if(exact_jacobian) sundialsAssert(CVDlsSetBandJacFn(cvode_mem[k], bjac_static));  
  } else if(checkOption("linear_solver","sparse")) {
    // Sparse jacobian
    // Preconditioning type
    int pretype;
    if(checkOption("pretype","none"))               pretype = PREC_NONE;
    else if(checkOption("pretype","left"))          pretype = PREC_LEFT;
    else if(checkOption("pretype","right"))         pretype = PREC_RIGHT;
    else if(checkOption("pretype","both"))          pretype = PREC_BOTH;
    else                                            throw "Unknown preconditioning type";

    // Attach the sparse solver  
    if(checkOption("sparse_solver","spgmr"))        sundialsAssert(CVSpgmr(cvode_mem[k], pretype, maxl));
    else if(checkOption("sparse_solver","spbcg"))   sundialsAssert(CVSpbcg(cvode_mem[k], pretype, maxl));
    else if(checkOption("sparse_solver","sptfqmr")) sundialsAssert(CVSptfqmr(cvode_mem[k], pretype, maxl));
    else                                            throw "Unknown sparse solver";
       
    // Attach functions for jacobian information
    if(exact_jacobian) sundialsAssert(CVSpilsSetJacTimesVecFn(cvode_mem[k], jtimes_static));
  } else throw "Unknown linear solver ";

  // set user data
  sundialsAssert(CVodeSetUserData(cvode_mem[k], this));

  // Initialize the quadrature
  if(quad){

    // Pass the quadrature function pointer
    sundialsAssert(CVodeQuadInit(cvode_mem[k], fQ_static, quad));  

    // Make sure that the quadrature calculation is used when determining the step sizes
    bool errconQ = true;
    sundialsAssert(CVodeSetQuadErrCon(cvode_mem[k], errconQ));

    // Set error tolerances
    double reltolQ = reltol;
    double abstolQ = abstol;
    sundialsAssert(CVodeQuadSStolerances(cvode_mem[k], reltolQ, abstolQ)); 
  }

  // Initialize forward sensitivity analysis
  if(!fsens.empty()){
    int Ns = fsens.size(); // number of sensitivities
  
    // method
    int ism = CV_SIMULTANEOUS;
    //  int ism = CV_STAGGERED;
    //  int ism = CV_STAGGERED1;
    sundialsAssert(CVodeSensInit1(cvode_mem[k], Ns, ism, rhs_fsens_static, &fsens[0]));
  
    // Set tolerances
    double reltolS = reltol;
    double abstolS[Ns];
    for(int i=0; i<Ns; ++i) abstolS[i] = abstol;
    sundialsAssert(CVodeSensSStolerances(cvode_mem[k], reltolS, &abstolS[0]));
  }

#if ASENS

  if(ocp.numFcn(ASENS_FUN) + ocp.numFcn(QSENS_FUN) > 0){

  int interpType;
  if(checkOption("interpolation_type","polynomial"))    interpType = CV_POLYNOMIAL;
  else if(checkOption("interpolation_type","hermite"))  interpType = CV_HERMITE;
  else throw "Unknown interpolation type";

  // Initialize the forward integration with checkpointing
  long Nd = 20; // option
  sundialsAssert(CVodeAdjInit(cvode_mem[k], Nd, CV_POLYNOMIAL));

  // Allocate id for adjoints
  adjointID[k].resize(lambda.size());

  // Allocate userdata structures
  adjointUD[k].resize(lambda.size());

  for(int i=0; i<lambda.size(); ++i){
    // Create a memory block for backward integration
    int lmmB = lmm;
    int iterB = iter;
    sundialsAssert(CVodeCreateB(cvode_mem[k], lmmB, iterB, &adjointID[k][i]));  

    // Allocate memory for backward integration
    sundialsAssert(CVodeInitB(cvode_mem[k], adjointID[k][i], rhs_asens_static, k, lambda[i]));

    // Set tolerances
    double reltolB = reltol;
    double abstolB = abstol;
    sundialsAssert(CVodeSStolerancesB(cvode_mem[k], adjointID[k][i], reltolB , abstolB ));

    adjointUD[k][i].this_ = this;
    adjointUD[k][i].id = i;

    // Pass the user data structure
    sundialsAssert(CVodeSetUserDataB(cvode_mem[k], adjointID[k][i], &adjointUD[k][i]));

  // Attach a linear solver 
  if(checkOption("linear_solver","dense")){
    // Dense jacobian
    sundialsAssert(CVDenseB(cvode_mem[k], adjointID[k][i], NV_LENGTH_S(xx)));
    // if(exact_jacobian) sundialsAssert(CVDlsSetDenseJacFnB(cvode_mem[k], adjointID[k][i], djac_static));
  } else if(checkOption("linear_solver","band")) {
    // Banded jacobian
    sundialsAssert(CVBandB(cvode_mem[k], adjointID[k][i], NV_LENGTH_S(xx), int(getOption("mupper")), int(getOption("mlower"))));
    // if(exact_jacobian) sundialsAssert(CVDlsSetBandJacFnB(cvode_mem[k], adjointID[k][i], bjac_static));  
  } else if(checkOption("linear_solver","sparse")) {
    // Sparse jacobian
    // Preconditioning type
    int pretype;
    if(checkOption("pretype","none"))               pretype = PREC_NONE;
    else if(checkOption("pretype","left"))          pretype = PREC_LEFT;
    else if(checkOption("pretype","right"))         pretype = PREC_RIGHT;
    else if(checkOption("pretype","both"))          pretype = PREC_BOTH;
    else                                            throw "Unknown preconditioning type";

    // Attach the sparse solver  
    if(checkOption("sparse_solver","spgmr"))        sundialsAssert(CVSpgmrB(cvode_mem[k], adjointID[k][i], pretype, maxl));
    else if(checkOption("sparse_solver","spbcg"))   sundialsAssert(CVSpbcgB(cvode_mem[k], adjointID[k][i], pretype, maxl));
    else if(checkOption("sparse_solver","sptfqmr")) sundialsAssert(CVSptfqmrB(cvode_mem[k], adjointID[k][i], pretype, maxl));
    else                                            throw "Unknown sparse solver";
       
    // Attach functions for jacobian information
    // if(exact_jacobian) sundialsAssert(CVSpilsSetJacTimesVecFnB(cvode_mem[k], adjointID[k][i], jtimes_static));
  } else throw "Unknown linear solver ";

    // Initialize the backward quadratures
    sundialsAssert(CVodeQuadInitB(cvode_mem[k],adjointID[k][i], fQB_static, bquad[i]));    
  }
}
#endif
  sundialsAssert(CVodeSetMaxStep(cvode_mem[k], 0.1));
  sundialsAssert(CVodeSetMaxNumSteps(cvode_mem[k], int(getOption("max_num_steps"))));  
}

std::ostream& operator<<(std::ostream &stream, const vector<long>& v){
    copy(v.begin(), v.end(), ostream_iterator<long>(stream, " "));
    return stream;
 }

void CvodesIntegratorNode::printStats(std::ostream &stream) const{

    int nm = cvode_mem.size();
    vector<long> nsteps(nm), nfevals(nm), nlinsetups(nm), netfails(nm);
    vector<int> qlast(nm), qcur(nm);
    vector<double> hinused(nm), hlast(nm), hcur(nm), tcur(nm);

  for(int k=0; k<nm; ++k){
    sundialsAssert(CVodeGetIntegratorStats(cvode_mem[k], &nsteps[k], &nfevals[k],
			  &nlinsetups[k], &netfails[k], &qlast[k], &qcur[k], 	 
			  &hinused[k], &hlast[k], &hcur[k], &tcur[k]));
  }

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


#if 0
  // Quadrature
   if(ops.quadrature && ocp.hasFunction(LTERM)){
      long nfQevals, nQetfails;
      flag = CVodeGetQuadStats(cvode_mem[k], &nfQevals, &nQetfails);  
      if(flag != CV_SUCCESS) throw "Error in CVodeGetQuadStats";

      stream << "Quadrature: " << std::endl;
      stream << "number of calls made to the user's quadrature right-hand side function: " << nfQevals << std::endl;
      stream << "number of local error test failures due to quadrature variables: " <<  nQetfails << std::endl;
      stream << std::endl;
 }
#endif
}


void CvodesIntegratorNode::ehfun_static(int error_code, const char *module, const char *function, char *msg, void *user_data){
  switch(error_code){
    case CV_MEM_NULL:
        std::cerr << "CV_MEM_NULL"  << std::endl; break;
    case CV_NO_ADJ:
	std::cerr << "CV_NO_ADJ" << std::endl; break;
    case CV_NO_FWD:
	std::cerr << "CV_NO_FWD" << std::endl; break;
    case CV_NO_BCK:
	std::cerr << "CV_NO_BCK" << std::endl; break;
    case CV_ILL_INPUT:
	std::cerr << "CV_ILL_INPUT" << std::endl; break;
    case CV_TOO_MUCH_WORK:
	std::cerr << "CV_TOO_MUCH_WORK" << std::endl; break;
    case CV_TOO_MUCH_ACC:
	std::cerr << "CV_TOO_MUCH_ACC" << std::endl; break;
    case CV_ERR_FAILURE:
	std::cerr << "CV_ERR_FAILURE" << std::endl; break;
    case CV_CONV_FAILURE:
	std::cerr << "CV_CONV_FAILURE" << std::endl; break;
    case CV_LSETUP_FAIL:
	std::cerr << "CV_LSETUP_FAIL" << std::endl; break;
    case CV_LSOLVE_FAIL:
	std::cerr << "CV_LSOLVE_FAIL" << std::endl; break;
    case CV_REIFWD_FAIL:
	std::cerr << "CV_REIFWD_FAIL" << std::endl; break;
    case CV_FWD_FAIL:
	std::cerr << "CV_FWD_FAIL" << std::endl; break;
    default:
        std::cerr << "Unknown error"  << std::endl;
        CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)user_data;
        this_->printStats();
        N_VPrint_Serial(this_->xx);
  }

  std::cerr << "Error code: " << error_code << std::endl;
  std::cerr << "Module: " << module << std::endl;
  std::cerr << "Function: " << function << std::endl;
}

int CvodesIntegratorNode::f_static(realtype t, N_Vector x, N_Vector xdot, void *user_data){
  try{
    CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)user_data;
    OCP_old& ocp = this_->ocp;

    // Pass the time to the ocp
    this_->ocp.setTime(t);

    // Pass the state vector to the ocp
    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(x));
    ocp.function(DE_F).eval(NV_DATA_S(xdot));

    // Scale time
    for(int i=0; i<NV_LENGTH_S(x); ++i)
      NV_DATA_S(xdot)[i] *= ocp.getDur();

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
}

int CvodesIntegratorNode::djac_static(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
try{
  CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)user_data;
  OCP_old &ocp = this_->ocp;

    // Retrieve the structure that holds the Jacobian
    _DlsMat *Jac_ = (_DlsMat *)Jac;

    // Get a pointer to the first element of the of the jacobian
    double *jdata = Jac->data;

    // Retrieve the leading dimension of the jacobian
    int ld_jac = Jac->ldim;


  // Pass the time to the ocp
  ocp.setTime(t);

  // Pass the state vector to the ocp
  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

  double temp[N][N];
  // Evaluate the jacobian
  ocp.function(DE_F).eval_dense(&temp[0][0], ld_jac, 1);
  for(int i=0; i<N; ++i)
    for(int j=0; j<N; ++j)
      jdata[j+i*ld_jac] = temp[j][i];

  // Scale time
  for(int i=0; i<N; ++i)
    for(int j=0; j<N; ++j)
      jdata[j+i*ld_jac] *= ocp.getDur();

  return 0;
} catch (const char * str){
  std::cerr << str << std::endl;
  return -1;
}
}

int CvodesIntegratorNode::jtimes_static( N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp){
try{
  CvodesIntegratorNode *this_ = (CvodesIntegratorNode*)user_data;
  OCP_old& ocp = this_->ocp;

  // Pass the time to the ocp
  ocp.setTime(t);

  // Pass the state vector to the ocp
  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

  // Evaluate
  ocp.function(DE_F).eval_times_vector(NV_DATA_S(Jv), NV_DATA_S(v), 1);

  // Scale time
  for(int i=0; i<NV_LENGTH_S(y); ++i)
    NV_DATA_S(Jv)[i] *= ocp.getDur();

  return 0;
} catch (const char * str){
  std::cerr << str << std::endl;
  return -1;
}
}

int CvodesIntegratorNode::g_static( realtype t, N_Vector y, realtype *gout, void *user_data){
  try{
    CvodesIntegratorNode *this_ = (CvodesIntegratorNode*)user_data;
    OCP_old& ocp = this_->ocp;

    // Pass the time to the ocp
    ocp.setTime(t);

    // Pass the state vector to the ocp
    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

    // Evaluate the root finding functions
    ocp.function(SWITCH).eval(gout);

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
}

int CvodesIntegratorNode::bjac_static(int N, int mupper, int mlower, 	 
 		 	realtype t, N_Vector y, N_Vector fy, 	 
 		 	DlsMat Jac, void *user_data, 	 
 		 	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  try{
    // Retrieve the pointer to the object and the ocp
    CvodesIntegratorNode *this_ = (CvodesIntegratorNode*)user_data;
    OCP_old &ocp = this_->ocp;

    // Retrieve the structure that holds the Jacobian
    _DlsMat *Jac_ = (_DlsMat *)Jac;

    // Get a pointer to the first element of the of the jacobian (disregarding the extra subdiagonals allocated for lapack factorizations)
    double *jdata = Jac->data + Jac->s_mu - Jac->mu;

    // Retrieve the leading dimension of the banded jacobian
    int ld_jac = Jac->ldim;

    // Pass the time to the ocp
    ocp.setTime(t);

    // Pass the state vector to the ocp
    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

    // Evaluate the jacobian
    ocp.function(DE_F).eval_band(jdata, mlower, mupper, ld_jac, 1);

    // Scale time
    for(int i=0; i<NV_LENGTH_S(y); ++i)
      for(int j=0; j<mlower+1+mupper; ++j)
	jdata[j+i*ld_jac] *= ocp.getDur();

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
}


void CvodesIntegratorNode::evaluate(int tape_order){

  setOption("use_events",  false);
  setOption("use_output",  true);
  setOption("use_fsens",   true);
  setOption("use_asens",   true);
  setOption("use_quad",    true);

  // Integrate forward once with root finding
  evaluateFwd();

  // Save the time points to the OCP
  // integrator.saveTimePoints();

  // Integrate forward again without root finding but with sensitivities
//  integrateForward(false, true, true, true, quadrature);
  
  // Integrate backward to get adjoints
  evaluateAdj();

}

void CvodesIntegratorNode::getInitialCondition(){

    // Start at the beginning of the first substage
    ocp.setTime(0,0);
  
    // Get a pointer to xx (current value of the state vector)
    double *x_data = NV_DATA_S(xx);

    // Explicitly defined initial conditions
    vector<int> ic_exp_ind = ocp.getVarInds(ocp.function(IC_EXPVAR));

    // Assert lengths
    if(ic_exp_ind.size() != NV_LENGTH_S(xx)) throw "Wrong number of initial conditions";

    // Evaluate
    double temp[ic_exp_ind.size()];
    ocp.function(IC_EXPFUN).eval(&temp[0]);


    // Pass state to the ocp
    ocp.setArg(ic_exp_ind, &temp[0]);

    // This shouldn't be necessary
    ocp.getArg(ocp.findVars(VAR_X),x_data);
}


void CvodesIntegratorNode::evaluateFwd(bool use_tape){
  bool use_events = getOption("use_events");
  bool use_fsens  = getOption("use_fsens");
  bool use_asens  = getOption("use_asens");

  if(use_events && use_asens) throw "The combination events and adjoint sensitivities is not possible!";

  // Modify the options
  if(!ocp.hasFcn(SWITCH)) use_events = false; // no switching functions ==> no events
  if(fsens.empty()) use_fsens = false; // no forward sensitivities

#if ASENS
  if(lambda.empty()) use_asens = false; // no adjoints
#else
  use_asens = false;
#endif

  int flag;

  // Evaluate time points
  ocp.evalTimepoints();

  // Prepare for saving outputs
  for(int i=0; i<i_out.size(); ++i){
      i_out[i] = 0;
      ocp.function(OUT_T,i).eval(&ocp.output(i).tdata[0]);
  }

  // Get consistent initial conditions
  getInitialCondition();

  // Get the initial conditions for the forward sensitivity equations
  if(use_fsens)
    for(int i=0; i<fsens.size(); ++i){
        N_VConst(0, fsens[i]);
//	ocp.eval(IC_IMPFUN,NV_DATA_S(fsens[i]),0, 1, i);
    }

  // Set the quadrature function to zero
  if(quad) N_VConst(0, quad);

  // New parameters with definitions
  tp_new.clear(); tp_new.resize(ocp.nsub());
  tp_def.clear(); tp_def.resize(ocp.nsub());
  tp_val.clear(); tp_val.resize(ocp.nsub());

  // Start with the first cvode memory block
  int mind = 0;

    // Integrate over all substages
    for(int sub=0; sub<ocp.nsub(); ++sub){
  
      double t = sub;
      ocp.setTime(t,sub);

      // stop at the end of the stage
      double t_stop = sub+1;

      // Initialize the solver
      if(sub==0 || use_asens){

	sundialsAssert(CVodeReInit(cvode_mem[mind], t, xx));
	// Initialize or turn off the forward sensititivities
	if(!fsens.empty()){
	  if(use_fsens) sundialsAssert(CVodeSensReInit(cvode_mem[mind], CV_STAGGERED1, &fsens[0]));
	  else          sundialsAssert(CVodeSensToggleOff(cvode_mem[mind]));
	}

	// Initialize or turn off the root finding
	sundialsAssert(CVodeRootInit(cvode_mem[mind], use_events ? ocp.function(SWITCH).size() : 0, g_static));

	// Re-initialize the quadratures finding
	if(quad) sundialsAssert(CVodeQuadReInit(cvode_mem[mind], quad));
      }

      // Set the solver to stop at the end of the substage
      sundialsAssert(CVodeSetStopTime(cvode_mem[mind], t_stop));

      while(1){

	// Integrate to this point (and possibly past)
	double t_next_out = nextOutput();
	double t_next = t_next_out < t_stop ? t_next_out : t_stop;

	// Try to integrate to the next time point
	if( fabs(t_next-t) < reltol){
	  flag = CV_SUCCESS; // "successful integration" if we are already at the desired point
	} else {

	  // Make an event iteration to determine the values of the switches
	  eventIteration(ocp);

	  // Integrate to the next point
	  if(use_asens){
	    int ncheck; // monitor this?
	    flag = CVodeF(cvode_mem[mind], t_next, xx, &t, CV_NORMAL, &ncheck);
	  } else {
	    flag = CVode(cvode_mem[mind], t_next, xx, &t, CV_NORMAL);
	  }

	  // Make sure that we got a successful return
	  if(!(flag == CV_ROOT_RETURN || flag == CV_SUCCESS || flag == CV_TSTOP_RETURN)) throw "CVode failed";

	  // Print progress
// 	  std::cout << "t = " << t << std::endl;

	  // Pass the new time to the ocp
	  ocp.setTime(t);

	  // Get the sensitivities
	  double tout_sens; // do we need this?
	  if(use_fsens) sundialsAssert(CVodeGetSens(cvode_mem[mind],&tout_sens,&fsens[0]));

	  // Pass the value of the state vector to the ocp
	  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(xx));

	  // Get the quadrature
	  double tout_quad; // do we need this?
	  if(quad) sundialsAssert(CVodeGetQuad(cvode_mem[mind], &tout_quad, quad));
	}

	// Check if we reached an output time
	if(fabs(t_next_out-t) <= reltol) saveOutput();

	// Break the inner loop if we are close enough to the end of the substage
	if(fabs(t_stop-t)<=reltol) break;

	// Check if a time event occured before reaching the output time
	//	if(flag == CV_TSTOP_RETURN) time_event();

	// Check if state event occured before the output time was reached
	if(flag == CV_ROOT_RETURN) state_event(cvode_mem[mind]);
    }

    // Next memory block
    if(use_asens) mind++;
  }




  // Print the result
//   N_VPrint_Serial(quad);

}

void CvodesIntegratorNode::evaluateAdj(){
return;
#if ASENS

  // Terminal conditions of the adjoints
  for(int i=0; i<lambda.size(); ++i){
    if(i<ocp.numFcn(ASENS_FUN)) // if integrated quantity
      ocp.eval(ASENS_FUN, NV_DATA_S(lambda[i]),0,1);
    else
      N_VConst(0, lambda[i]);
  
    // Set the quadratures to zero
     N_VConst(0, bquad[i]);
  }


  // Integrate backward over all substages
  for(int sub=ocp.nsub()-1; sub>=0; --sub){
      // Start at the end of the substage
      double t = sub+1;
      ocp.setTime(t,sub);
	  
      // Event iteration
      eventIteration(ocp,false);

      // Stop at the beginning of the substage
      double t_stop = sub;

      for(int i=0; i<lambda.size(); ++i){
	// Re-initialize the integration
	sundialsAssert(CVodeReInitB(cvode_mem[sub], adjointID[sub][i], t, lambda[i]));

	// Re-initialize the root-finding
	sundialsAssert(CVodeQuadReInitB(cvode_mem[sub], adjointID[sub][i], bquad[i]));

      }

      // Integrate backwards
      sundialsAssert(CVodeB(cvode_mem[sub], t_stop, CV_NORMAL));

      // Obtain the solution
      for(int i=0; i<lambda.size(); ++i){
	sundialsAssert(CVodeGetB(cvode_mem[sub], adjointID[sub][i], &t, lambda[i]));

	double tret_quad;
	sundialsAssert(CVodeGetQuadB(cvode_mem[sub], adjointID[sub][i], &tret_quad, bquad[i]));

	 std::cout << "quad = " << std::endl;
	 N_VPrint_Serial(bquad[i]);
	 N_VConst(0, bquad[i]);

      }
      std::cout << "t = " << t << std::endl;

      assert(fabs(t-t_stop) < reltol);
    }

#endif

  }


#if 0
assert(0);

  int flag;
  // Turn off root-finding
  flag = CVodeRootInit(cvode_mem, 0, g_static);
  if(flag != CV_SUCCESS) throw "Error in CVodeRootInit";

  // Set the intiial condition
  flag = CVodeReInit(cvode_mem, 0, xx);
  if(flag != CV_SUCCESS) throw "Error in CVodeReInit";

  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(xx));
  double t0 = 0;
  ocp.scaleBackTime(t0);
  ocp.setArg(t_ind,&t0);
  ocp.eventIteration();

  // Stop times
  double stoptimes[] = {1, 6-sqrt(7.2), 6+sqrt(7.2), 11, 12, 24};
  double tret;

  // Integrate over time points
  for(int i=0; i<6; ++i){

    // Set stop time
    flag = CVodeSetStopTime(cvode_mem, stoptimes[i]);
    if(flag != CV_SUCCESS) throw "Error in CVodeSetStopTime"; 
  
    int ncheck;
    flag = CVodeF(cvode_mem, stoptimes[i], xx, &tret, CV_NORMAL, &ncheck);
    std::cout << flag << std::endl;


    if(flag != CV_SUCCESS && flag != CV_TSTOP_RETURN) throw "Error in CVodeF";
  
    printStats();

    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(xx));
    double tret_real = tret;
    ocp.scaleBackTime(tret_real);

    ocp.setArg(t_ind,&tret_real);
    ocp.eventIteration();
  }

  N_VPrint_Serial(xx);

  // Create a memory structure
  int lmmB = ops.lmm;
  int iterB = iter;
  flag = CVodeCreateB(cvode_mem, lmmB, iterB, &adjointID[0]);
  if(flag != CV_SUCCESS) throw "Error in CVodeCreateB";

  std::cout << "adjointID[0]: " << adjointID[0] << std::endl;

  // Initialize (just to allocate memory) - data and time shouldn't matter
//  double tB0 = ocp.tf().getValue();

/*  N_Vector yB0 = N_VNew_Serial(ocp.nx());*/
  N_Vector yB0 = xx;
  N_VConst(0, yB0);
//  x[0] = 1;

  // Initialize the adjoint problem
  double tB0 = stoptimes[5];

  flag = CVodeInitB(cvode_mem, adjointID[0], rhs_asens_static, tB0, yB0);
  if(flag != CV_SUCCESS) throw "Error in CVodeInitB";
  // Set tolerances
  double reltolB = ops.reltol;
  double abstolB = ops.abstol;
  flag = CVodeSStolerancesB(cvode_mem, adjointID[0], ops.reltol, ops.abstol);
  if(flag != CV_SUCCESS) throw "Error in CVodeSStolerancesB";

  // Attach a linear solver  
  int nB = 10;
  flag = CVBandB(cvode_mem, adjointID[0], nB,1,1);
  if(flag != CV_SUCCESS) throw "Error in CVDenseB"; 

  // Set user data
  adjointUD[0].this_ = this;
  adjointUD[0].id = 0;
  flag = CVodeSetUserDataB(cvode_mem, adjointID[0], &adjointUD[0]);
  if(flag != CV_SUCCESS) throw "Error in CVodeSetUserDataB";

  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(xx));
    double tf = 24;
    ocp.setArg(t_ind,&tf);
    ocp.eventIteration();

    flag = CVodeSetStopTime(cvode_mem, stoptimes[4]);
    if(flag != CV_SUCCESS) throw "Error in CVodeSetStopTime"; 

  // Solve the adjoint problem
  flag = CVodeB(cvode_mem, 12.01, CV_NORMAL);


assert(0);


  if(flag != CV_SUCCESS){

  std::cout << "CV_MEM_NULL"  << CV_MEM_NULL << std::endl;
  std::cout << "CV_NO_ADJ"  << CV_NO_ADJ << std::endl;
  std::cout << "CV_NO_BCK"  << CV_NO_BCK << std::endl;
  std::cout << "CV_NO_FWD"  << CV_NO_FWD << std::endl;
  std::cout << "CV_ILL_INPUT"  << CV_ILL_INPUT << std::endl;
/*  std::cout << "CV_BAD_ITASK"  << CV_BAD_ITASK << std::endl;*/
  std::cout << "CV_TOO_MUCH_WORK"  << CV_TOO_MUCH_WORK << std::endl;
  std::cout << "CV_TOO_MUCH_ACC"  << CV_TOO_MUCH_ACC << std::endl;
  std::cout << "CV_ERR_FAILURE"  << CV_ERR_FAILURE << std::endl;
  std::cout << "CV_CONV_FAILURE"  << CV_CONV_FAILURE << std::endl;
  std::cout << "CV_LSETUP_FAIL"  << CV_LSETUP_FAIL << std::endl;
  std::cout << "CV_LSOLVE_FAIL"  << CV_LSOLVE_FAIL << std::endl;
/*  std::cout << "CV_BCKMEM_NULL"  << CV_BCKMEM_NULL << std::endl;*/
/*  std::cout << "CV_BAD_TBOUT"  << CV_BAD_TBOUT << std::endl;*/
  std::cout << "CV_REIFWD_FAIL"  << CV_REIFWD_FAIL << std::endl;
  std::cout << "CV_FWD_FAIL"  << CV_FWD_FAIL << std::endl;

  std::cout << "flag = " << flag << std::endl;
  throw "Error in CVodeB"; 
  }

  // Get the solution
  flag = CVodeGetB(cvode_mem, adjointID[0], &tret, xx);
  std::cout << "flag = " << flag << std::endl;
  if(flag != CV_SUCCESS) throw "Error in CVodeGetB";

  // Print solution
  std::cout << "tret = " << tret << std::endl;
  N_VPrint_Serial(xx);

#endif


void CvodesIntegratorNode::time_event(){
 std::cout << "Time event at t = " << ocp.getConTime() << std::endl;
 throw "Time events not implemented";
}

void CvodesIntegratorNode::state_event(void *mem){

  const std::vector<int> &v_ind = ocp.findVars(VAR_VD);


  // Scale back time
  double t = ocp.getConTime();
  ocp.scaleBackTime(t);

  std::cout << "hello";
  std::cout << ocp.tpval << std::endl;
  std::cout << "state event occured at t = " << t << std::endl;
  std::cout << "hello, sub = " << ocp.getDisTime() << std::endl;


  // Find out which switches fired
  int rootsfound[int(ocp.function(SWITCH).size())];
  sundialsAssert(CVodeGetRootInfo(mem, rootsfound));  

  for(int j=0; j<ocp.function(SWITCH).size(); ++j)
    if(rootsfound[j]){
      // Add the unscaled time
      tp_val[ocp.getDisTime()] << t;

      // Give a suitable name to the variable
      std::ostringstream name;
      name << "Tp_" << int(1000*t);
      Matrix Tp(name.str());
      tp_new[ocp.getDisTime()] << Tp;

      // Evaluate the zero-crossing function to get the definition
      ocp.setArgS(ocp.findVars(VAR_P),ocp.getVars(VAR_P));
      std::vector<double> vvv(v_ind.size());
      ocp.getArg(v_ind,&vvv[0]);
      ocp.setArgS(v_ind,toMatrix(vvv));
      
      Matrix def = ocp.eval(SWITCH)[j];

      tp_def[ocp.getDisTime()] << def;
      break;
    }

#if 0
	// Change the corresponding switches
	for(int j=0; j<ocp.g.size(); ++j)
	  if(rootsfound[j])
	    v[j] = 1-v[j];
#endif

//         ocp.eventIteration(sub,t, &x[0], 0, &ocp.ppp[0], &ocp.vvv[0]); 

std::cout << "hello2";

}


int CvodesIntegratorNode::fQ_static( realtype t, N_Vector y, N_Vector yQdot, void *user_data){
  try{
    // Retrieve the pointer to the object
    CvodesIntegratorNode *this_ = (CvodesIntegratorNode*)user_data;
    OCP_old& ocp = this_->ocp;

    // Pass the time to the ocp
    ocp.setTime(t);

    // Pass the state vector to the ocp
    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

    for(int ord=0; ord<3; ++ord)
      for(int i=0; i<this_->output_ind[ord].size(); ++i){
        int oind = this_->output_ind[ord][i];
        if(ocp.output(oind).has_q)
          ocp.function(OUT_Q,oind).eval(this_->quad_ind[ord][i] + NV_DATA_S(yQdot));
      }

    // Scale time
    for(int i=0; i<NV_LENGTH_S(yQdot); ++i){
      NV_DATA_S(yQdot)[i] *= ocp.getDur();    
    }

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
}

int CvodesIntegratorNode::rhs_asens_static(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB){
#if 0
  try{
    // Retrieve the pointer to the userdata structure
    assert(user_dataB != 0);
    UserDataB* ud = (UserDataB*)user_dataB;
    CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)ud->this_;
    OCP_old& ocp = this_->ocp;

  // Pass the time to the ocp
  ocp.setTime(t);

  // Pass the state vector to the ocp
  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

  // Evaluate
  ocp.eval_times_vector(DE_F, NV_DATA_S(yBdot), NV_DATA_S(yB), 0, 1);

  // Scale time
  for(int i=0; i<NV_LENGTH_S(y); ++i)
    NV_DATA_S(yBdot)[i] *= ocp.getDur();

  // Negate
  for(int i=0; i<NV_LENGTH_S(y); ++i)
    NV_DATA_S(yBdot)[i] *= -1;

  // The adjoints corresponding to integrated quantities get an extra term (equation 2.19 in CVODES user's guide)
  if(ud->id>=ocp.numFcn(ASENS_FUN)){ 
    
    double temp[NV_LENGTH_S(y)];
    ocp.eval(QSENS_FUN, &temp[0],0, 1, ud->id-ocp.numFcn(ASENS_FUN));

    for(int i=0; i<NV_LENGTH_S(y); ++i)
      NV_DATA_S(yBdot)[i] -= temp[i];
   }

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
#endif
}


int CvodesIntegratorNode::rhs_fsens_static(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2){
  try{
    CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)user_data;
    OCP_old& ocp = this_->ocp;

    // This function calculates the right hand side of the sensitivity ODE, i.e. (df/dy)*yS[iS] + (df/dp[iS])
    // Pass the time to the ocp
    ocp.setTime(t);

    // Pass the state vector to the ocp
    ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

    // Evaluate
    ocp.function(DE_F).eval_times_vector(NV_DATA_S(ySdot), NV_DATA_S(yS), 1);

    // Scale time
    for(int i=0; i<NV_LENGTH_S(y); ++i)
      NV_DATA_S(ySdot)[i] *= ocp.getDur();

    ocp.function(DE_F).eval(NV_DATA_S(tmp2), 1, this_->jacf_ind[iS]);

    // Add to the first term
    for(int i=0; i<NV_LENGTH_S(y); ++i) NV_DATA_S(ySdot)[i] += NV_DATA_S(tmp2)[i];

    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }
}

int CvodesIntegratorNode::fQB_static( realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB){
#if ASENS


  try{
    // Retrieve the pointer to the userdata structure
    assert(user_dataB != 0);
    UserDataB* ud = (UserDataB*)user_dataB;
    CvodesIntegratorNode* this_ = (CvodesIntegratorNode*)ud->this_;
    OCP_old& ocp = this_->ocp;

  // Pass the time to the ocp
  ocp.setTime(t);

  // Pass the state vector to the ocp
  ocp.setArg(ocp.findVars(VAR_X),NV_DATA_S(y));

  // Evaluate the gradient of the adjoint function with respect to the parameters
  ocp.eval(QSENS_FUN, NV_DATA_S(qBdot), 0, 1, 1);


  // Evaluate the gradient of of the de rhs with respect to the parameters
  int ind = ocp.findAdjQuad(ud->id);

   double temp[3];
   ocp.eval(DE_F, &temp[0], 0, 1, ind);
  
  NV_DATA_S(qBdot)[0] += temp[0]*NV_DATA_S(yB)[0] + temp[1]*NV_DATA_S(yB)[1] + temp[2]*NV_DATA_S(yB)[2];
  NV_DATA_S(qBdot)[0] *= ocp.getDur();

  //  ocp.eval(DE_F, xdot_data);
    return 0;
  } catch (const char * str){
    std::cerr << str << std::endl;
    return -1;
  }

#endif
}

double CvodesIntegratorNode::nextOutput(){
  double t_next = numeric_limits<double>::infinity();
  
  for(int i=0; i<i_out.size(); ++i){
    // Get a reference
    const vector<double>& t_out = ocp.output(i).tdata;
    
    if(i_out[i]<t_out.size()){
      // candidate for the next output time
      double t_next_test = t_out[i_out[i]];

      // accept if smaller
      if(t_next_test < t_next)
	t_next = t_next_test;
    }
  }

  // Scale to the interval [0,1]
  ocp.scaleTime(t_next);

  return t_next;
}

void CvodesIntegratorNode::saveOutput(){
  // Get the current time and unscale it
  double t_now = ocp.getConTime();
  ocp.scaleBackTime(t_now);

  // loop over all output functions
  for(int ord=0; ord<3; ++ord) // order of sensitivities
    for(int i=0; i<output_ind[ord].size(); ++i){

  
        // The index of the output function
        int o_ind = output_ind[ord][i];

        // The candidate output function
        int t_ind = i_out[o_ind];

        // Get dimensions
        int nf = ocp.function(OUT_G,o_ind).size();
        int nt = ocp.function(OUT_T,o_ind).size();
        int nd = ocp.function(OUT_D,o_ind).size().ncol;
        int np = ocp.output(o_ind).p.size();

        // Get a reference to the output times
        const vector<double>& t_out = ocp.output(o_ind).tdata;

        if(t_ind<t_out.size()){ // if we have not reached the last point
          // candidate for the output time
          double t_test = t_out[t_ind];

          // save if close enough to the current time
          if(fabs(t_test - t_now) <= reltol){

            // Function evaluation
            if(ocp.output(o_ind).has_q){
              double *quad_data = NV_DATA_S(quad);


           double *out_data = &output[o_ind].data(0)[t_ind*ocp.output(o_ind).nf];
/*           double *out_data = &ocp.output(o_ind).fdata[t_ind*ocp.output(o_ind).nf];*/
   //           double *out_data = ocp.output(o_ind).data(t_ind);
              int i_first = quad_ind[ord][i];
              int i_last = quad_ind[ord][i+1];
              for(int j=0; j<i_last-i_first; j++)
                 out_data[j] = quad_data[i_first+j];
            } else {
              // Evaluate and save the output function
              OCPOutput& oo = ocp.output(o_ind);              
              ocp.function(OUT_G,o_ind).eval(&output[o_ind].data()[t_ind*oo.nf]);
            }

            // Save sensitivity information
            if(ord>0){
              int i_first = fsens_ind[ord][i];
              int i_last = fsens_ind[ord][i+1];

              for(int j=0; j<i_last-i_first; j++){
                OCPOutput& oo = ocp.output(o_ind);              
                double *out_data = &output[o_ind].data(1)[t_ind*oo.nf];
                  // Dg/Dp = dg/dp + dg/dx*dx/dp + dg/ddep*ddep/dp

              ocp.function(OUT_G,o_ind).eval(out_data,1,jacoutg_ind[o_ind]);

              const Matrix& g= ocp.function(OUT_G,o_ind);
                  cout << "g = " << g << endl;
                  cout << "o_ind = " << o_ind << endl;

          cout << "j = " << j << endl;
        N_VPrint_Serial(fsens[i_first+j]);


            vector<double> outtmp(nf);
            ocp.function(OUT_G,o_ind).eval_times_vector(&outtmp[0], NV_DATA_S(fsens[i_first+j]), 1, 0);
            for(int k=0; k<outtmp.size(); ++k)
                out_data[i] += outtmp[i];

              }


/*                for(int k=0; k<NV_LENGTH_S(fsens[i_first+j]); ++k){*/
/*                   *out_data = fsens_data[k];
                    out_data++;*/

//           cout << "datasize = " << ocp.getFcnData(OUT_G,o_ind).size() << endl;
//           cout << "order0 = " << nf*nt << endl;
//           cout << "k = " << k << endl;
//           cout << "j = " << k << endl;
//           cout << "sizef = " << NV_LENGTH_S(fsens[0]) << endl;
//           cout << "i_first  = " << i_first << endl;
//           cout << "i_last = " << i_last << endl;

//             assert(nf*nt+nf*nd*t_ind + (j+i_first)*NV_LENGTH_S(fsens[0]) + k <=  ocp.getFcnData(OUT_G,o_ind).size());


/*                }*/
//           }


            }

            // Continue to next output
            i_out[o_ind]++;
          }
       }
    }
}

CvodesIntegrator::CvodesIntegrator(){
}

CvodesIntegrator::CvodesIntegrator(OCP_old &ocp){
  assignNode(new CvodesIntegratorNode(ocp));
}

const CvodesIntegratorNode* CvodesIntegrator::operator->() const{
  return (const CvodesIntegratorNode*)FX::operator->();
}

CvodesIntegratorNode* CvodesIntegrator::operator->(){
  return (CvodesIntegratorNode*)FX::operator->();
}

} // namespace OPTICON

