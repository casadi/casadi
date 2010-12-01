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

#ifndef CVODES_INTEGRATOR_HPP
#define CVODES_INTEGRATOR_HPP

#include "sundials_integrator.hpp"

#include <cvodes/cvodes.h>            /* prototypes for CVODE fcts. and consts. */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h> 
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>

namespace OPTICON{

class CvodesIntegratorNode : public SundialsIntegrator{
friend class CvodesIntegrator;
public:
/** \brief  Create function (calls the constructor) */
  virtual CvodesIntegratorNode* create(OCP_old& ocp) const;

/** \brief  Static crator function (calls the constructor) */
  static IntegratorNode_old* creator(OCP_old& ocp);
  
/** \brief  Destructor */
  virtual ~CvodesIntegratorNode();

/** \brief  Initialize stage */
  virtual void init();

/** \brief  Integrate over the stage - supports zero-crossing functions and user output */
  virtual void evaluate(int tape_order=0);

/** \brief  Integrate the problem */
  virtual void evaluateFwd(bool use_tape=false);
  virtual void evaluateAdj();

/** \brief  Get the initial condition for the state given parameter values */
  void getInitialCondition();

/** \brief  Print solver statistics */
  void printStats(std::ostream &stream=std::cout) const;

/** \brief  Pointers to CVode memory blocks */
  std::vector<void *> cvode_mem;

/** \brief  Vectors */
  N_Vector xx;

/** \brief  quadratures */
  N_Vector quad;
  std::vector<int> quad_ind[3]; 

/** \brief  forward sensitivities */
  std::vector<N_Vector> fsens;
  std::vector<int> fsens_ind[3]; 

/** \brief  index of the jacobian associated with each forward sensitivity */
  std::vector<int> jacf_ind;
  std::vector<int> jacic_ind;
  std::vector<int> jacoutg_ind;

#if 0
  std::vector<N_Vector> bquad;

/** \brief  Lagrange multipliers */
  std::vector<N_Vector> lambda;

/** \brief  For adjoint sensitivities */
  std::vector<N_Vector> asens;
  
//   // For adjoint sensitivities of quadrature
/** \brief    std::vector<N_Vector> asensq; */

#endif


/** \brief  CVODES user functions (explainations needed) */
  static int f_static(realtype t, N_Vector nv_x, N_Vector nv_xdot, void *user_data);
  static void ehfun_static(int error_code, const char *module, const char *function, char *msg, void *eh_data);
  static int djac_static(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 
  static int bjac_static(int N, int mupper, int mlower, 	 
 		 	realtype t, N_Vector y, N_Vector fy, 	 
 		 	DlsMat Jac, void *user_data, 	 
 		 	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3); 
  static int jtimes_static( N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);
  static int g_static( realtype t, N_Vector y, realtype *gout, void *user_data);
  static int fQ_static( realtype t, N_Vector y, N_Vector yQdot, void *user_data);
  static int rhs_fsens_static(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2); 
  static int rhs_asens_static(realtype t, N_Vector y, N_Vector yB, N_Vector yBdot, void *user_dataB);
  static int fQB_static( realtype t, N_Vector y, N_Vector yB, N_Vector qBdot, void *user_dataB);

/** \brief  handle time event */
  void time_event(); 

/** \brief  handle state event */
  void state_event(void *mem); 

/** \brief  These functions should be removed!  */
  std::vector<Matrix> tp_new, tp_def, tp_val;


protected:
/** \brief  Constructor (only to be called from the CvodesIntegrator class) */
  explicit CvodesIntegratorNode(OCP_old &ocp);

/** \brief  Create and intiialize a CVode memory block */
  void initCvodes(int k);

/** \brief  Get the time of the next output */
  double nextOutput();  
  
/** \brief  Save output to ocp */
  void saveOutput();

};

/** \brief  Pointer class */
class CvodesIntegrator : public Integrator_old{
public:
   
/** \brief  Constructors */
  explicit CvodesIntegrator(); // default constuctor
  explicit CvodesIntegrator(OCP_old &ocp);

/** \brief  Access functions of the node */
  CvodesIntegratorNode* operator->();
  const CvodesIntegratorNode* operator->() const;
};

} // namespace OPTICON


#endif //CVODES_INTEGRATOR_HPP
