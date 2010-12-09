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

#ifndef SUNDIALS_INTEGRATOR_HPP
#define SUNDIALS_INTEGRATOR_HPP

#include "../integrator.hpp"

#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h>  /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>  /* definition of type realtype */

namespace OPTICON{

/** \brief  Forward declaration */
class SundialsIntegrator;

/** \brief  Userdata structure to identify the different backward problems */
struct UserDataB{
  SundialsIntegrator* this_;
  int id;
}; 

class SundialsIntegrator : public IntegratorNode_old{
public:
/** \brief  Constructor */
  SundialsIntegrator(OCP_old &ocp);

/** \brief  Destructor */
  virtual ~SundialsIntegrator() = 0;

/** \brief  Initialize stage */
  virtual void init();

/** \brief  Integrate over the stage - supports zero-crossing functions and user output */
  virtual void evaluate(int tape_order=0) = 0;

/** \brief  Integrate the problem - possibly prepare for solving the adjoint problem, no root-finding not supported! */
  virtual void evaluateFwd(bool use_tape=false) = 0;
  virtual void evaluateAdj() = 0;

/** \brief  Print solver statistics */
  virtual void printStats(std::ostream &stream=std::cout) const = 0 ;

/** \brief  Id for the adjoint problems */
  std::vector<std::vector<int> > adjointID;

/** \brief  Userdata for the adjoint problems */
  std::vector<std::vector<UserDataB> > adjointUD;


protected:

/** \brief  Check error flags of SUNDIALS functions */
  void sundialsAssert(int flag) const;


};

} // namespace OPTICON


#endif //SUNDIALS_INTEGRATOR_HPP
