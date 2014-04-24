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

#include "ocp_solver_internal.hpp"
#include "integrator.hpp"


INPUTSCHEME(OCPInput)
OUTPUTSCHEME(OCPOutput)

using namespace std;

namespace casadi {

OCPSolverInternal::OCPSolverInternal(const Function& ffcn, const Function& mfcn,
                                     const Function& cfcn, const Function& rfcn) :
    ffcn_(ffcn), mfcn_(mfcn), cfcn_(cfcn), rfcn_(rfcn) {
  addOption("number_of_parameters",  OT_INTEGER,                0);
  addOption("number_of_grid_points", OT_INTEGER,               20);
  addOption("final_time",            OT_REAL,                 1.0);

  input_.scheme = SCHEME_OCPInput;
  output_.scheme = SCHEME_OCPOutput;
}

OCPSolverInternal::~OCPSolverInternal() {
}


void OCPSolverInternal::init() {
  // Initialize the functions
  ffcn_.init();
  mfcn_.init();
  if (!cfcn_.isNull()) cfcn_.init();
  if (!mfcn_.isNull()) mfcn_.init();

  // Get the number of grid points
  nk_ = getOption("number_of_grid_points");

  // Read final time
  tf_ = getOption("final_time");

  // Get the number of differential states
  nx_ = ffcn_.input(DAE_X).size();

  // Get the number of parameters
  np_ = getOption("number_of_parameters");

  // Get the number of controls
  nu_ = ffcn_.input(DAE_P).size() - np_;

  // Number of point constraints
  nh_ = cfcn_.isNull() ? 0 : cfcn_.output().size();

  // Number of point coupling constraints
  ng_ = 0;

  casadi_assert_message(mfcn_.getNumInputs()<=2,
                        "Mayer term must map endstate [ (nx x 1) , (np x 1) ] to cost (1 x 1). "
                        "So it needs to accept 2 matrix-valued inputs. You supplied "
                        << mfcn_.getNumInputs());

  if (mfcn_.getNumInputs()==2) {
      casadi_assert_message(mfcn_.input(1).size()==np_,
                            "Mayer term must map endstate [ (nx x 1) , (np x 1) ] "
                            "to cost (1 x 1). Shape of the second input "
                            << mfcn_.input(1).dimString()
                            << " must match the number of parameters np " << np_);
  }


  // Specify the inputs
  setNumInputs(OCP_NUM_IN);
  input(OCP_LBX) = input(OCP_UBX) = input(OCP_X_INIT) = Matrix<double>::zeros(nx_, nk_+1);
  input(OCP_LBU) = input(OCP_UBU) = input(OCP_U_INIT) = Matrix<double>::zeros(nu_, nk_);
  input(OCP_LBP) = input(OCP_UBP) = input(OCP_P_INIT) = Matrix<double>::zeros(np_);
  input(OCP_LBH) = input(OCP_UBH) = Matrix<double>::zeros(nh_, nk_+1);
  input(OCP_LBG) = input(OCP_UBG) = Matrix<double>::zeros(ng_);

  // Specify the outputs
  setNumOutputs(OCP_NUM_OUT);
  output(OCP_X_OPT) = input(OCP_X_INIT);
  output(OCP_U_OPT) = input(OCP_U_INIT);
  output(OCP_P_OPT) = input(OCP_P_INIT);
  output(OCP_COST) = 0.;

  // Call the init function of the base class
  FunctionInternal::init();


}

} // namespace casadi

