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

#ifndef FX_HPP
#define FX_HPP

#include "../mx/mx.hpp"
#include <vector>
#include <string>
#include "../options_functionality.hpp"
#include "function_io.hpp"

namespace CasADi{
  
/** Forward declaration of internal class */
class FXInternal;

/** \brief General function

  A general function \f$f\f$ in casadi can be multi-input, multi-output.\n
  Number of inputs:  nin    getNumInputs()\n
  Number of outputs: nout   getNumOutputs()\n
  
  We can view this function as a being composed of a (nin, nout) grid of single-input, single-output primitive functions.\n
  Each such primitive function \f$f_{i,j} \forall i \in [0,nin-1], j \in [0,nout-1]\f$ can map as \f$\mathbf{R}^{n,m}\to\mathbf{R}^{p,q}\f$, 
  in which n,m,p,q can take different values for every (i,j) pair.\n
  
  When passing input, you specify which partition i is active.     You pass the numbers flattened, as a vector of size \f$(n*m)\f$.\n
  When requesting output, you specify which partition j is active. You get the numbers flattened, as a vector of size \f$(p*q)\f$.\n
  
  To calculate jacobians, you need to have \f$(m=1,q=1)\f$.
  
  Write the jacobian as \f$J_{i,j} = \nabla f_{i,j} = \frac{\partial f_{i,j}(\vec{x})}{\partial \vec{x}}\f$.
  
  Using \f$\vec{v} \in \mathbf{R}^n\f$ as a forward seed:  setFwdSeed(v,i)\n
  Retrieving \f$\vec{s}_f \in \mathbf{R}^p\f$ from:        getFwdSens(sf,j)\n
  
  Using \f$\vec{w} \in \mathbf{R}^p\f$ as a forward seed:  setAdjSeed(w,j)\n
  Retrieving \f$\vec{s}_a \in \mathbf{R}^n \f$ from:        getAdjSens(sa,i)\n
  
  We have the following relationships:
  
  \f$ \vec{s}_f = \nabla f_{i,j} . \vec{v}\f$ \n
  \f$ \vec{s}_a = (\nabla f_{i,j})^T . \vec{w}\f$
  
  \section Notes for developers
  
  Each function consists of 4 files:\n
  1. public class header file: imported in python\n
  2. public class implementation\n
  3. internal class header file: should only be used by derived classes\n
  4. internal class implementation\n
  
  python and c++ should be 1-to-1\n
  There should be no extra features in 1.\n
  All the functionality should exist in 1.\n
  If it means that c++ will be more "pythonic", so be it.
  
  \author Joel Andersson 
  \date 2010
*/
class FX : public OptionsFunctionality{

  public:
  /** \brief  default constructor */
  FX(); 

  /** \brief  Destructor */
  ~FX();

  /** \brief  Get number of inputs */
  int getNumInputs() const;

  /** \brief  Get number of outputs */
  int getNumOutputs() const;

  /** \brief  Set number of inputs (invoked autmatically) */
  void setNumInputs(int num_in);

  /** \brief  Set number of outputs  (invoked automatically) */
  void setNumOutputs(int num_out);

  /** \brief  Access an input */
  FunctionIO& input(int i=0);

  /** \brief  Const access an input */
  const FunctionIO& input(int i=0) const;
  
  /** \brief  Access an output*/
  FunctionIO& output(int i=0);

  /** \brief  Const access an output*/
  const FunctionIO& output(int i=0) const;
  
  /** \brief  Evaluate */
  void evaluate(int fsens_order=0, int asens_order=0);
  
  /// the same as evaluate(0,0)
  void solve(); 
    
  /** \brief Jacobian of output oind with respect to input iind */
  FX jacobian(int iind=0, int oind=0);
  
  /** \brief Hessian of output oind with respect to input iind */
  FX hessian(int iind=0, int oind=0);

#ifndef SWIG
  /** \brief  Create a function call (evaluation mx node), single input */
  std::vector<MX> call(const MX &x) const;
#endif // SWIG

  /** \brief  Create a function call (evaluation mx node) */
  std::vector<MX> call(const std::vector<MX> &x) const;

  // Legacy code: change for something else, but what??
#ifndef USE_FUNCTORS
  /** \brief  Create a function call (generate mx node), single input: DEPRECIATED, USE "call" instead */
  MX operator()(const MX &x, int ind=0) const;

  /** \brief  Create a function call (generate mx node): DEPRECIATED, USE "call" instead. */
  MX operator()(const std::vector<MX> &x, int ind=0) const;
#endif // USE_FUNCTORS
  
  /** \brief  Access functions of the node */
  FXInternal* operator->();

  /** \brief  Const access functions of the node */
  const FXInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Is initialized?
  bool isInit() const;
  
  /// Assert that the function has been initialized
  void assertInit() const;
  
  /// Set an input, output, forward seed/sensitivity or adjoint seed/sensitivity
#define SETTERS(T)\
  void setInput(T val, int ind=0)             { assertInit(); input(ind).set(val);  } \
  void setOutput(T val, int ind=0)            { assertInit(); output(ind).set(val); } \
  void setFwdSeed(T val, int ind=0, int dir=0){ assertInit(); input(ind).set(val,1+dir); } \
  void setFwdSens(T val, int ind=0, int dir=0){ assertInit(); output(ind).set(val,1+dir); } \
  void setAdjSeed(T val, int ind=0, int dir=0){ assertInit(); output(ind).set(val,-1-dir); } \
  void setAdjSens(T val, int ind=0, int dir=0){ assertInit(); input(ind).set(val,-1-dir); }

SETTERS(double);
SETTERS(const double*);
SETTERS(const std::vector<double>&);
#undef SETTERS

#ifdef SWIG
// Forward renaming declarations
%rename(getInput) getInputData;
%rename(getOutput) getOutputData;
%rename(getFwdSeed) getFwdSeedData;
%rename(getFwdSens) getFwdSensData;
%rename(getAdjSeed) getAdjSeedData;
%rename(getAdjSens) getAdjSensData;
#else // SWIG
#define GETTERS(T)\
    void getInput(T val, int ind=0) const             { assertInit(); input(ind).get(val);} \
    void getOutput(T val, int ind=0) const            { assertInit(); output(ind).get(val);} \
    void getFwdSeed(T val, int ind=0, int dir=0) const{ assertInit(); input(ind).get(val,1+dir);} \
    void getFwdSens(T val, int ind=0, int dir=0) const{ assertInit(); output(ind).get(val,1+dir);} \
    void getAdjSeed(T val, int ind=0, int dir=0) const{ assertInit(); output(ind).get(val,-1-dir);} \
    void getAdjSens(T val, int ind=0, int dir=0) const{ assertInit(); input(ind).get(val,-1-dir);}
GETTERS(double&);
GETTERS(double*);
GETTERS(std::vector<double>&);
#undef GETTERS
#endif // SWIG

/** \brief  Get the input data */
const std::vector<double>& getInputData(int ind=0) const;

/** \brief  Get the output data */
const std::vector<double>& getOutputData(int ind=0) const;

/** \brief  Get the forward seed data  */
const std::vector<double>& getFwdSeedData(int ind=0, int dir=0) const;

/** \brief  Get the forward sensitivity data */
const std::vector<double>& getFwdSensData(int ind=0, int dir=0) const;

/** \brief  Get an adjoint seed data */
const std::vector<double>& getAdjSeedData(int ind=0, int dir=0) const;

/** \brief  Get an adjoint sensitivity data */
const std::vector<double>& getAdjSensData(int ind=0, int dir=0) const;

// Not in the SWIG interface since it causes problems with the typemap
#ifndef SWIG

/** \brief  Get the input data */
std::vector<double>& getInputData(int ind=0);

/** \brief  Get the output data */
std::vector<double>& getOutputData(int ind=0);

/** \brief  Get the forward seed data  */
std::vector<double>& getFwdSeedData(int ind=0, int dir=0);

/** \brief  Get the forward sensitivity data */
std::vector<double>& getFwdSensData(int ind=0, int dir=0);

/** \brief  Get an adjoint seed data */
std::vector<double>& getAdjSeedData(int ind=0, int dir=0);

/** \brief  Get an adjoint sensitivity data */
std::vector<double>& getAdjSensData(int ind=0, int dir=0);

#endif // SWIG


};
} // namespace CasADi


#endif // FX_HPP
