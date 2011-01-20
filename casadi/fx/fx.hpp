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

  /// Access input argument
  Matrix<double>& input(int iind=0);
    
  /// Const access input argument
  const Matrix<double>& input(int iind=0) const;

  /// Access input argument
  Matrix<double>& output(int oind=0);
    
  /// Const access input argument
  const Matrix<double>& output(int oind=0) const;

  /// Access forward seed
  Matrix<double>& fwdSeed(int iind=0, int dir=0);
    
  /// Const access forward seed
  const Matrix<double>& fwdSeed(int iind=0, int dir=0) const;

  /// Access forward sensitivity
  Matrix<double>& fwdSens(int oind=0, int dir=0);
    
  /// Const access forward sensitivity
  const Matrix<double>& fwdSens(int oind=0, int dir=0) const;

  /// Access adjoint seed
  Matrix<double>& adjSeed(int oind=0, int dir=0);
    
  /// Const access adjoint seed
  const Matrix<double>& adjSeed(int oind=0, int dir=0) const;

  /// Access forward sensitivity
  Matrix<double>& adjSens(int iind=0, int dir=0);
    
  /// Const access forward sensitivity
  const Matrix<double>& adjSens(int iind=0, int dir=0) const;

  /// Add modules to be monitored
  void addMonitor(const std::string& mon);
  
  /// Remove modules to be monitored
  void removeMonitor(const std::string& mon);
  
  
#ifdef DOXYGENPROC
/// \name Setters
/// Set an input, output, forward seed/sensitivity or adjoint seed/sensitivity\n
/// T can be double&, double*, std::vector<double>&, Matrix<double> &\n
/// Assumes a properly allocated val.\n
/// @{
/** 
    \brief Reads in the input argument from val.
*/
void setInput(T val, int ind=0) const;
/** 
    \brief Reads in the output argument from val.
*/
void setOutput(T val, int ind=0) const;
/** 
    \brief Reads in the forward seed from val.
*/
void setFwdSeed(T val,  int ind=0, int dir=0) const;
/** 
    \brief Reads in the forward sensitivity from val.
*/
void setFwdSens(T val, int ind=0, int dir=0) const ;
/** 
    \brief Reads in the adjoint seed from val.
*/
void setAdjSeed(T val,  int ind=0, int dir=0) const;
/** 
    \brief Reads in the adjoint sensitivity from val.
*/
void setAdjSens(T val, int ind=0, int dir=0) const ;
/// @}

#endif

#define SETTERS(T)\
  void setInput(T val, int ind=0)             { assertInit(); input(ind).set(val);  } \
  void setOutput(T val, int ind=0)            { assertInit(); input(ind).set(val); } \
  void setFwdSeed(T val, int ind=0, int dir=0){ assertInit(); fwdSeed(ind,dir).set(val); } \
  void setFwdSens(T val, int ind=0, int dir=0){ assertInit(); fwdSeed(ind,dir).set(val); } \
  void setAdjSeed(T val, int ind=0, int dir=0){ assertInit(); adjSeed(ind,dir).set(val); } \
  void setAdjSens(T val, int ind=0, int dir=0){ assertInit(); adjSens(ind,dir).set(val); }

/// \cond
SETTERS(double);
SETTERS(const double*);
SETTERS(const std::vector<double>&);
/// \endcond
#undef SETTERS

#define GETTERS(T)\
    void getInput(T val, int ind=0) const             { assertInit(); input(ind).get(val);} \
    void getOutput(T val, int ind=0) const            { assertInit(); output(ind).get(val);} \
    void getFwdSeed(T val, int ind=0, int dir=0) const{ assertInit(); fwdSeed(ind,dir).get(val);} \
    void getFwdSens(T val, int ind=0, int dir=0) const{ assertInit(); fwdSens(ind,dir).get(val);} \
    void getAdjSeed(T val, int ind=0, int dir=0) const{ assertInit(); adjSeed(ind,dir).get(val);} \
    void getAdjSens(T val, int ind=0, int dir=0) const{ assertInit(); adjSens(ind,dir).get(val);}
/// \cond
GETTERS(double&);
GETTERS(double*);
GETTERS(std::vector<double>&);
/// \endcond
#undef GETTERS

#ifdef DOXYGENPROC
/// \name Getters
/// A group of accessor for numerical data that operate on preallocated data.\n
/// get an input, output, forward seed/sensitivity or adjoint seed/sensitivity\n
/// T can be double&, double*, std::vector<double>&, Matrix<double> &\n
/// Assumes a properly allocated val.\n
/// @{

/** \brief Writes out the input argument into val.
*/
void getInput(T val, int ind=0) const;
 
/** 
    \brief Writes out the output argument into val.
*/
void getOutput(T val, int ind=0) const;

/** 
    \brief Writes out the forward seed into val.
*/
void getFwdSeed(T val,  int ind=0, int dir=0) const;

/**  
    \brief Writes out the forward sensitivity into val.
*/
void getFwdSens(T val, int ind=0, int dir=0) const;
/** 
    \brief Writes out the adjoint seed into val.
*/
void getAdjSeed(T val,  int ind=0, int dir=0) const ;

/** 
    \brief Writes out the adjoint sensitivity into val.
*/
void getAdjSens(T val, int ind=0, int dir=0) const;
/// @}
#endif


};
} // namespace CasADi


#endif // FX_HPP
