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
class FXNode;

/** General function
  \author Joel Andersson 
  \date 2010
*/
class FX : public OptionsFunctionality{
  friend class FXNode;

  public:
  /** \brief  default constructor */
  FX(); 

  /** \brief  Destructor */
  ~FX();

  /** \brief  Copy constructor */
  FX(const FX& fx);

  /** \brief  Get number of inputs */
  int getNumInputs() const;

  /** \brief  Get number of outputs */
  int getNumOutputs() const;

  /** \brief  Set number of inputs (invoked automatically) */
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
  void solve(); // the same as evaluate(0,0)
    
  /** \brief Jacobian of output oind with respect to input iind */
  FX jacobian(int iind=0, int oind=0);
  
  /** \brief Hessian of output oind with respect to input iind */
  FX hessian(int iind=0, int oind=0);

  /** \brief  Evaluate symbolically (generate mx node) */
  MX operator()(const MX &x, int ind=0) const;
  MX operator()(const std::vector<MX> &x, int ind=0) const;

  /** \brief  Access functions of the node */
  FXNode* operator->();
  const FXNode* operator->() const;

  /// Assert that the node is pointing to the right type of object
  void assertNode() const;
  
#if 0
// SWIG doesn't support nested template instantiations -- WORKAROUND BELOW
  
  /** \brief  Set an input */
  template <class T> void setInput(T val, int ind=0){input(ind).set(val); }

  /** \brief  Set an output */
  template <class T> void setOutput(T val, int ind=0){output(ind).set(val); }

  /** \brief  Set a forward seed */
  template <class T> void setFwdSeed(T val, int ind=0, int dir=0){input(ind).setF(val,dir); }

  /** \brief  Set a forward sensitivity */
  template <class T> void setFwdSens(T val, int ind=0, int dir=0){output(ind).setF(val,dir); }

  /** \brief  Set an adjoint seed */
  template <class T> void setAdjSeed(T val, int ind=0, int dir=0){output(ind).setA(val,dir); }

  /** \brief  Set an adjoint sensitivity */
  template <class T> void setAdjSens(T val, int ind=0, int dir=0){input(ind).setA(val,dir); }
#else 
// WORKAROUND
#define SETTERS(T)\
  void setInput(T val, int ind=0){input(ind).set(val); } \
  void setOutput(T val, int ind=0){output(ind).set(val); } \
  void setFwdSeed(T val, int ind=0, int dir=0){input(ind).setF(val,dir); } \
  void setFwdSens(T val, int ind=0, int dir=0){output(ind).setF(val,dir); } \
  void setAdjSeed(T val, int ind=0, int dir=0){output(ind).setA(val,dir); } \
  void setAdjSens(T val, int ind=0, int dir=0){input(ind).setA(val,dir); }

SETTERS(double);
SETTERS(const double*);
SETTERS(const std::vector<double>&);
#undef SETTERS
#endif

#if 0
// TODO: THIS DOESN'T WORK, FIND OUT WHY! CARLO?

  /** \brief  Get a input */
  template <class T> void getInput(T val, int ind=0) const {input(ind).get(val); }

  /** \brief  Get a output */
  template <class T> void getOutput(T val, int ind=0) const {output(ind).get(val); }

  /** \brief  Get a forward seed */
  template <class T> void getFwdSeed(T val, int ind=0, int dir=0) const {input(ind).getF(val,dir); }

  /** \brief  Get a forward sensitivity */
  template <class T> void getFwdSens(T val, int ind=0, int dir=0) const {output(ind).getF(val,dir); }

  /** \brief  Get an adjoint seed */
  template <class T> void getAdjSeed(T val, int ind=0, int dir=0) const {output(ind).getA(val,dir); }

  /** \brief  Get an adjoint sensitivity */
  template <class T> void getAdjSens(T val, int ind=0, int dir=0) const {input(ind).getA(val,dir); }
#else

// WORKAROUND
#define GETTERS(T)\
    void getInput(T val, int ind=0) const{input(ind).get(val);} \
    void getOutput(T val, int ind=0) const{output(ind).get(val);} \
    void getFwdSeed(T val, int ind=0, int dir=0) const{input(ind).getF(val,dir);} \
    void getFwdSens(T val, int ind=0, int dir=0) const{output(ind).getF(val,dir);} \
    void getAdjSeed(T val, int ind=0, int dir=0) const{input(ind).getA(val,dir);} \
    void getAdjSens(T val, int ind=0, int dir=0) const{output(ind).getA(val,dir);}

GETTERS(double&);
GETTERS(double*);
GETTERS(std::vector<double>&);
#undef GETTERS
#endif

/** \brief  Get the input data */
const vector<double>& getInputData(int ind=0) const {return input(ind).data(); }

/** \brief  Get the output data */
const vector<double>& getOutputData(int ind=0) const {return output(ind).data(); }

/** \brief  Get the forward seed data  */
const vector<double>& getFwdSeedData(int ind=0, int dir=0) const {return input(ind).dataF(dir); }

/** \brief  Get the forward sensitivity data */
const vector<double>& getFwdSensData(int ind=0, int dir=0) const {return output(ind).dataF(dir); }

/** \brief  Get an adjoint seed data */
const vector<double>& getAdjSeedData(int ind=0, int dir=0) const {return output(ind).dataA(dir); }

/** \brief  Get an adjoint sensitivity data */
const vector<double>& getAdjSensData(int ind=0, int dir=0) const {return input(ind).dataA(dir); }


/** \brief  Get the input data */
vector<double>& getInputData(int ind=0) {return input(ind).data(); }

/** \brief  Get the output data */
vector<double>& getOutputData(int ind=0) {return output(ind).data(); }

/** \brief  Get the forward seed data  */
vector<double>& getFwdSeedData(int ind=0, int dir=0) {return input(ind).dataF(dir); }

/** \brief  Get the forward sensitivity data */
vector<double>& getFwdSensData(int ind=0, int dir=0) {return output(ind).dataF(dir); }

/** \brief  Get an adjoint seed data */
vector<double>& getAdjSeedData(int ind=0, int dir=0) {return output(ind).dataA(dir); }

/** \brief  Get an adjoint sensitivity data */
vector<double>& getAdjSensData(int ind=0, int dir=0) {return input(ind).dataA(dir); }

};

/** \brief Internal class for FX
  \author Joel Andersson 
  \date 2010
 A regular user should never work with any Node class. Use FX directly.
*/
class FXNode : public OptionsFunctionalityNode{
  friend class FX;

  protected:
  /** \brief  Default constructor (accessable from the FX class and derived classes) */
  FXNode();

  public:
  /** \brief  Destructor */
  virtual ~FXNode() = 0;

  /** \brief  Evaluate */
  virtual void evaluate(int fsens_order, int asens_order) = 0;

  /** Initialize and make the object ready for setting arguments and evaluation. This method is typically called after setting options but before evaluating. 
      If passed to another class (in the constructor), this class should invoke this function when initialized. */
  virtual void init();

  /** \brief Jacobian of output oind with respect to input iind */
  virtual FX jacobian(int iind=0, int oind=0);

  /** \brief Hessian of output oind with respect to input iind */
  virtual FX hessian(int iind=0, int oind=0);

  /** \brief  Access an input */
  FunctionIO& input(int i=0);

  /** \brief  Const access an input */
  const FunctionIO& input(int i=0) const;
  
  /** \brief  Access an output*/
  FunctionIO& output(int i=0);

  /** \brief  Const access an output*/
  const FunctionIO& output(int i=0) const;
    
  /** \brief  Print */
  virtual void print(std::ostream &stream) const;
 
  /** \brief  Inputs of the function */
  std::vector<FunctionIO> input_;

  /** \brief  Output of the function */
  std::vector<FunctionIO> output_;

  /// Assert that the function has been initialized
  void assertInit() const;
  
  protected:

  /** \brief  Has the function been initialized? */
  bool is_init_;
    
  int ad_order_;

  /** \brief  Number of forward and adjoint derivatives */
  int nfdir_, nadir_;

};


} // namespace CasADi


#endif // FX_HPP
