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

#ifndef FX_INTERNAL_HPP
#define FX_INTERNAL_HPP

#include "fx.hpp"
#include <set>

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

namespace CasADi{
  
/** \brief Internal class for FX
  \author Joel Andersson 
  \date 2010
 A regular user should never work with any Node class. Use FX directly.
*/
class FXInternal : public OptionsFunctionalityNode{
  friend class FX;

  protected:
    /** \brief  Default constructor (accessable from the FX class and derived classes) */
    FXInternal();
  
  public:
    /** \brief  Destructor */
    virtual ~FXInternal() = 0;

    /** \brief  Evaluate switch*/
    void evaluate_switch(int nfdir, int nadir);

    /** \brief  Evaluate */
    virtual void evaluate(int nfdir, int nadir) = 0;

    /** \brief Initialize
      Initialize and make the object ready for setting arguments and evaluation. This method is typically called after setting options but before evaluating. 
      If passed to another class (in the constructor), this class should invoke this function when initialized. */
    virtual void init();

    /** \brief  Update the number of sensitivity directions during or after initialization, 
        if recursive==true, updateNumSens is also invoked for the baseclass. */
    virtual void updateNumSens(bool recursive);
    
    /** \brief Calculate the jacobian of a number of function outputs with respect to a number of function inputs, optionally include the function outputs */
    virtual FX jacobian(const std::vector<std::pair<int,int> >& jblocks);

    /** \brief Switch between numeric and symbolic jacobian */
    FX jacobian_switch(const std::vector<std::pair<int,int> >& jblocks);

    /** \brief Numeric Jacobian */
    FX numeric_jacobian(const std::vector<std::pair<int,int> >& jblocks);

    /** \brief Hessian of output oind with respect to input iind */
    virtual FX hessian(int iind=0, int oind=0);

    /** \brief  Propagate the sparsity pattern through a set of directional derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /** \brief  Is the class able to propate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd){ return false;}

    /** \brief  Reset the sparsity propagation */
    virtual void spInit(bool fwd){}
    
    /** \brief  Evaluate symbolically, SX type */
    virtual void evalSX(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
                        const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
                        const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
                        bool output_given, bool eliminate_constants);

    /** \brief  Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& input, std::vector<MX>& output, 
                        const std::vector<std::vector<MX> >& fwdSeed, std::vector<std::vector<MX> >& fwdSens, 
                        const std::vector<std::vector<MX> >& adjSeed, std::vector<std::vector<MX> >& adjSens,
                        bool output_given, bool eliminate_constants);

    /** \brief  Evaluate symbolically, SX type (overloaded)*/
    void eval(const std::vector<SXMatrix>& input, std::vector<SXMatrix>& output, 
              const std::vector<std::vector<SXMatrix> >& fwdSeed, std::vector<std::vector<SXMatrix> >& fwdSens, 
              const std::vector<std::vector<SXMatrix> >& adjSeed, std::vector<std::vector<SXMatrix> >& adjSens,
              bool output_given, bool eliminate_constants){
      evalSX(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,eliminate_constants);
    }

    /** \brief  Evaluate symbolically, MX type (overloaded)*/
    void eval(const std::vector<MX>& input, std::vector<MX>& output, 
              const std::vector<std::vector<MX> >& fwdSeed, std::vector<std::vector<MX> >& fwdSens, 
              const std::vector<std::vector<MX> >& adjSeed, std::vector<std::vector<MX> >& adjSens,
              bool output_given, bool eliminate_constants){
      evalMX(input,output,fwdSeed,fwdSens,adjSeed,adjSens,output_given,eliminate_constants);
    }
    
    /** \brief  Access an input */
    FunctionIO& inputStruct(int i=0);

    /** \brief  Const access an input */
    const FunctionIO& inputStruct(int i=0) const;
    
    /** \brief  Access an output*/
    FunctionIO& outputStruct(int i=0);

    /** \brief  Const access an output*/
    const FunctionIO& outputStruct(int i=0) const;
      
    /** \brief  Print */
    virtual void print(std::ostream &stream) const;
    
    /** \brief  Print */
    virtual void repr(std::ostream &stream) const;
  
    /** \brief  Inputs of the function */
    std::vector<FunctionIO> input_;

    /** \brief  Output of the function */
    std::vector<FunctionIO> output_;

    /** \brief Get the unidirectional or bidirectional partition */
    void getPartition(int iind, int oind, CRSSparsity& D1, CRSSparsity& D2, bool compact, bool symmetric);

    /// Verbose mode?
    bool verbose() const;
    
    /// Is function fcn being monitored
    bool monitored(const std::string& mod) const;
    
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

    /// Set the number of function inputs
    void setNumInputs(int num_in);

    /// Set the number of function outputs
    void setNumOutputs(int num_out);

    /// Get the number of function inputs
    int getNumInputs() const;

    /// Get the number of function outputs
    int getNumOutputs() const;
    
    /// Get all statistics obtained at the end of the last evaluate call
    const Dictionary & getStats() const;

    /// Get single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string & name) const;
    
    /// Generate the sparsity of a Jacobian block
    virtual CRSSparsity getJacSparsity(int iind, int oind);
    
    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const CRSSparsity& sp, int iind, int oind, bool compact);
    
    /// Get, if necessary generate, the sparsity of a Jacobian block
    CRSSparsity& jacSparsity(int iind, int oind, bool compact);
    
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<MX> symbolicInput() const;
  
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<SXMatrix> symbolicInputSX() const;
  
    /// Get the Jacobian of all outputs with respect to all inputs
    void getFullJacobian();

    /** \brief  Number of forward and adjoint derivatives */
    int nfdir_, nadir_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;
    
    /** \brief  Log the status of the solver */
    void log(const std::string& msg) const;

    /** \brief  Log the status of the solver, function given */
    void log(const std::string& fcn, const std::string& msg) const;

    /// Set of module names which are extra monitored
    std::set<std::string> monitors_;
    
    /** \brief  Dictionary of statistics (resulting from evaluate) */
    Dictionary stats_;

    /** \brief  Stored Jacobians */
    bool store_jacobians_;
    std::vector<std::vector<FX> > jacs_;
    
    /// Sparsity of the Jacobian blocks
    std::vector<std::vector<CRSSparsity> > jac_sparsity_, jac_sparsity_compact_;

    /// Use numeric jacobian instead of symbolic
    bool numeric_jacobian_;
    
    /// User-provided Jacobian generator function
    JacobianGenerator jacgen_;

    /// User-provided sparsity generator function
    SparsityGenerator spgen_;

    /// User-set field
    void* user_data_;
    
    /// Full jacobian function used to calculate directional derivatives instead of the using directional derivatives
    bool jac_for_sens_;
    FX full_jacobian_;
    
    bool monitor_inputs_;
    
    bool monitor_outputs_;
};


} // namespace CasADi


#endif // FX_INTERNAL_HPP
