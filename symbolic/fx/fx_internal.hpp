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
#include "../weak_ref.hpp"
#include <set>
#include "code_generator.hpp"

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
  class FXInternal : public OptionsFunctionalityNode, public IOInterface<FXInternal>{
    friend class FX;

  protected:
    /** \brief  Default constructor (accessable from the FX class and derived classes) */
    FXInternal();
  
  public:
    /** \brief  Destructor */
    virtual ~FXInternal() = 0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief  Evaluate */
    virtual void evaluate(int nfdir, int nadir) = 0;
  
    /** \brief  Evaluate with directional derivative compression */
    void evaluateCompressed(int nfdir, int nadir);

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation. This method is typically called after setting options but before evaluating. 
        If passed to another class (in the constructor), this class should invoke this function when initialized. */
    virtual void init();

    /** \brief  Update the number of sensitivity directions during or after initialization, 
        if recursive==true, updateNumSens is also invoked for the baseclass. */
    virtual void updateNumSens(bool recursive);
    
    /** \brief Request a number of forward/adjoint derivative directions */
    void requestNumSens(int nfwd, int nadj);
  
    /** \brief  Propagate the sparsity pattern through a set of directional derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /** \brief  Is the class able to propate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd){ return false;}

    /** \brief  Reset the sparsity propagation */
    virtual void spInit(bool fwd){}
    
    /** \brief  Evaluate symbolically, SX type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
                        const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
                        const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens);

    /** \brief  Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
                        const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                        const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens);

    /** \brief  Evaluate symbolically, SX type, matching sparsity patterns */
    virtual void evalSXsparse(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
                              const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
                              const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens);
    
    /** \brief Call a function, MX type (overloaded) */
    void call(const MXVector& arg, MXVector& res, 
              const MXVectorVector& fseed, MXVectorVector& fsens, 
              const MXVectorVector& aseed, MXVectorVector& asens,
              bool always_inline, bool never_inline);
    
    /** \brief Call a function, SX type (overloaded) */
    void call(const std::vector<SXMatrix>& arg, std::vector<SXMatrix>& res, 
              const std::vector<std::vector<SXMatrix> >& fseed, std::vector<std::vector<SXMatrix> >& fsens, 
              const std::vector<std::vector<SXMatrix> >& aseed, std::vector<std::vector<SXMatrix> >& asens,
              bool always_inline, bool never_inline);
        
    //@{
    /** \brief Return Hessian function */
    FX hessian(int iind, int oind);
    virtual FX getHessian(int iind, int oind);
    //@}

    //@{
    /** \brief Return gradient function */
    FX gradient(int iind, int oind);
    virtual FX getGradient(int iind, int oind);
    //@}

    //@{
    /** \brief Return Jacobian function */
    FX jacobian(int iind, int oind, bool compact, bool symmetric);
    virtual FX getJacobian(int iind, int oind, bool compact, bool symmetric);
    virtual FX getNumericJacobian(int iind, int oind, bool compact, bool symmetric);
    //@}
    
    //@{
    /** \brief Return Jacobian of all input nonzeros with respect to all output nonzeros */
    FX fullJacobian();
    virtual FX getFullJacobian();
    //@}

    //@{
    /** \brief Return function that calculates forward derivatives 
     *    This method returns a cached instance if available, and calls FX getDerivative(int nfwd, int nadj) if no cached version is available.
     */
    FX derivative(int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives */
    virtual FX getDerivative(int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives by creating the Jacobian then multiplying */
    virtual FX getDerivativeViaJac(int nfwd, int nadj);

    //@}

    /** \brief  Print to a c file */
    virtual void generateCode(const std::string& filename);

    /** \brief Generate code for function inputs and outputs */
    void generateIO(CodeGenerator& gen);

    /** \brief Generate code the functon */
    virtual void generateFunction(std::ostream &stream, const std::string& fname, const std::string& input_type, const std::string& output_type, const std::string& type, CodeGenerator& gen) const;
    
    /** \brief Generate code for the declarations of the C function */
    virtual void generateDeclarations(std::ostream &stream, const std::string& type, CodeGenerator& gen) const;

    /** \brief Generate code for the function body */
    virtual void generateBody(std::ostream &stream, const std::string& type, CodeGenerator& gen) const;

    /** \brief  Print */
    virtual void print(std::ostream &stream) const;
    
    /** \brief  Print */
    virtual void repr(std::ostream &stream) const;
            
    /** \brief Get the unidirectional or bidirectional partition */
    void getPartition(int iind, int oind, CRSSparsity& D1, CRSSparsity& D2, bool compact, bool symmetric);

    /// Verbose mode?
    bool verbose() const;
    
    /// Is function fcn being monitored
    bool monitored(const std::string& mod) const;
        
    /// Get total number of scalar inputs (i.e. the number of nonzeros in all of the matrix-valued inputs)
    int getNumScalarInputs() const;

    /// Get total number of scalar outputs (i.e. the number of nonzeros in all of the matrix-valued outputs)
    int getNumScalarOutputs() const;
    
    /// Get all statistics obtained at the end of the last evaluate call
    const Dictionary & getStats() const;

    /// Get single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string & name) const;
    
    /// Generate the sparsity of a Jacobian block
    virtual CRSSparsity getJacSparsity(int iind, int oind, bool symmetric);
    
    /// A flavour of getJacSparsity without any magic
    CRSSparsity getJacSparsityPlain(int iind, int oind);
    
    /// A flavour of getJacSparsity that does hierachical block structure recognition
    CRSSparsity getJacSparsityHierarchical(int iind, int oind);
    
    /// A flavour of getJacSparsity that does hierachical block structure recognition for symmetric jacobians
    CRSSparsity getJacSparsityHierarchicalSymm(int iind, int oind);
    
    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const CRSSparsity& sp, int iind, int oind, bool compact);
    
    /// Get, if necessary generate, the sparsity of a Jacobian block
    CRSSparsity& jacSparsity(int iind, int oind, bool compact, bool symmetric);
    
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<MX> symbolicInput() const;
  
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<SXMatrix> symbolicInputSX() const;
  
    // Workaround helper functions: assign nonzeros but ignore all -1
    static void assignIgnore(MX& y, const MX& x, const std::vector<int>& nz);
    static void assignIgnore(SXMatrix& y, const SXMatrix& x, const std::vector<int>& nz);

    //@{
    /** \brief Access input/output scheme */
    inline const InputOutputScheme& inputScheme() const{ return inputScheme_;}
    inline const InputOutputScheme& outputScheme() const{ return outputScheme_;}
    inline InputOutputScheme& inputScheme(){ return inputScheme_;}
    inline InputOutputScheme& outputScheme(){ return outputScheme_;}
    //@}

    //@{
    /// Input/output structures of the function */
    inline const std::vector<FunctionIO>& input_struct() const{ return input_;}
    inline const std::vector<FunctionIO>& output_struct() const{ return output_;}
    inline std::vector<FunctionIO>& input_struct(){ return input_;}
    inline std::vector<FunctionIO>& output_struct(){ return output_;}
    //@}

    //@{
    /// Input/output access without checking (faster, but unsafe)
    inline const Matrix<double>& inputNoCheck(int iind=0) const{ return inputS<false>(iind).data;}
    inline const Matrix<double>& outputNoCheck(int oind=0) const{ return outputS<false>(oind).data;}
    inline const Matrix<double>& fwdSeedNoCheck(int iind=0, int dir=0) const{ return const_cast<FXInternal*>(this)->fwdSeedNoCheck(iind,dir); }
    inline const Matrix<double>& fwdSensNoCheck(int oind=0, int dir=0) const{ return const_cast<FXInternal*>(this)->fwdSensNoCheck(oind,dir);}    
    inline const Matrix<double>& adjSeedNoCheck(int oind=0, int dir=0) const{ return const_cast<FXInternal*>(this)->adjSeedNoCheck(oind,dir);}
    inline const Matrix<double>& adjSensNoCheck(int iind=0, int dir=0) const{ return const_cast<FXInternal*>(this)->adjSensNoCheck(iind,dir);}

    inline Matrix<double>& inputNoCheck(int iind=0){ return inputS<false>(iind).data;}
    inline Matrix<double>& outputNoCheck(int oind=0){ return outputS<false>(oind).data;}
    inline Matrix<double>& fwdSeedNoCheck(int iind=0, int dir=0){ return inputS<false>(iind).dataF[dir]; }
    inline Matrix<double>& fwdSensNoCheck(int oind=0, int dir=0){ return outputS<false>(oind).dataF[dir]; }
    inline Matrix<double>& adjSeedNoCheck(int oind=0, int dir=0){ return outputS<false>(oind).dataA[dir];}
    inline Matrix<double>& adjSensNoCheck(int iind=0, int dir=0){ return inputS<false>(iind).dataA[dir];}
    //@}

    /** \brief  Log the status of the solver */
    void log(const std::string& msg) const;

    /** \brief  Log the status of the solver, function given */
    void log(const std::string& fcn, const std::string& msg) const;

    // Codegen function
    FX dynamicCompilation(FX f, std::string fname, std::string fdescr, std::string compiler);

    /** \brief  Inputs of the function */
    std::vector<FunctionIO> input_;

    /** \brief  Output of the function */
    std::vector<FunctionIO> output_;

    /** \brief  Number of forward and adjoint derivatives */
    int nfdir_, nadir_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;
    
    /// Set of module names which are extra monitored
    std::set<std::string> monitors_;
    
    /** \brief  Dictionary of statistics (resulting from evaluate) */
    Dictionary stats_;
    
    /** \brief  Flag to indicate wether statistics must be gathered */
    bool gather_stats_;

    /// Cache for functions to evaluate directional derivatives
    std::vector<std::vector<FX> > derivative_fcn_; // NOTE: This can result in circular dependencies!

    /// Cache for full Jacobian
    WeakRef full_jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    Matrix<CRSSparsity> jac_sparsity_, jac_sparsity_compact_;

    /// Which derivative directions are currently being compressed
    std::vector<bool> compressed_fwd_, compressed_adj_;

    /// User-provided Jacobian generator function
    JacobianGenerator jacgen_;

    /// User-provided sparsity generator function
    SparsityGenerator spgen_;

    /// User-set field
    void* user_data_;
    
    bool monitor_inputs_, monitor_outputs_;
    
    /// The name of the input scheme of this function
    InputOutputScheme inputScheme_;
    
    /// The name of the output scheme of this function
    InputOutputScheme outputScheme_;
    
    /// Errors are thrown when NaN is produced
    bool regularity_check_;
    
  };


} // namespace CasADi


#endif // FX_INTERNAL_HPP
