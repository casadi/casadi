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
#include "../functor.hpp"
#include "../weak_ref.hpp"
#include <set>
#include "code_generator.hpp"

// This macro is for documentation purposes
#define INPUTSCHEME(name)

// This macro is for documentation purposes
#define OUTPUTSCHEME(name)

namespace CasADi{

  class MXFunction;
  
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
    virtual void evaluate() = 0;

    /** \brief Initialize
        Initialize and make the object ready for setting arguments and evaluation. This method is typically called after setting options but before evaluating. 
        If passed to another class (in the constructor), this class should invoke this function when initialized. */
    virtual void init();

    /** \brief  Propagate the sparsity pattern through a set of directional derivatives forward or backward */
    virtual void spEvaluate(bool fwd);

    /** \brief  Propagate the sparsity pattern through a set of directional derivatives forward or backward, using the sparsity patterns */
    virtual void spEvaluateViaJacSparsity(bool fwd);

    /** \brief  Is the class able to propate seeds through the algorithm? */
    virtual bool spCanEvaluate(bool fwd){ return false;}

    /** \brief  Reset the sparsity propagation */
    virtual void spInit(bool fwd){}
    
    /** \brief  Evaluate symbolically, SXElement type, possibly nonmatching sparsity patterns */
    virtual void evalSX(const std::vector<SX>& arg, std::vector<SX>& res, 
                        const std::vector<std::vector<SX> >& fseed, std::vector<std::vector<SX> >& fsens, 
                        const std::vector<std::vector<SX> >& aseed, std::vector<std::vector<SX> >& asens);

    /** \brief  Evaluate symbolically, MX type */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res, 
                        const std::vector<std::vector<MX> >& fseed, std::vector<std::vector<MX> >& fsens, 
                        const std::vector<std::vector<MX> >& aseed, std::vector<std::vector<MX> >& asens);

    /** \brief  Evaluate symbolically, SXElement type, matching sparsity patterns */
    virtual void evalSXsparse(const std::vector<SX>& arg, std::vector<SX>& res, 
                              const std::vector<std::vector<SX> >& fseed, std::vector<std::vector<SX> >& fsens, 
                              const std::vector<std::vector<SX> >& aseed, std::vector<std::vector<SX> >& asens);

    /** \brief  Create function call node */
    virtual void createCall(const std::vector<MX> &arg, std::vector<MX> &res, 
                            const std::vector<std::vector<MX> > &fseed, std::vector<std::vector<MX> > &fsens, 
                            const std::vector<std::vector<MX> > &aseed, std::vector<std::vector<MX> > &asens);

    /** \brief  Create derivative node */
    virtual void createCallDerivative(const std::vector<MX> &arg, std::vector<MX> &res, 
                                      const std::vector<std::vector<MX> > &fseed, std::vector<std::vector<MX> > &fsens, 
                                      const std::vector<std::vector<MX> > &aseed, std::vector<std::vector<MX> > &asens);

    /** \brief  Create a call to this */
    std::vector<MX> callSelf(const std::vector<MX> &arg);
    
    /** \brief Call a function, MX type (overloaded) */
    void call(const MXVector& arg, MXVector& res, 
              const MXVectorVector& fseed, MXVectorVector& fsens, 
              const MXVectorVector& aseed, MXVectorVector& asens,
              bool always_inline, bool never_inline);
    
    /** \brief Call a function, SXElement type (overloaded) */
    void call(const std::vector<SX>& arg, std::vector<SX>& res, 
              const std::vector<std::vector<SX> >& fseed, std::vector<std::vector<SX> >& fsens, 
              const std::vector<std::vector<SX> >& aseed, std::vector<std::vector<SX> >& asens,
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
    /** \brief Return tangent function */
    FX tangent(int iind, int oind);
    virtual FX getTangent(int iind, int oind);
    //@}

    //@{
    /** \brief Return Jacobian function */
    FX jacobian(int iind, int oind, bool compact, bool symmetric);
    void setJacobian(const FX& jac, int iind, int oind, bool compact);
    virtual FX getJacobian(int iind, int oind, bool compact, bool symmetric);
    virtual FX getNumericJacobian(int iind, int oind, bool compact, bool symmetric);
    //@}
    
    //@{
    /** \brief Return Jacobian of all input elements with respect to all output elements */
    FX fullJacobian();
    virtual FX getFullJacobian();
    //@}

    //@{
    /** \brief Return function that calculates forward derivatives 
     *    This method returns a cached instance if available, and calls FX getDerivative(int nfwd, int nadj) if no cached version is available.
     */
    FX derivative(int nfwd, int nadj);

    /** Set a function that calculates nfwd forward dedrivatives and nadj adjoint derivatives */
    void setDerivative(const FX& fcn, int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives */
    virtual FX getDerivative(int nfwd, int nadj);

    /** \brief Constructs and returns a function that calculates forward derivatives by creating the Jacobian then multiplying */
    virtual FX getDerivativeViaJac(int nfwd, int nadj);

    //@}
    
    
    /** \brief Create a helper MXFunction with some properties copied
    *
    * Copied properties: 
    *
    *    input/outputscheme
    *    ad_mode
    *   
    *  The function is not initialized
    */
    MXFunction wrapMXFunction();

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
    
    /** \brief Check if the numerical values of the supplied bounds make sense */
    virtual void checkInputs() const {};
            
    /** \brief Get the unidirectional or bidirectional partition */
    void getPartition(int iind, int oind, Sparsity& D1, Sparsity& D2, bool compact, bool symmetric);

    /// Verbose mode?
    bool verbose() const;
    
    /// Is function fcn being monitored
    bool monitored(const std::string& mod) const;
        
    /** \brief  Get total number of nonzeros in all of the matrix-valued inputs */
    int getNumInputNonzeros() const;

    /** \brief  Get total number of nonzeros in all of the matrix-valued outputs */
    int getNumOutputNonzeros() const;

    /** \brief  Get total number of elements in all of the matrix-valued inputs */
    int getNumInputElements() const;

    /** \brief  Get total number of elements in all of the matrix-valued outputs */
    int getNumOutputElements() const;
    
    /// Get all statistics obtained at the end of the last evaluate call
    const Dictionary & getStats() const;

    /// Get single statistic obtained at the end of the last evaluate call
    GenericType getStat(const std::string & name) const;
    
    /// Generate the sparsity of a Jacobian block
    virtual Sparsity getJacSparsity(int iind, int oind, bool symmetric);
    
    /// A flavour of getJacSparsity without any magic
    Sparsity getJacSparsityPlain(int iind, int oind);
    
    /// A flavour of getJacSparsity that does hierachical block structure recognition
    Sparsity getJacSparsityHierarchical(int iind, int oind);
    
    /// A flavour of getJacSparsity that does hierachical block structure recognition for symmetric jacobians
    Sparsity getJacSparsityHierarchicalSymm(int iind, int oind);
    
    /// Generate the sparsity of a Jacobian block
    void setJacSparsity(const Sparsity& sp, int iind, int oind, bool compact);
    
    /// Get, if necessary generate, the sparsity of a Jacobian block
    Sparsity& jacSparsity(int iind, int oind, bool compact, bool symmetric);
    
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<MX> symbolicInput() const;

    /// Get a vector of symbolic variables corresponding to the outputs
    virtual std::vector<MX> symbolicOutput(const std::vector<MX>& arg);
  
    /// Get a vector of symbolic variables with the same dimensions as the inputs
    virtual std::vector<SX> symbolicInputSX() const;
  
    // Workaround helper functions: assign nonzeros but ignore all -1
    static void assignIgnore(MX& y, const MX& x, const std::vector<int>& nz);
    static void assignIgnore(SX& y, const SX& x, const std::vector<int>& nz);

    //@{
    /** \brief Access input/output scheme */
    inline const IOScheme& inputScheme() const{ return input_.scheme;}
    inline const IOScheme& outputScheme() const{ return output_.scheme;}
    inline IOScheme& inputScheme(){ return input_.scheme;}
    inline IOScheme& outputScheme(){ return output_.scheme;}
    //@}

    //@{
    /// Input/output structures of the function */
    inline const IOSchemeVector<DMatrix>& input_struct() const{ return input_;}
    inline const IOSchemeVector<DMatrix>& output_struct() const{ return output_;}
    inline IOSchemeVector<DMatrix>& input_struct(){ return input_;}
    inline IOSchemeVector<DMatrix>& output_struct(){ return output_;}
    //@}

    //@{
    /// Input/output access without checking (faster, but unsafe)
    inline const Matrix<double>& inputNoCheck(int iind=0) const{ return inputS<false>(iind);}
    inline const Matrix<double>& outputNoCheck(int oind=0) const{ return outputS<false>(oind);}

    inline Matrix<double>& inputNoCheck(int iind=0){ return inputS<false>(iind);}
    inline Matrix<double>& outputNoCheck(int oind=0){ return outputS<false>(oind);}
    //@}

    /** \brief  Log the status of the solver */
    void log(const std::string& msg) const;

    /** \brief  Log the status of the solver, function given */
    void log(const std::string& fcn, const std::string& msg) const;

    // Codegen function
    FX dynamicCompilation(FX f, std::string fname, std::string fdescr, std::string compiler);

    // The following functions are called internally from EvaluateMX. For documentation, see the MXNode class
    //@{
    virtual void evaluateD(MXNode* node, const DMatrixPtrV& arg, DMatrixPtrV& res, std::vector<int>& itmp, std::vector<double>& rtmp);
    virtual void evaluateSX(MXNode* node, const SXPtrV& arg, SXPtrV& res, std::vector<int>& itmp, std::vector<SXElement>& rtmp);
    virtual void evaluateMX(MXNode* node, const MXPtrV& arg, MXPtrV& res, const MXPtrVV& fseed, MXPtrVV& fsens, const MXPtrVV& aseed, MXPtrVV& asens, bool output_given);
    virtual void propagateSparsity(MXNode* node, DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp, std::vector<double>& rtmp, bool fwd);
    virtual void nTmp(MXNode* node, size_t& ni, size_t& nr);
    virtual void generateOperation(const MXNode* node, std::ostream &stream, const std::vector<std::string>& arg, const std::vector<std::string>& res, CodeGenerator& gen) const;
    virtual void printPart(const MXNode* node, std::ostream &stream, int part) const;
    //@}

    /** \brief  Inputs of the function */
    IOSchemeVector<DMatrix> input_;

    /** \brief  Output of the function */
    IOSchemeVector<DMatrix> output_;

    /** \brief  Verbose -- for debugging purposes */
    bool verbose_;
    
    /// Set of module names which are extra monitored
    std::set<std::string> monitors_;
    
    /** \brief  Dictionary of statistics (resulting from evaluate) */
    Dictionary stats_;
    
    /** \brief  Flag to indicate wether statistics must be gathered */
    bool gather_stats_;

    /// Cache for functions to evaluate directional derivatives
    std::vector<std::vector<WeakRef> > derivative_fcn_;

    /// Cache for full Jacobian
    WeakRef full_jacobian_;

    /// Cache for sparsities of the Jacobian blocks
    Matrix<Sparsity> jac_sparsity_, jac_sparsity_compact_;

    /// Cache for Jacobians
    Matrix<WeakRef> jac_, jac_compact_;

    /// User-set field
    void* user_data_;
    
    bool monitor_inputs_, monitor_outputs_;
        
    /// Errors are thrown when NaN is produced
    bool regularity_check_;
    
    /// Errors are thrown if numerical values of inputs look bad
    bool inputs_check_;
    
  };


} // namespace CasADi


#endif // FX_INTERNAL_HPP
