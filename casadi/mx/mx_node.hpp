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

#ifndef MX_NODE_HPP
#define MX_NODE_HPP

#include "mx.hpp"
#include "../fx/fx.hpp"
#include <vector>

namespace CasADi{
  //@{
  /** \brief Convenience function, convert vectors to vectors of pointers */
  template<class T>
  std::vector<T*> ptrVec(std::vector<T>& v){
    std::vector<T*> ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = &v[i];
    return ret;
  }

  template<class T>
  const std::vector<T*> ptrVec(const std::vector<T>& v){
    std::vector<T*> ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = const_cast<T*>(&v[i]);
    return ret;
  }
  
  template<class T>
  std::vector<std::vector<T*> > ptrVec(std::vector<std::vector<T> >& v){
    std::vector<std::vector<T*> > ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = ptrVec(v[i]);
    return ret;
  }
  
  template<class T>
  const std::vector<std::vector<T*> > ptrVec(const std::vector<std::vector<T> >& v){
    std::vector<std::vector<T*> > ret(v.size());
    for(int i=0; i<v.size(); ++i) 
      ret[i] = ptrVec(v[i]);
    return ret;
  }
  //@}

  
  /** \brief Node class for MX objects
    \author Joel Andersson 
    \date 2010
    Internal class.
*/
class MXNode : public SharedObjectNode{
  friend class MX;
  friend class MXFunctionInternal;
  
  public:
    /// Constructor
    MXNode();
  
    /** \brief  Destructor */
    virtual ~MXNode()=0;

    /** \brief  Clone function */
    virtual MXNode* clone() const = 0;

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
    
    /** \brief  Print expression */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const=0;
    
    /** \brief  Print expression */
    virtual void print(std::ostream &stream) const;

    /** \brief  Evaluate the function */
    virtual void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output, 
                          const DMatrixPtrVV& fwdSeed, DMatrixPtrVV& fwdSens, 
                          const DMatrixPtrVV& adjSeed, DMatrixPtrVV& adjSens) = 0;

    /** \brief  Evaluate the function, no derivatives*/
    void evaluate(const DMatrixPtrV& input, DMatrixPtrV& output);

    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output, 
                            const SXMatrixPtrVV& fwdSeed, SXMatrixPtrVV& fwdSens, 
                            const SXMatrixPtrVV& adjSeed, SXMatrixPtrVV& adjSens) = 0;

    /** \brief  Evaluate symbolically (SX), no derivatives */
    void evaluateSX(const SXMatrixPtrV& input, SXMatrixPtrV& output);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, 
                            const MXPtrVV& fwdSeed, MXPtrVV& fwdSens, 
                            const MXPtrVV& adjSeed, MXPtrVV& adjSens, bool output_given)=0;

    /** \brief  Evaluate symbolically (MX), no derivatives */
    void evaluateMX(const MXPtrV& input, MXPtrV& output);
    
    /** \brief  Get the name */
    virtual const std::string& getName() const;
    
    /** \brief  Check if symbolic */
    virtual bool isSymbolic() const;

    /** \brief  Check if constant */
    virtual bool isConstant() const;

    /** \brief  Check if mapping */
    virtual bool isMapping() const{return false;}

    /** \brief  Check if evaluation */
    virtual bool isEvaluation() const{return false;}

    /** \brief  Check if evaluation output */
    virtual bool isOutputNode() const{return false;}

    /** \brief  Check if jacobian reference */
    virtual bool isJacobian() const{return false;}

    /** \brief  Get function reference */
    virtual FX& getFunction();

    /** \brief  Get function input */
    virtual int getFunctionInput() const;

    /** \brief  Get function output */
    virtual int getFunctionOutput() const;

    /** \brief  dependencies - functions that have to be evaluated before this one */
    const MX& dep(int ind=0) const;
    MX& dep(int ind=0);
    
    /** \brief  Number of dependencies */
    int ndep() const;
    
    /** \brief  Does the node depend on other nodes*/
    virtual bool hasDep() const{return ndep()>0; }
    
    /** \brief  Number of outputs */
    virtual int getNumOutputs() const{ return 1;}

    /// Get the sparsity
    const CRSSparsity& sparsity() const;
    
    /** \brief Is the node nonlinear */
    virtual bool isNonLinear(){return false;}
    
    /// Set the sparsity
    void setSparsity(const CRSSparsity& sparsity);
    
    /// Set unary dependency
    void setDependencies(const MX& dep);
    
    /// Set binary dependencies
    void setDependencies(const MX& dep1, const MX& dep2);
    
    /// Set ternary dependencies
    void setDependencies(const MX& dep1, const MX& dep2, const MX& dep3);
    
    /// Set multiple dependencies
    void setDependencies(const std::vector<MX>& dep);
        
    /// Add a dependency
    int addDependency(const MX& dep);

    /// Add a dependency (index given)
    virtual void addDependency(int depind, const std::vector<int>& nz_d, const std::vector<int>& nz);

    /// Add a dependency (mapping matrix)
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d, const std::vector<int>& nz);
    
    /// Add a dependency (mapping matrix)
    virtual void addDependency(const MX& d, const std::vector<int>& nz_d);

    /// Is it a certain operation
    virtual bool isOperation(int op) const{ return false;}
    
    /** \brief  Get the jacobian of an function evaluation with respect to the iind-th argument */
    virtual MX jac(int iind);
    
    /// Number of elements
    int numel() const;
    
    /// Get size
    int size() const;
    
    /// Get size
    int size1() const;
    
    /// Get size
    int size2() const;
    
    /** Temporary variables to be used in user algorithms like sorting, 
    the user is resposible of making sure that use is thread-safe
    The variable is initialized to zero
    */
    int temp;
    
    /// Numeric evaluation
    virtual Matrix<double> eval(const std::vector<DMatrix>& x){
      std::vector<DMatrix> ret(1,DMatrix(sparsity_));
      const DMatrixPtrV mx_input = ptrVec(x);
      DMatrixPtrV mx_output = ptrVec(ret);
      evaluate(mx_input,mx_output);
      return ret[0];
    }

    /// Symbolic evaluation (scalar graph)
    virtual Matrix<SX> eval(const std::vector<Matrix<SX> >& x){return Matrix<SX>();}

    /// Symbolic evaluation (matrix graph)
    virtual MX eval(const std::vector<MX>& x){return MX();}
    
    /** \brief  dependencies - functions that have to be evaluated before this one */
    std::vector<MX> dep_;
    
    /** \brief  The sparsity pattern */
    CRSSparsity sparsity_;
};


} // namespace CasADi


#endif // MX_NODE_HPP
