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

#ifndef SX_FUNCTION_HPP
#define SX_FUNCTION_HPP

#include "fx.hpp"

namespace CasADi{

#ifndef SWIG

  /** \brief  An atomic operation for the SX virtual machine */
  struct ScalarAtomic{
    int op;     /// Operator index
    int i0;
    union{
      double d;
      struct{ int i1,i2; };
    };
  };

#endif // SWIG

  /// Forward declaration of internal class
  class SXFunctionInternal;

  /// Forward declaration of MXFunction
  class MXFunction;

  /**   \brief Dynamically created function that can be expanded into a series of scalar operations.
        \author Joel Andersson 
        \date 2010-2013
  */

  class SXFunction : public FX{

  public:
    /// Default constructor
    SXFunction();
  
    /// Expand an MXFunction
    explicit SXFunction(const MXFunction &f);

    /// Expand an FX
    explicit SXFunction(const FX &f);
  
    /// Multiple (matrix valued) input, multiple (matrix valued) output 
    SXFunction(const std::vector< SXMatrix>& arg, const std::vector<SXMatrix>& res);

    /// Multiple (matrix valued) input, multiple (matrix valued) output 
    SXFunction(const std::vector< SXMatrix>& arg, const IOSchemeVector< SXMatrix >& res);

    /// Multiple (matrix valued) input, multiple (matrix valued) output 
    SXFunction(const IOSchemeVector< SXMatrix >& arg, const std::vector< SXMatrix>& res);
  
    /// Multiple (matrix valued) input, multiple (matrix valued) output 
    SXFunction(const IOSchemeVector< SXMatrix >& arg, const IOSchemeVector< SXMatrix >& res);
  
#ifndef SWIG

    /// Multiple (vector valued) input, multiple (vector valued) output 
    SXFunction(const std::vector< std::vector<SX> >& arg, const std::vector< std::vector<SX> >& res);

    /// Single (scalar/matrix/vector valued) input, single (scalar/matrix/vector valued) output  
    SXFunction(const SXMatrix& arg, const SXMatrix& res);

    /// Multiple (vector valued) input, single (scalar/vector/matrix valued) output 
    SXFunction(const std::vector< std::vector<SX> >& arg, const SXMatrix& res);

    /// Multiple (matrix valued) input, single (scalar/vector/matrix valued) output 
    SXFunction(const std::vector< SXMatrix>& arg, const SXMatrix& res);

    /// Single (scalar/vector/matrix valued) input, multiple (vector valued) output 
    SXFunction(const SXMatrix& arg, const std::vector< std::vector<SX> >& res);

    /// Single (scalar/vector/matrix valued) input, multiple (matrix valued) output 
    SXFunction(const SXMatrix& arg, const std::vector< SXMatrix>& res);
#endif // SWIG

    /// Access functions of the node 
    SXFunctionInternal* operator->();

    /// Const access functions of the node 
    const SXFunctionInternal* operator->() const;

    //@{
    /** \brief Jacobian via source code transformation
     *
     * \see CasADi::Jacobian for an AD approach
     */
    SXMatrix jac(int iind=0, int oind=0, bool compact=false, bool symmetric=false);
    SXMatrix jac(const std::string& iname, int oind=0, bool compact=false, bool symmetric=false) { return jac(inputSchemeEntry(iname),oind,compact,symmetric); } 
    SXMatrix jac(int iind, const std::string& oname, bool compact=false, bool symmetric=false) { return jac(iind,outputSchemeEntry(oname),compact,symmetric); } 
    SXMatrix jac(const std::string& iname, const std::string& oname, bool compact=false, bool symmetric=false) { return jac(inputSchemeEntry(iname),outputSchemeEntry(oname),compact,symmetric); } 
    //@}
   
    //@{
    /// Gradient via source code transformation
    SXMatrix grad(int iind=0, int oind=0);
    SXMatrix grad(const std::string& iname, int oind=0) { return grad(inputSchemeEntry(iname),oind); }
    SXMatrix grad(int iind, const std::string& oname) { return grad(iind,outputSchemeEntry(oname)); }
    SXMatrix grad(const std::string& iname, const std::string& oname) { return grad(inputSchemeEntry(iname),outputSchemeEntry(oname)); }
    //@}
  
    //@{
    /// Hessian (forward over adjoint) via source code transformation
    SXMatrix hess(int iind=0, int oind=0);
    SXMatrix hess(const std::string& iname, int oind=0) { return hess(inputSchemeEntry(iname),oind); }
    SXMatrix hess(int iind, const std::string& oname) { return hess(iind,outputSchemeEntry(oname)); }
    SXMatrix hess(const std::string& iname, const std::string& oname) { return hess(inputSchemeEntry(iname),outputSchemeEntry(oname)); }
    //@}
  
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
    
    /** \brief Get function input */
    const SXMatrix& inputExpr(int iind) const;
    const SXMatrix& inputExpr(const std::string& iname) const { return inputExpr(inputSchemeEntry(iname)); }
  
    /** \brief Get function output */
    const SXMatrix& outputExpr(int oind) const;
    const SXMatrix& outputExpr(const std::string& oname) const { return outputExpr(outputSchemeEntry(oname)); }
  
    /** \brief Get all function inputs */
    const std::vector<SXMatrix>& inputExpr() const;
  
    /** \brief Get all function outputs */
    const std::vector<SXMatrix> & outputExpr() const;
    
#ifndef SWIG
    /** \brief Access the algorithm directly */
    const std::vector<ScalarAtomic>& algorithm() const;
#endif // SWIG
  
    /** \brief Get the number of atomic operations */
    int getAlgorithmSize() const{ return algorithm().size();}

    /** \brief Get the length of the work vector */
    int getWorkSize() const;

    /** \brief Get an atomic operation operator index */
    int getAtomicOperation(int k) const{ return algorithm().at(k).op;}

    /** \brief Get the (integer) input arguments of an atomic operation */
    std::pair<int,int> getAtomicInput(int k) const{ 
      const ScalarAtomic& atomic = algorithm().at(k);
      return std::pair<int,int>(atomic.i1,atomic.i2);}

    /** \brief Get the floating point output argument of an atomic operation */
    double getAtomicInputReal(int k) const{ return algorithm().at(k).d;}

    /** \brief Get the (integer) output argument of an atomic operation */
    int getAtomicOutput(int k) const{ return algorithm().at(k).i0;}

    /** \brief Number of nodes in the algorithm */
    int countNodes() const;
  
    /** \brief Clear the function from its symbolic representation, to free up memory, no symbolic evaluations are possible after this */
    void clearSymbolic();
 
    /** \brief Get all the free variables of the function */
    std::vector<SX> getFree() const;

    /** \brief Get the corresponding matrix type */
    typedef SXMatrix MatType;  
  
#ifndef SWIG 
    /// Construct a function that has only the k'th output
    SXFunction operator[](int k) const;
#endif //SWIG 

    SXFunction indexed_one_based(int k) const{ return operator[](k-1);}
    SXFunction indexed_zero_based(int k) const{ return operator[](k);}

  };

} // namespace CasADi

#endif // SX_FUNCTION_HPP
