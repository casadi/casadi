/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef FUNCTION_IO_HPP
#define FUNCTION_IO_HPP

#include "../matrix_size.hpp"
#include <vector>

// TODO: Remove this class and replace with functions for the sparsity pattern as well as a small structure to be placed inside the FX class

namespace CasADi{

/** Sparsity format for getting and setting inputs and outputs */
enum Sparsity{SPARSE,SPARSESYM,DENSE,DENSESYM};
  
/** \brief  An output or input of the function (forward or adjoint)
  \author Joel Andersson 
  \date 2010
  This class is bound to disappear in the future.
*/
class FunctionIO{
  public:
/** \brief  Constructor */
    explicit FunctionIO();

    /** \brief  Set size */
    void setSize(int nrow, int ncol=1);

    /** \brief  Get sparsity pattern, CRS format */
    void getSparsityCRS(std::vector<int>& rowind, std::vector<int> &col) const; // Compressed row storage

    /** \brief  Get sparsity pattern, general sparse format */
    void getSparsity(std::vector<int>& row, std::vector<int> &col) const; // General sparse
    
    /** \brief  Set the sparsity pattern, CRS format */
    void setSparsityCRS(const std::vector<int>& rowind, const std::vector<int> &col); // Compressed row storage

    /** \brief  Set the sparsity pattern, general sparse format */
    void setSparsity(const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& elind);

    /** \brief  Set the non-zero elements */
    template <class T> void set(const T& val, Sparsity sp=SPARSE){setv(val,data(),sp);}

    /** \brief  Set the adjoint sensitivity non-zero elements */
    template <class T> void setA(const T& val, int dir=0, Sparsity sp=SPARSE) {setv(val,dataA(dir),sp);}

    /** \brief  Set the forward sensitivity non-zero elements */
    template <class T> void setF(const T& val, int dir=0, Sparsity sp=SPARSE) {setv(val,dataF(dir),sp);}

    
#if 0
// THIS IMPLEMENTATION SEEMS TO BE INCORRECT FOR SCALAR ARGUMENTS! WORKAROUND BELOW

    /** \brief  Get the non-zero elements */
    template <class T> void get(T& val) const{getv(val,data());}
    
    /** \brief  Get the forward sensitivity non-zero elements */
    template <class T> void getF(T& val, int dir=0) const {getv(val,dataF(dir));}

    /** \brief  Get the adjoint sensitivity non-zero elements */
    template <class T> void getA(T& val, int dir=0) const {getv(val,dataA(dir));}
#else
// WORKAROUND
#define GETTERS(T)\
    void get(T val, Sparsity sp=SPARSE) const{getv(val,data(),sp);} \
    void getF(T val, int dir=0, Sparsity sp=SPARSE) const {getv(val,dataF(dir),sp);} \
    void getA(T val, int dir=0, Sparsity sp=SPARSE) const {getv(val,dataA(dir),sp);}

GETTERS(double&)
GETTERS(double*)
GETTERS(std::vector<double>&)
#undef GETTERS
#endif

    
    /** \brief  Get the result */
    void getDense(double *x, int ord=0) const;
    void getSparseSym(double *res, int ord=0) const;    // general sparse, symmetric matrix
    void getTimesVector(const double *v, double *res, int ord=0) const;

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation */
    /** \brief  kl:  The number of subdiagonals in res */
    /** \brief  ku:  The number of superdiagonals in res */
    /** \brief  ldres:  The leading dimension in res */
    /** \brief  res:  The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, double *res, int ord=0) const;
        
    /// Set the number of derivative directions
    void setNumFwdDir(int ndir);
    void setNumAdjDir(int ndir);

    /// Get the number of derivative directions    
    int numFwdDir() const;
    int numAdjDir() const;
    
  /** \brief  Access the non-zero entries */
    std::vector<double>& data();
    const std::vector<double>& data() const;
    
  /** \brief  Access the non-zero entries of the forward sensitivities */
    std::vector<double>& dataF(int dir=0);
    const std::vector<double>& dataF(int dir=0) const;

  /** \brief  Access the non-zero entries of the adjoint sensitivities */
    std::vector<double>& dataA(int dir=0);
    const std::vector<double>& dataA(int dir=0) const;
        
  /** \brief  Get the size */
    int size1() const;
    int size2() const;

    /** \brief  Number of non-zeros */
    int size() const;

    /** \brief  Number of non-zeros in the upper triangular half */
    int sizeU() const;

    /** \brief  Number of non-zeros in the lower triangular half */
    int sizeL() const;

    /** \brief  Get the number of elements */
    int numel() const;
    
    // (Re)allocate the data
    void init();

    /** \brief  Sparsity (CRS format): index of the first element on each row */
    std::vector<int> rowind_;

    /** \brief  Sparsity (CRS format): column for each non-zero element */
    std::vector<int> col_;
    
    /** \brief  Size */
    int nrow_, ncol_;

  protected:
    
    /** \brief  Non-zero elements for the input/output */
    std::vector<double> data_;

    /** \brief  Non-zero elements for the first order forward seeds/derivatives */
    std::vector< std::vector<double> > dataF_;

    /** \brief  Non-zero elements for the first order adjoint seeds/derivatives */
    std::vector< std::vector<double> > dataA_;

    /** \brief Assert that the number of nonzeros is correct */
    void assertNNZ(int sz) const;
    
    /** \brief Assert that the number of elements is correct */
    void assertNumEl(int sz) const;
    
    // Legacy
    std::vector<double>& data(int ord);
    const std::vector<double>& data(int ord) const;

    /** \brief  Set the non-zero elements, scalar */
    void setv(double val, std::vector<double>& v, Sparsity sp) const;
    
    /** \brief  Get the non-zero elements, scalar */
    void getv(double& val, const std::vector<double>& v, Sparsity sp) const;

    /** \brief  Set the non-zero elements, vector */
    void setv(const std::vector<double>& val, std::vector<double>& v, Sparsity sp) const;

    /** \brief  Get the non-zero elements, vector */
    void getv(std::vector<double>& val, const std::vector<double>& v, Sparsity sp) const;

    /** \brief  Set the non-zero elements, array */
    void setv(const double* val, std::vector<double>& v, Sparsity sp) const;

    /** \brief  Get the non-zero elements, array */
    void getv(double* val, const std::vector<double>& v, Sparsity sp) const;

    
    
};

} // namespace CasADi

#endif // FUNCTION_IO_HPP
