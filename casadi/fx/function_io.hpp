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

#ifndef FUNCTION_IO_HPP
#define FUNCTION_IO_HPP

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
    FunctionIO();

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
            
    /// Set the number of derivative directions
    void setNumFwdDir(int ndir);
    void setNumAdjDir(int ndir);

    /// Get the number of derivative directions    
    int numFwdDir() const;
    int numAdjDir() const;
    
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

    /** \brief  Access the non-zero entries */
    std::vector<double>& data(int dir=0);

    /** \brief  Const access the non-zero entries */
    const std::vector<double>& data(int dir=0) const;

    /** \brief  Set the non-zero elements, scalar */
    void set(double val, int dir=0, Sparsity sp=SPARSE);
    
    /** \brief  Get the non-zero elements, scalar */
    void get(double& val, int dir=0, Sparsity sp=SPARSE) const;

    /** \brief  Set the non-zero elements, vector */
    void set(const std::vector<double>& val, int dir=0, Sparsity sp=SPARSE);

    /** \brief  Get the non-zero elements, vector */
    void get(std::vector<double>& val, int dir=0, Sparsity sp=SPARSE) const;

    /** \brief  Get the non-zero elements, array */
    void get(double* val, int dir=0, Sparsity sp=SPARSE) const;    

#ifndef SWIG
    /** \brief  Set the non-zero elements, array */
    void set(const double* val, int dir=0, Sparsity sp=SPARSE);
#endif
    
    /** \brief  Get the result */
    void getSparseSym(double *res, int dir=0) const;    // general sparse, symmetric matrix

    /** \brief  Get the result times a vector */
    void getTimesVector(const double *v, double *res, int dir=0) const;

    /** \brief  Save the result to the LAPACK banded format -- see LAPACK documentation 
    kl:    The number of subdiagonals in res 
    ku:    The number of superdiagonals in res 
    ldres: The leading dimension in res 
    res:   The number of superdiagonals */
    void getBand(int kl, int ku, int ldres, double *res, int dir=0) const;
        
  protected:

    /** \brief  Size */
    int nfdir_, nadir_;
    
    /** \brief  Non-zero elements for the input/output or forward/adjoint seeds/derivatives */
    std::vector< std::vector<double> > data_;

    /** \brief Assert that the number of nonzeros is correct */
    void assertNNZ(int sz, Sparsity sp) const;
    
    /** \brief Assert that the number of elements is correct */
    void assertNumEl(int sz) const;
        
};

} // namespace CasADi

#endif // FUNCTION_IO_HPP
