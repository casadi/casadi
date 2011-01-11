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
#include "../matrix/matrix.hpp"

// TODO: Remove this class and replace with functions for the sparsity pattern as well as a small structure to be placed inside the FX class

namespace CasADi{

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
    /*    void setSparsity(const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& elind);*/
            
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

    /** \brief  Dense? */
    bool dense_;
    
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

#ifndef SWIG
    /** \brief  Get the non-zero elements, array */
    void get(double* val, int dir=0, Sparsity sp=SPARSE) const;    

    /** \brief  Set the non-zero elements, array */
    void set(const double* val, int dir=0, Sparsity sp=SPARSE);
#endif
    
    const std::vector<int>& rowind() const;
    std::vector<int>& rowind();
    const std::vector<int>& col() const;
    std::vector<int>& col();
    int rowind(int i) const;
    int col(int el) const;

    Matrix<double>& mat(int dir=0);
    const Matrix<double>& mat(int dir=0) const;
    Matrix<double>& matF(int dir=0);
    const Matrix<double>& matF(int dir=0) const;
    Matrix<double>& matA(int dir=0);
    const Matrix<double>& matA(int dir=0) const;

    
  protected:

    /** \brief  Input/output data */
    Matrix<double> mat_;
    
    /** \brief  Forward sensitivity data */
    std::vector< Matrix<double> > matF_;
    
    /** \brief  Adjoint sensitivity data */
    std::vector< Matrix<double> > matA_;
            
};

} // namespace CasADi

#endif // FUNCTION_IO_HPP
