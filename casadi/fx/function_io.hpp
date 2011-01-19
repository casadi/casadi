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

namespace CasADi{

/** \brief  Structure that contains the numerical values for the inputs or outputs of a function
    This class is bound to disappear in the future.
  \author Joel Andersson 
  \date 2010
*/
class FunctionIO{
  public:
    /// Constructor
    FunctionIO();

    /// Get a reference to the data
    Matrix<double>& get();
    
    /// Get a const reference to the data
    const Matrix<double>& get() const;

    /// Get a reference to the forward derivative data
    Matrix<double>& getFwd(int dir=0);
    
    /// Get a const reference to the forward derivative data
    const Matrix<double>& getFwd(int dir=0) const;

    /// Get a reference to the adjoint derivative data
    Matrix<double>& getAdj(int dir=0);
    
    /// Get a const reference to the adjoint derivative data
    const Matrix<double>& getAdj(int dir=0) const;

    /// Set the number of forward derivative directions
    void setNumFwdDir(int ndir);
    
    /// Set the number of adjoint derivative directions
    void setNumAdjDir(int ndir);

    /// Get the number of forward derivative directions    
    int numFwdDir() const;

    /// Get the number of adjoint derivative directions    
    int numAdjDir() const;
    
    /// (Re)allocate the data
    void init();
    
    /// The input is sparse
    void setSparse();
    
  protected:
    
    /// Input/output data
    Matrix<double> mat_;
    
    /// Forward derivative data
    std::vector< Matrix<double> > matF_;
    
    /// Adjoint derivative data
    std::vector< Matrix<double> > matA_;

    /// Is dense?
    bool dense_;
    
};

} // namespace CasADi

#endif // FUNCTION_IO_HPP
