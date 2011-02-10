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
  \author Joel Andersson 
  \date 2010-2011
*/
class FunctionIO{
  public:
    /// Constructor
    FunctionIO();

    /// (Re)allocate the data
    void init();
    
    /// Input/output data
    Matrix<double> data;
    
    /// Forward derivative data
    std::vector< Matrix<double> > dataF;
    
    /// Adjoint derivative data
    std::vector< Matrix<double> > dataA;
};

} // namespace CasADi

#endif // FUNCTION_IO_HPP
