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

#ifndef FX_TOOLS_HPP
#define FX_TOOLS_HPP

#include "fx.hpp"
#include "mx_function.hpp"

namespace CasADi{

    /** \brief Prints out a human readable report about possible constraint violations - specific constraints
    *
    * Constraint visualiser strip:
    *  o-------=-------o   Indicates that the value is nicely inbetween the bounds
    *  o-=-------------o   Indicates that the value is closer to the lower bound
    *  X---------------o   Indicates that the lower bound is active
    *  8---------------o   Indicates that the lower bound is -inifinity
    *  o------------=--o   Indicates that the value is closer to the upper bound
    *  o---------------X   Indicates that the upper bound is active
    *  o---------------8   Indicates that the upper bound is inifinity
    *     VIOLATED         Indicates constraint violation
    *
    */
    void reportConstraints(std::ostream &stream, const Matrix<double> &v, const Matrix<double> &lb, const Matrix<double> &ub, const std::string &name, double tol=1e-8);

    /** \brief make integration start and end time a parameter
    Applies the conversion  t = t0 + (tf-t0)*tau to the supplied dae.
    with tau dimensionless time.
    The Input/OuputScheme of the result is the same as the scheme of the dae, except for input(DAE_P), which is extended by t0 and tf at the top.
    */
    FX parameterizeTime(FX dae);
    
    /** \brief adapts an output function such that start and end time are parameters
    Applies the conversion  t = t0 + (tf-t0)*tau to the supplied dae.
    with tau dimensionless time.
    The InputScheme of the result is the same as the scheme of the dae, except for input(DAE_P), which is extended by t0 and tf at the top.
    */
    FX parameterizeTimeOutput(FX outputfcn);
    
    /** \brief sample a function on a 1D grid
    * \param fx an initialized function mapping from single p-by-1 to single m-by-n
    * \param a grid of numbers p-by-N
    *
    *  For each row in the grid, fx is numerically evaluated and the output is put in a resulting matrix of size m-by-n*p
    *
    * If your fx is really multiple output, and you wish to use a particular output, use the slice operator on the fx.
    *
    * @see evalf
    */
    Matrix<double> numSample1D(FX &fx, const Matrix<double> &grid);
    
    /** \brief sample a function on a 1D grid
    * \param fx an initialized function mapping from single p-by-1 to single m-by-n
    * \param a grid of numbers p-by-N
    *
    *  For each row in the grid, fx is numerically evaluated and the output is put in a resulting matrix of size m*p-by-n
    *
    * If your fx is really multiple output, and you wish to use a particular output, use the slice operator on the fx.
    */
    Matrix<double> numSample1DT(FX &fx, const Matrix<double> &grid);
    
    /** \brief sample a function on a 2D grid
    * \param fx an initialized function mapping from (1-by-1,1-by-1) to single m-by-n
    * \param a grid of numbers
    *
    *  For each point (i,j) in the grid, fx is numerically evaluated and the output is put in a matrix
    *
    *  If your fx is really multiple output, and you wish to use a particular output, use the slice operator on the fx.
    */
    Matrix<double> numSample2D(FX &fx, const Matrix<double> &grid);
                        
} // namespace CasADi


#endif // FX_TOOLS_HPP
