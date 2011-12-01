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

                        
} // namespace CasADi


#endif // FX_TOOLS_HPP
