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

#ifndef OCP_TOOLS_HPP
#define OCP_TOOLS_HPP

#include "symbolic_ocp.hpp"

namespace CasADi{
  
  /// Update dependent variables in an OCP
  void updateDependent(SymbolicOCP& ocp);

  // Type of collocation points
  enum CollocationPoints{LEGENDRE,RADAU};


#ifndef SWIG
  // Legendre collocation points
  static double legendre_points1[] = {0,0.500000};
  static double legendre_points2[] = {0,0.211325,0.788675};
  static double legendre_points3[] = {0,0.112702,0.500000,0.887298};
  static double legendre_points4[] = {0,0.069432,0.330009,0.669991,0.930568};
  static double legendre_points5[] = {0,0.046910,0.230765,0.500000,0.769235,0.953090};
  static double* legendre_points[] = {0,legendre_points1,legendre_points2,legendre_points3,legendre_points4,legendre_points5};

  // Radau collocation points
  static double radau_points1[] = {0,1.000000};
  static double radau_points2[] = {0,0.333333,1.000000};
  static double radau_points3[] = {0,0.155051,0.644949,1.000000};
  static double radau_points4[] = {0,0.088588,0.409467,0.787659,1.000000};
  static double radau_points5[] = {0,0.057104,0.276843,0.583590,0.860240,1.000000};
  static double* radau_points[] = {0,radau_points1,radau_points2,radau_points3,radau_points4,radau_points5};

  static double** collocation_points[] = {legendre_points,radau_points};
#endif // SWIG 
  
} // namespace CasADi

#endif // OCP_TOOLS_HPP

