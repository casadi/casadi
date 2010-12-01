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

#ifndef OCP_TOOLS_HPP
#define OCP_TOOLS_HPP

#include "casadi/sx/sx_matrix.hpp"

namespace OPTICON{
  using namespace CasADi;

  // Go through the equations of an ocp and make them explicit, if possible
/*  void makeExplicit(AcadoOCP &ocp);*/
  
  // Make all equations implicit
/*  void makeImplicit(AcadoOCP &ocp);*/
  
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

  // Type of collocation points
  enum CollocationPoints{LEGENDRE,RADAU};
  static double** collocation_points[] = {legendre_points,radau_points};

  // Get the coefficeints for the collocation and continuity equations
  void get_collocation_coeff(int K, std::vector<std::vector<double> >& C, std::vector<double>& D, CollocationPoints cp);

  // Collocate a variable (one variable per finite element)
  void collocate(const SXMatrix& var, std::vector< SXMatrix >& VAR, int N);

  // Collocate a variable (K variables per finite element)
  void collocate(const SXMatrix& var, std::vector< std::vector< SXMatrix > >& VAR, int N, int K);
  
  // Collocate a variable (K variables per finite element)
  void collocate_final(const SXMatrix& var, SXMatrix &VARF);


  
  
  
  
  
  
  
  
  
  
  
  
#if 0
/** \brief  Eliminate time derivatives from the dynamic equations and replace them by state derivatives */
/** \brief   void eliminateTimeDerivatives(OCP_old &ocp); */

/** \brief  Make smooth by eliminating switches */
  void makeSmooth(OCP_old &ocp);

/** \brief  Eliminate a lagrange objective term by adding an additional state */
  void eliminateLagrangeTerm(OCP_old &ocp);

/** \brief  Convert implicit ic's into explicit ic's (when possible) and write the remaining ic's in a compact form to allow for newton iterations */
  void sortInitialConditions(OCP_old &ocp, bool make_explicit=false);

/** \brief  Reformulate the differential equation */
  void makeImplicit(OCP_old &ocp); // Rewrites the differential equation in fully implicit form: 0 = f(xdot,x,p,u)
  void makeSemiExplicit(OCP_old &ocp); // Rewrites the differential equation in semi-explicit form: [xdot,0] = [f(x,y,p,u),g(x,y,p,u)]
  void makeExplicit(OCP_old &ocp); // Rewrites the differential equation in explicit form: xdot = f(x,p,u)
  bool isExplicit(const OCP_old &ocp); // checks if the differential equation is of explicit form

  /** 
    Parametrize the control u by introducing n_disc parameters for each control
    The argument met specifies the discretization method 
      (by piecewise constant discretization on a uniform grid)
    Returns the the discretized controls (which are parameters of the ocp)
  */
  Matrix parametrizeControls(OCP_old &ocp, const Matrix &u, int n_disc, const Matrix& met=Matrix());

/** \brief  Event iteration - see Principles of Object-Oriented Modeling and Simulation with Modelica 2.1, chapter 18 */
  void eventIteration(OCP_old& ocp, bool forward=true); // should not have any argument

/** \brief  Generate Lagrange discretization matrices */
  void generateLegendreMatrices(int n, Matrix &D, Matrix &w, Matrix &tau);

/** \brief  Generate Lagrange polynomials */
  void generateLagrangePolynomials(const Matrix &t, const Matrix &tau, Matrix &L);

/** \brief  Print to screen */
  std::ostream& operator<<(std::ostream &stream, const OCP_old& ocp);

/** \brief  Eliminate all dependent variables in the functions (this function should be made unecessary by means of a smarter ocp class) */
  void eliminateDependent(OCP_old& ocp);

#endif
  
} // namespace OPTICON

#endif // OCP_TOOLS_HPP

