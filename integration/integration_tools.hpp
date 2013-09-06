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

#ifndef INTEGRATION_TOOLS_HPP
#define INTEGRATION_TOOLS_HPP

#include <vector>
#include <algorithm>
#include "symbolic/casadi_exception.hpp"
#include "symbolic/options_functionality.hpp"
#include "symbolic/mx/mx.hpp"

namespace CasADi{

  /** \brief Obtain collocation points of specific order and scheme
  \param scheme  'radau' or 'legendre'
  **/
  std::vector<double> collocationPoints(int order, const std::string& scheme="radau");
  
  
  /** \brief Obtain collocation interpolating matrices
  \param tau_root  location of collocation points, as obtained from collocationPoints
  \param[out] C interpolating coefficients to obtain derivatives
      Length: order+1, order + 1
      
      dX/dt @collPoint(j) ~ Sum_i C[j][i]*X@collPoint(i)
      
  \param[out] D interpolating coefficients to obtain end state
      Length: order+1
  */
#ifndef SWIG
  void collocationInterpolators(const std::vector<double> & tau_root, std::vector< std::vector<double> > &C, std::vector< double > &D);
#else // SWIG
  void collocationInterpolators(const std::vector<double> & tau_root, std::vector< std::vector<double> > &OUTPUT, std::vector< double > &OUTPUT);
#endif // SWIG
  
#ifndef SWIG
  static double legendre_points1[] = { 0.00000000000000000000, 0.50000000000000000000 };
  static double legendre_points2[] = { 0.00000000000000000000, 0.21132486540518713447, 0.78867513459481286553 };
  static double legendre_points3[] = { 0.00000000000000000000, 0.11270166537925824235, 0.50000000000000000000, 0.88729833462074170214 };
  static double legendre_points4[] = { 0.00000000000000000000, 0.06943184420297354720, 0.33000947820757187134, 0.66999052179242823968, 0.93056815579702623076 };
  static double legendre_points5[] = { 0.00000000000000000000, 0.04691007703066807366, 0.23076534494715861268, 0.49999999999999994449, 0.76923465505284149835, 0.95308992296933192634 };
  static double legendre_points6[] = { 0.00000000000000000000, 0.03376524289842347537, 0.16939530676686742616, 0.38069040695840172805, 0.61930959304159849399, 0.83060469323313235179, 0.96623475710157580298 };
  static double legendre_points7[] = { 0.00000000000000000000, 0.02544604382862047931, 0.12923440720030288098, 0.29707742431130129690, 0.50000000000000000000, 0.70292257568869853657, 0.87076559279969734106, 0.97455395617137896558 };
  static double legendre_points8[] = { 0.00000000000000000000, 0.01985507175123157886, 0.10166676129318691357, 0.23723379504183561561, 0.40828267875217505445, 0.59171732124782483453, 0.76276620495816449541, 0.89833323870681347501, 0.98014492824876797705 };
  static double legendre_points9[] = { 0.00000000000000000000, 0.01591988024618706810, 0.08198444633668211523, 0.19331428364970504319, 0.33787328829809543107, 0.49999999999999988898, 0.66212671170190451342, 0.80668571635029517886, 0.91801555366331766272, 0.98408011975381259884 };
  static double* legendre_points[] =  { 0, legendre_points1, legendre_points2, legendre_points3, legendre_points4, legendre_points5, legendre_points6, legendre_points7, legendre_points8, legendre_points9};
  
  // Radau collocation points
  static double radau_points1[] = { 0.00000000000000000000, 1.00000000000000000000 };
  static double radau_points2[] = { 0.00000000000000000000, 0.33333333333333337034, 1.00000000000000000000 };
  static double radau_points3[] = { 0.00000000000000000000, 0.15505102572168222297, 0.64494897427831787695, 1.00000000000000000000 };
  static double radau_points4[] = { 0.00000000000000000000, 0.08858795951270420632, 0.40946686444073465694, 0.78765946176084700170, 1.00000000000000000000 };
  static double radau_points5[] = { 0.00000000000000000000, 0.05710419611451822419, 0.27684301363812369168, 0.58359043236891683382, 0.86024013565621926247, 1.00000000000000000000 };
  static double radau_points6[] = { 0.00000000000000000000, 0.03980985705146905529, 0.19801341787360787761, 0.43797481024738621480, 0.69546427335363603106, 0.90146491420117347282, 1.00000000000000000000 };
  static double radau_points7[] = { 0.00000000000000000000, 0.02931642715978521885, 0.14807859966848435640, 0.33698469028115418666, 0.55867151877155019069, 0.76923386203005450490, 0.92694567131974103802, 1.00000000000000000000 };
  static double radau_points8[] = { 0.00000000000000000000, 0.02247938643871305597, 0.11467905316090415413, 0.26578982278458951338, 0.45284637366944457959, 0.64737528288683043876, 0.81975930826310761113, 0.94373743946307731001, 1.00000000000000000000 };
  static double radau_points9[] = { 0.00000000000000000000, 0.01777991514736393386, 0.09132360789979432347, 0.21430847939563035798, 0.37193216458327238438, 0.54518668480342658000, 0.71317524285556954666, 0.85563374295785443735, 0.95536604471003006012, 1.00000000000000000000 };
  static double* radau_points[] =  { 0, radau_points1, radau_points2, radau_points3, radau_points4, radau_points5, radau_points6, radau_points7, radau_points8, radau_points9};

  static double** collocation_points[] = {legendre_points,radau_points};
#endif // SWIG 

  // Type of collocation points
  enum CollocationPoints{LEGENDRE,RADAU};
  
  /** \brief Construct an explicit Runge-Kutta integrator
  * \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  * \param tf    Integration end time
  * \param order Order of integration
  * \param ne    Number of times the RK primitive is repeated over the integration interval
  */
  FX explicitRK(FX& f, const MX &tf=1, int order=4, int ne = 1);
  
  /** \brief Construct an implicit Runge-Kutta integrator
  * \param f dynamical system
  * \copydoc scheme_DAEInput
  * \copydoc scheme_DAEOutput
  * \param tf    Integration end time
  * \param order Order of integration
  * \param scheme Collocation scheme, as excepted by collocationPoints function.
  * \param ne    Number of times the RK primitive is repeated over the integration interval
  */
  FX implicitRK(FX& f, implicitFunctionCreator impl, const Dictionary& dict = Dictionary(), const MX &tf=1, int order=4, const std::string& scheme="radau", int ne = 1);
    
} // namespace CasADi

#endif // INTEGRATION_TOOLS_HPP

