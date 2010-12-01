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

#ifndef MODEL_SIMULATOR
#define MODEL_SIMULATOR

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


#include <math.h>

namespace KINEMATICS{
using namespace CasADi;
using namespace boost::numeric;


class  ModelSimulator
{
public:
  ModelSimulator();
  ~ModelSimulator();
  
  ublas::vector<double> evaluate(double t,const ublas::vector<double> &q,const ublas::vector<double> &u);
  ublas::vector<double> evaluatelin(double t,const ublas::vector<double> &q);
  ublas::vector<double> output();

  void integrateOnce(double dt);
  void integrate(double dt);
  void integratelin(double dt);
  /** \brief initialize the state */
  void init(const ublas::vector<double> &q);
  void setInput(ublas::vector<double> input);
  void setdeltat(double t);
  void linearize(double t,const ublas::vector<double> &q,const ublas::vector<double> &u);

void initDumb(int i,double ini_);
double outputDumb(int i);
void setInputDumb(int i,double ini_);

void test();

private:
  SXFunction f;
  SXFunction g;
  SXFunction j;
  ublas::vector<double> state;
  ublas::vector<double> input;
  double t;
  double deltat; /** smallest timestep for integration */
  double Deltat; /** smallest timestep for linearisation */
  ublas::matrix<double> A;
  ublas::vector<double> b;
  ublas::vector<double> q0;
  ublas::vector<double> u0;

  ublas::vector<double> k1;
  ublas::vector<double> k2;
  ublas::vector<double> k3;
  ublas::vector<double> k4;
  ublas::vector<double> state_;

};

} // namespace KINEMATICS

#endif //MODEL_SIMULATOR
