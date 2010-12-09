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

#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>

#define PI 3.141592653589793238462643383279

using namespace std;
using namespace KINEMATICS;

// CART PENDULUM
//
// This is a trivial example to get stated with the kinematics toolbox.
// pendulum on a cart

int main(int argc, char *argv[]){

try{
  
SXMatrix t("t");			// Time
SXMatrix x("x"),theta("theta");	// The states we will use
SXMatrix dx("dx"),dtheta("dtheta");

SXMatrix L("L");			// Length of the cart

SXMatrix q; q<< x << theta;
SXMatrix dq; dq<< dx << dtheta;

Frame f0("world frame",q,dq,t); // The world frame is an inertial frame
				// It has a notion of time.

// We chain transformations together as follows:

Frame f1("cart frame",f0,tr(x,0,0)); // The frame attached to the cart
Frame f2("pendulum anchor frame",f1,TRzp(-1)*TRz(theta)); // The frame attached to the pendulum's anchorpoint
Frame f3("pendulum CM frame",f2,tr(L/2,0,0)); // The frame attached to the pendulum's CM



cout << "We express all in the inertial frame" << endl;

cout << "== positions ==" << endl;
KinVec v=pos(f1,f0);
cout << "The position of the cart is: " << v.getCoords() << endl;
v=pos(f3,f0);
cout << "The position of the pendulum's CM is: " << v.getCoords() << endl;

cout << "== velocities ==" << endl;
v=vel(f1,f0,f0);
cout << "The velocity of the cart is: " << v.getCoords() << endl;
v=vel(f2,f0,f0);
cout << "The velocity of the pendulum's anchorpoint is: " << v.getCoords() << endl;
v=vel(f3,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;

cout << "== accelerations ==" << endl;
v=acc(f3,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;

cout << "== rotations ==" << endl;
v=rotVel(f2,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;
v=rotVel(f3,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;

cout << "== rotations ==" << endl;
v=rotAcc(f2,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;
v=rotAcc(f3,f0,f0);
cout << "The velocity of the pendulum's CM is: " << v.getCoords() << endl;


} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
