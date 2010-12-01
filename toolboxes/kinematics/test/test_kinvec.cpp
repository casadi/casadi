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

#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>

using namespace KINEMATICS;
using namespace std;

int main(int argc, char *argv[]){

try{

  SXMatrix alpha("alpha");
  SXMatrix L("L");
  SXMatrix H("H");
  SXMatrix dalpha("dalpha");
  SXMatrix dL("dL");
  SXMatrix ddalpha("ddalpha");
  SXMatrix ddL("ddL");
  
  SXMatrix q(2,1);q[0]=L;q[1]=alpha;
  SXMatrix dq(2,1);dq[0]=dL;dq[1]=dalpha;
  SXMatrix ddq(2,1);ddq[0]=ddL;ddq[1]=ddalpha;

  Frame f0("World Frame",q,dq,ddq);
  Frame f1("f1",f0,TRz(alpha)*tr(0,0,H));
  Frame fa("fa",f1,tr(L,0,0));
  
 
  cout << "Test case #1 - pure kinvec transformations" << endl << endl;
  KinVec p(1,2,3,0,f1);
  cout << p.expressedIn(f0) << p << p.expressedIn(fa) << endl;
  p=KinVec(1,2,3,0,f0);
  cout << p << p.expressedIn(f1) << p.expressedIn(fa) << endl;
  p=KinVec(1,2,3,0,fa);
  
  cout << p.expressedIn(f0) << p.expressedIn(f1) <<  p << endl;
  
  cout << "Test case #2 - pos constructor" << endl << endl;
    
  p=pos(fa,f0);  
  cout << p << p.expressedIn(f1) << p.expressedIn(fa) << endl;
  p=pos(f0,fa);  
  cout << p << p.expressedIn(f1) << p.expressedIn(f0) << endl;
  
  cout << "Test case #3 - vel constructor" << endl << endl; 
  
  KinVec v=vel(fa,f0,f0);  
  cout << v << endl;

  v=vel(fa,f0,f1);  
  cout << v << endl;
  
  v=vel(fa,f1,f0);  
  cout << v << endl;
  
  v=vel(fa,f1,f1);  
  cout << v << endl;
  
  cout << "Test case #4 - accel constructor" << endl << endl;

  
  KinVec a=acc(fa,f1,f1);  
  cout << a << endl;
  
  a=acc(fa,f0,f0);  
  cout << a << endl;
  
  
  cout << "Test case #5 - omega constructor" << endl << endl;
  
  KinVec w=rotVel(fa,f0,f0);  
  cout << w << endl;
  
  w=rotVel(fa,f0,fa);   
  cout << w << endl;
  
    
  cout << "Test case #6 - alpha constructor" << endl << endl;
  
  KinVec al=rotAcc(fa,f0,f0);  
  cout << al << endl;
  
  al=rotAcc(fa,f0,fa);   
  cout << al << endl;
  
  cout << "Test case #7 - Operations on KinVecs" << endl << endl;
  
  KinVec n(1,2,3,0,f0);
  a=vel(fa,f1,f1);  
  KinVec b=vel(fa,f1,fa); 
  cout << a;
  cout << b;
  cout << a+b;
  KinVec na=n-a;
  cout << na << endl;
  SXMatrix I=eye(3);I(2,2)=5;
  cout << I*na;
    
  cout << "Test case #8 - Operations on KinVecs" << endl << endl;
  v=vel(fa,f1,f0); 
  a=KinVec(0,0,0,0,f0); 
  cout << v;
  cout << a;
  cout << v+a;
  cout << norm(a);
  
  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
