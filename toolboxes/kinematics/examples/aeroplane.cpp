#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>



using namespace std;
//using namespace KinVec;

int main(int argc, char *argv[]){

  
// ------------------------------
// Symbolic variable definitions|
// ------------------------------

Expression	t("t");
Expression      x("x"), y("y"), z("z"), R("R"), P("P"), Y("Y");
Expression      dx("dx"), dy("dy"), dz("dz"), dR("dR"), dP("dP"), dY("dY");
  
Expression      alpha("alpha"),beta("beta");

Expression q;   q   << x  << y  << z  << R  << P  <<  Y;
Expression dq;  dq  <<dx  << dy << dz << dR << dP << dY;
Expression ddq; ddq <<ddx <<ddy <<ddz <<ddR <<ddP <<ddY;

Expression      mu("mu"),nu("nu");
Expression	dmu("dmu"),dnu("dnu");  
  
// -------------------
// Frame definitions |
// -------------------

Frame f0("Airport frame",q,dq,t);
Frame f1("plane CM frame",f0,tr(x,y,z));
Frame f2("proper plane CM frame",f1,TRperm(1,-2,-3));
Frame f3("plane body frame",f2,TRz(Y)*TRy(P)*TRx(R));
Frame f4("wind frame",f3,TRz(beta)*TRy(-alpha));

// -------------------
// Model parameters  |
// -------------------

const double m             =       0.8; //APS.defparameter(plane/inertia/mk,m)
const double A_w           = 0.1012901; //APS.defparameter(plane/wing/A)
const double b_w           =      1.00; //APS.defparameter(plane/wing/b)
const double c_w           = 0.1012901; //APS.defparameter(plane/wing/c)
const double Ixx           =      0.08; //APS.defparameter(plane/inertia/Ixx)
const double Iyy           =      0.05; //APS.defparameter(plane/inertia/Iyy)
const double Izz           =      0.15; //APS.defparameter(plane/inertia/Izz)
const double g             =      9.81; //APS.defparameter(universal/g)
const double rho           =      1.23; //APS.defparameter(universal/rho)
const double PI            = 3.1415926; //APS.defparameter(universal/PI)

Expression I(3,3);I(0,0)=Ixx;I(1,1)=Iyy;I(2,2)=Izz;

// FCS.load(parameters),FCS.settype(cppjoel)
// Matrix A in the expression (y=A.x + b), which gives the force coefficients y in terms of inputs x   [FCS.descr] FCS.declare(Av,CA) 
Matrix CA(5,6);  // FCS.auto 
CA(0,0)=+5.393613e+00;CA(0,1)=+0.000000e+00;CA(0,2)=+0.000000e+00;CA(0,3)=+2.857637e-01;CA(0,4)=+0.000000e+00;CA(0,5)=+0.000000e+00;  // FCS.auto 
CA(1,0)=+5.065224e-02;CA(1,1)=+0.000000e+00;CA(1,2)=+0.000000e+00;CA(1,3)=+0.000000e+00;CA(1,4)=+1.013045e+00;CA(1,5)=+0.000000e+00;  // FCS.auto 
CA(2,0)=+0.000000e+00;CA(2,1)=-1.570796e-02;CA(2,2)=+2.215820e-01;CA(2,3)=+0.000000e+00;CA(2,4)=+0.000000e+00;CA(2,5)=+0.000000e+00;  // FCS.auto 
CA(3,0)=-6.104714e-01;CA(3,1)=+0.000000e+00;CA(3,2)=+0.000000e+00;CA(3,3)=-1.213133e+00;CA(3,4)=+0.000000e+00;CA(3,5)=+0.000000e+00;  // FCS.auto 
CA(4,0)=+0.000000e+00;CA(4,1)=+1.608902e-01;CA(4,2)=+0.000000e+00;CA(4,3)=+0.000000e+00;CA(4,4)=+0.000000e+00;CA(4,5)=-3.760209e-03; // FCS.auto


// Matrix b in the expression (y=A.x + b), which gives the force coefficients y in terms of inputs x   [FCS.descr] FCS.declare(bv,Cb) 
Matrix Cb(5,1);  // FCS.auto 
Cb(0,0)=+1.253403e-01;  // FCS.auto 
Cb(1,0)=+2.063315e-02;  // FCS.auto 
Cb(2,0)=+0.000000e+00;  // FCS.auto 
Cb(3,0)=+2.250681e-01;  // FCS.auto 
Cb(4,0)=+0.000000e+00; // FCS.auto

double Rdamp = 1;
double Pdamp = 1;
double Ydamp = 1;

// -------------------
// Forces and moments|
// -------------------

KinVec We(0,0,0,0,f0); //external wind velocity    [ m/s   ]
KinVec W=We-vel(f3,f0,f3); // Wind in frame 3 coordinates 
  
beta =  atan(W[1]/W[0]);
alpha=  -atan(W[2]/sqrt(W[1]*W[1]+W[0]*W[0]));

Expression Cv=CA*ax+Cb;

Expression CL=Cv[0];Expression CD=Cv[1];
Expression CR=Cv[2]-dR*Rdamp;
Expression CP=Cv[3]-dP*Pdamp;
Expression CY=Cv[4]-dY*Ydamp;

KinVec Fg=-m_*g*ez(f0);
KinVec Faer=-CL*A_w*qaer*ez(f4)+CD*A_w*qaer*ex(f4);
KinVec F = Fg + Faer; // total force

KinVec T(CR*qaer*A_w*b_w,CP*qaer*A_w*c_w,CY*qaer*A_w*b_w,0,f3);

// ---------------------
// Equations of motion |
// ---------------------

F==m*acc(f3,f0,f0);

KinVec w=omega(f3,f0,f3);
KinVec a=omegad(f3,f0,f3);

T == (cross(w,I*w) + I*a);



}
