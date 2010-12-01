#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>

using namespace KINEMATICS;
using namespace std;

void utest(const std::string& name,const KinVec &v, const SXMatrix &q, const SXMatrix &dq, const SXMatrix &ddq) {
  SXMatrix arg;arg << q << dq << ddq;
  SXFunction temp(arg,v.getCoords());

  double x0[18]={0.6784619950185145,0.0363933158999219,0.6252975183418072,.9937314578694818,.04158924615444981,.6326253893666012,
		.5931521815169765,.6546577661145074,.4707106740336644,.05192355366451751,.002327880375313951,.3156076962859116,
		.8602295262448114,.04206547654085679,.8486699428493936,.2329775607444364,.6911141946221595,0.6115558081936561};
  vector<double> res(3);
  for(int ii=0; ii<temp->input[0].data().size(); ++ii)
    temp->input[0].data()[ii] = x0[ii];
  temp.evaluate();
  temp.getOutputDense(&res[0]);
  cout << "Unit test " << name << ": ["<< res[0] << "," << res[1] << "," << res[2] << "]"<< endl;
}

int main(int argc, char *argv[]){

try{


  SX a("a"),b("b"),c("c");
  SX x("x"),y("y"),z("z");
  
  
  SXMatrix u=14;


  SX da("da"),db("db"),dc("dc");
  SX dx("dx"),dy("dy"),dz("dz");

  SX dda("dda"),ddb("ddb"),ddc("ddc");
  SX ddx("ddx"),ddy("ddy"),ddz("ddz");

  SXMatrix q;q << a << b << c << x << y << z;
  SXMatrix dq;dq << da << db << dc << dx << dy << dz;
  SXMatrix ddq;ddq << dda << ddb << ddc << ddx << ddy << ddz;
   

  Frame f0("World Frame",q,dq,ddq);
  Frame f1("f1",f0,TRx(a)*tr(x,0,0)*TRy(b));
  Frame f2("f2",f1,tr(0,y,z));
  Frame f3("f3",f2,TRz(c));
  Frame f4("f4",f3,TRz(u));
  
	
  //utest("q0",KinVec(a,b,c,0,f0),q,dq,t);utest("q1",KinVec(x,y,z,0,f0),q,dq,t);
  //utest("dq0",KinVec(da,db,dc,0,f0),q,dq,t);utest("dq1",KinVec(dx,dy,dz,0,f0),q,dq,t);
  //utest("ddq0",KinVec(da.der(t),db.der(t),dc.der(t),0,f0),q,dq,t);utest("ddq1",KinVec(dx.der(t),dy.der(t),dz.der(t),0,f0),q,dq,t);

	KinVec v = pos( f3, f0 );

		
	utest("p303",v.expressedIn(f3),q,dq,ddq);


	
	v = vel( f3, f0, f0, q, dq, ddq );

	utest("v300",v,q,dq,ddq);

	utest("v301",vel(f3,f0,f1),q,dq,ddq);

	utest("a100",acc(f1,f0,f0),q,dq,ddq);
	//cout << acc(f1,f0,f0) << endl;
	utest("a300",acc(f3,f0,f0),q,dq,ddq);
	utest("a301",acc(f3,f0,f1),q,dq,ddq);




	KinVec w=rotVel(f1,f0,f0);
	utest("w100",w,q,dq,ddq);
	utest("w300",rotVel(f3,f0,f0),q,dq,ddq);
	utest("w103",rotVel(f1,f0,f3),q,dq,ddq);

	utest("alph103",rotAcc(f1,f0,f3),q,dq,ddq);


   // Operations on KinVecs

   KinVec A=rotVel(f1,f0,f0)+rotVel(f2,f1,f0)+rotVel(f3,f2,f0);
   KinVec B=rotVel(f3,f0,f0); // a-c is zero (relation from mechanics)
   KinVec C=rotVel(f3,f0,f1); // b-c is zero

   utest("A",A,q,dq,ddq);
   utest("0",A-B,q,dq,ddq);
   utest("0",A-C,q,dq,ddq);

   utest("0",cross(A,B*2),q,dq,ddq); // parallel vectors have zero crossproduct

   KinVec x0(1,0,0,0,f0);
   KinVec x1(0,1,0,0,f0);

   utest("0",cross(x0,x1)-ez(f0),q,dq,ddq); // (1,0,0) ^ (0,1,0) = (0,0,1)

   utest("A",A,q,dq,ddq);
   A[2]=3.14;
   utest("A3.14",A,q,dq,ddq);
  // raw constructors
  KinVec x2(a,b,c,0,f0);
  KinVec x3(da,db,dc,0,f0);
  //KinVec x4=x2.der();
  utest("0",x2.der()-x3,q,dq,ddq);
  
  // Mixed Datatypes
  SXMatrix	    c0=a;
  SXMatrix c1=b;
  KinVec x4(c0,c1,c,0,f0);
  utest("0",x4.der()-x3,q,dq,ddq);

  KinVec a0=acc(f4,f0,f0);


   //

	SXMatrix phi("phi"),theta("theta"),delta("delta");	// The states we will use
	SXMatrix dphi("dphi"),dtheta("dtheta"),ddelta("ddelta");
	SXMatrix ddphi("ddphi"),ddtheta("ddtheta"),dddelta("dddelta");

	SXMatrix q_; q_<< phi << theta << delta;      // The time changing variables of our system
	SXMatrix dq_; dq_<< dphi << dtheta << ddelta; // The derivatives of these variables
	SXMatrix ddq_; ddq_<< ddphi << ddtheta << dddelta;

	SXMatrix r("r");

	Frame F0("world frame",q_,dq_,ddq_);
	Frame F1("CM frame",F0,TRz(phi)*TRy(-theta)*tr(a,0,0)); 
	Frame F2("rotating frame",F1,TRx(delta)); 

	
  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
