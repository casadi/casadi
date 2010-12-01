#include <iostream>
#include <fstream>
#include <string>

#include <kinetics.hpp>

using namespace std;
using namespace KINEMATICS;

SXMatrix ff1(const SXMatrix &a, const SXMatrix &b) {
  // FCS.load(pointmodel),FCS.settype(cppjoel)
  // cable drag formula f1: int x*sqrt(a+x*b+x^2*c)   [FCS.descr] FCS.declare(ff1,C) 
  SXMatrix C=(3*b*(pow(b,2)-4*a)*log(2*(pow(b+a+1,0.5)+1)+b)+2*pow(b+a+1,0.5)*(-3*pow(b,2)+2*b+8*a+8))/48-(3*b*(pow(b,2)-4*a)*log(b+2*pow(a,0.5))+2*pow(a,0.5)*(8*a-3*pow(b,2)))/48; // FCS.auto
  return C;
}

SXMatrix ff2(const SXMatrix  &a, const SXMatrix  &b) {
  // FCS.load(pointmodel),FCS.settype(cppjoel)
  // cable drag formula f2: int x^2*sqrt(a+x*b+x^2*c)   [FCS.descr] FCS.declare(ff2,C) 
  SXMatrix C=(2*pow(b+a+1,0.5)*(15*pow(b,3)-10*pow(b,2)-52*a*b+8*b+24*a+48)-3*(5*pow(b,4)-24*a*pow(b,2)+16*pow(a,2))*log(2*(pow(b+a+1,0.5)+1)+b))/384-(2*pow(a,0.5)*(15*pow(b,3)-52*a*b)-3*(5*pow(b,4)-24*a*pow(b,2)+16*pow(a,2))*log(b+2*pow(a,0.5)))/384; // FCS.auto
  return C;
}

void splitdep(const KinVec &k,  int i,
                       KinVec           &v1,
                       KinVec           &v2  ){
    

    SXMatrix x=zeros(k.q.size(),1);
    SXMatrix dx=zeros(k.dq.size(),1);
    x[i]=k.q[i];
    dx[i]=k.dq[i];
	
    SXMatrix t  = jacobian(k.v,k.q)*x+jacobian(k.v,k.dq)*dx;  // NO !  SHOULD BE FORWARD DERIVATIVE WITH SEED x or dx !!!
    SXMatrix t2 = k.v-t;

    v2 = KinVec( t , k.type, k.ref, k.q, k.dq, k.ddq, k.order );
    v1 = KinVec( t2, k.type, k.ref, k.q, k.dq, k.ddq, k.order );
}

void numeval(const KinVec &v, const SXMatrix &q, const SXMatrix &dq, const SXMatrix &ddq) {
  SXMatrix arg;arg << q << dq;

  
  SXFunction temp(arg,v.getCoords());

  //double x0[17]={0,1,0.1,0.2,0.3,0.4,0.5,10,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,0,0.1,0.05};
  double x0[20]={0,1,0.1,0.2,0.3,0.4,0.5,0,0.1,0.05,10,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,0,0,0};
  vector<double> res(3);
  for(int ii=0; ii<temp->input[0].data().size(); ++ii)
    temp->input[0].data()[ii] = x0[ii];

  temp.evaluate();
  temp.getOutputDense(&res[0]);
  cout << "Num eval ["<< res[0] << "," << res[1] << "," << res[2] << "]"<< endl;
}

void numeval(const SXMatrix &v, const SXMatrix &q, const SXMatrix &dq, const SXMatrix &ddq) {
  SXMatrix arg;arg << q << dq;
  
  SXFunction temp(arg,v);

  double x0[20]={0,1,0.1,0.2,0.3,0.4,0.5,0,0.1,0.05,10,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,0,0,0};
  double res;
  for(int ii=0; ii<temp->input[0].data().size(); ++ii)
    temp->input[0].data()[ii] = x0[ii];

  temp.evaluate();
  temp.getOutputDense(&res);
  cout << "Num eval "<< res << endl;
}

int main(int argc, char *argv[]){

try{
  
  // <<DEFINE STATES
  SXMatrix	  t("t");
  SXMatrix      delta("delta"), r("r"),  phi("phi"), theta("theta"), R("R"), P("P"), Y("Y");
  SXMatrix      ddelta("ddelta"), dr("dr"),  dphi("dphi"), dtheta("dtheta"), dR("dR"), dP("dP"), dY("dY");
  SXMatrix      E("E"),dE("dE"),ddE("ddE");
  SXMatrix      mu("mu"),nu("nu");
  SXMatrix	  dmu("dmu"),dnu("dnu");
  SXMatrix	  ddmu("ddmu"),ddnu("ddnu");
  SXMatrix      alpha("alpha"),beta("beta");
  
  SXMatrix      dddelta("dddelta"), ddr("ddr"),  ddphi("ddphi"), ddtheta("ddtheta"), ddR("ddR"), ddP("ddP"), ddY("ddY");
    
  SXMatrix q; q << delta << r <<phi <<theta << R << P <<Y << E << mu << nu;
  SXMatrix dq; dq << ddelta << dr <<dphi <<dtheta << dR << dP <<dY << dE << dmu << dnu;
  
  dddelta=0;
  ddr=0;
  
  SXMatrix ddq; ddq << dddelta << ddr <<ddphi <<ddtheta << ddR << ddP <<ddY << ddE << ddmu << ddnu ;
  
  SXMatrix ddqe=zeros(ddq.size());
  ddqe[0]=dddelta;
  ddqe[1]=ddr;
  
  Frame f0("world frame",q,dq,ddq);
  
  // >>
  

  
  // <<DEFINE PARAMETERS
  // APS.load(simple)
  const double L               =       1.1; //APS.defparameter(carrousel/dimensions/L)
  const double H               =       2.5; //APS.defparameter(carrousel/dimensions/H)
  const double mk              =       0.5; //APS.defparameter(plane/inertia/mk)
  const double A_w             = 0.1012901; //APS.defparameter(plane/wing/A)
  const double b_w             =      1.00; //APS.defparameter(plane/wing/b)
  const double c_w             = 0.1012901; //APS.defparameter(plane/wing/c)
  const double x_a             =         0; //APS.defparameter(plane/anchorpoint/x)
  const double y_a             =         0; //APS.defparameter(plane/anchorpoint/y)
  const double z_a             =     0.025; //APS.defparameter(plane/anchorpoint/z)
  const double Ixx             =     0.008; //APS.defparameter(plane/inertia/Ixx)
  const double Iyy             =     0.010; //APS.defparameter(plane/inertia/Iyy)
  const double Izz             =     0.016; //APS.defparameter(plane/inertia/Izz)
  const double Rd               =     0.005; //APS.defparameter(plane/Rd)
  const double Pd               =     0.005; //APS.defparameter(plane/Pd)
  const double Yd               =     0.005; //APS.defparameter(plane/Yd)
  const double g               =      9.81; //APS.defparameter(universal/g)
  const double rho             =      1.23; //APS.defparameter(universal/rho)
  const double rho_c           =      1450; //APS.defparameter(cable/rho)
  const double c_c             =         1; //APS.defparameter(cable/c)
  const double d_c             =     0.003; //APS.defparameter(cable/d)
  const double PI              = 3.1415926; //APS.defparameter(universal/PI)
  const double mu_c            = 0.0102494; //APS.defparameter(cable/mu)
  
  SXMatrix mc      =  mu_c*r        ;   // mass of the cable
  
  SXMatrix m_      =  mk + mc     / 2.0;   // effective gravitational mass  
    
  SXMatrix dmc = mu_c*dr;   // time derivative of the mass

  // FCS.load(parameters),FCS.settype(cppjoel)
  // SXMatrix A in the expression (y=A.x + b), which gives the force coefficients y in terms of inputs x   [FCS.descr] FCS.declare(Av,CA) 
  SXMatrix CA(5,6);  // FCS.auto 
  CA(0,0)=+5.393613e+00;CA(0,1)=+0.000000e+00;CA(0,2)=+0.000000e+00;CA(0,3)=+2.857637e-01;CA(0,4)=+0.000000e+00;CA(0,5)=+0.000000e+00;  // FCS.auto 
  CA(1,0)=+5.065224e-02;CA(1,1)=+0.000000e+00;CA(1,2)=+0.000000e+00;CA(1,3)=+0.000000e+00;CA(1,4)=+1.013045e+00;CA(1,5)=+0.000000e+00;  // FCS.auto 
  CA(2,0)=+0.000000e+00;CA(2,1)=-1.570796e-02;CA(2,2)=+2.215820e-01;CA(2,3)=+0.000000e+00;CA(2,4)=+0.000000e+00;CA(2,5)=+0.000000e+00;  // FCS.auto 
  CA(3,0)=-6.104714e-01;CA(3,1)=+0.000000e+00;CA(3,2)=+0.000000e+00;CA(3,3)=-1.213133e+00;CA(3,4)=+0.000000e+00;CA(3,5)=+0.000000e+00;  // FCS.auto 
  CA(4,0)=+0.000000e+00;CA(4,1)=+1.608902e-01;CA(4,2)=+0.000000e+00;CA(4,3)=+0.000000e+00;CA(4,4)=+0.000000e+00;CA(4,5)=-3.760209e-03; // FCS.auto
  
  
  // SXMatrix b in the expression (y=A.x + b), which gives the force coefficients y in terms of inputs x   [FCS.descr] FCS.declare(bv,Cb) 
  SXMatrix Cb(5,1);  // FCS.auto 
  Cb(0,0)=+1.253403e-01;  // FCS.auto 
  Cb(1,0)=+2.063315e-02;  // FCS.auto 
  Cb(2,0)=+0.000000e+00;  // FCS.auto 
  Cb(3,0)=+2.250681e-01;  // FCS.auto 
  Cb(4,0)=+0.000000e+00; // FCS.auto
      
  KinVec We(0,0,0,0,f0);//external wind velocity    [ m/s   ]

  SXMatrix I(3,3);I(0,0)=Ixx;I(1,1)=Iyy;I(2,2)=Izz;
  //>>

  // << DEFINE FRAMES

  Frame f1("anchorpoint frame",f0,TRz(delta)*tr(L,0,H));
  Frame f2("spherical frame",f1,TRz(phi)*TRy(-theta)*tr(r,0,0));
  Frame f3("ideal-pose frame",f2,TRperm(2,-3,-1));
  Frame f4("body frame",f3,TRz(Y)*TRy(P)*TRx(R));

  
  // >>
   
  //<<CALCULATE FORCES AND MOMENTS

  KinVec W=We.expressedIn(f4)-vel(f4,f0,f4); // Wind in frame 4 coordinates
  
  numeval(W,q,dq,ddq);
  
  SXMatrix normW=norm(W);KinVec nW=W/normW;
 
  KinVec Wv1,Wv2;splitdep(W.expressedIn(f3),1,Wv1,Wv2); //r dependent and independant part
  Wv1[2]=0;Wv2[2]=0;
  
  SXMatrix normWv2=norm(Wv2);
           
  SXMatrix Wa=pow(norm(Wv1),2)/pow(normWv2,2);
  SXMatrix Wb=2*Wv1*Wv2/pow(normWv2,2);

  beta =  atan(W[1]/W[0]);
  alpha=  -atan(W[2]/sqrt(W[1]*W[1]+W[0]*W[0]));
  
  Frame f5("wind frame",f4,TRz(beta)*TRy(-alpha));

  SXMatrix Cf      =  (rho*d_c/2.0)*r;
  SXMatrix qaer    =  rho*normW*normW/2;
  
  SXMatrix ax;ax << alpha << beta << mu << nu << alpha*alpha << alpha*beta;
  SXMatrix Cv=CA*ax+Cb;

  SXMatrix CL=Cv[0];SXMatrix CD=Cv[1];
  SXMatrix CR=Cv[2]-dR*Rd;
  SXMatrix CP=Cv[3]-dP*Pd;
  SXMatrix CY=Cv[4]-dY*Yd;
  

 
  KinVec Fg=-m_*g*ez(f0,f3);
  KinVec Faer=-CL*A_w*qaer*ez(f5,f3)-CD*A_w*qaer*ex(f5,f3);
  KinVec Ff=Cf  *  c_c*normWv2 * (ff1(Wa,Wb)* Wv1 + ff2(Wa,Wb)* Wv2);
  KinVec F = Fg + Faer + Ff; // total force
  
  

	  
  cout << "forces" << endl;
        numeval(Faer,q,dq,ddq);

  KinVec T(CR*qaer*A_w*b_w,CP*qaer*A_w*c_w,CY*qaer*A_w*b_w,0,f4);	// total torque
  

  //>>
  
  numeval(CR,q,dq,ddq);
  
  numeval(qaer,q,dq,ddq);
 numeval(A_w,q,dq,ddq);
  numeval(b_w,q,dq,ddq);

  KinVec v=vel(f3,f0,f0);
  KinVec v1,v2;
  splitdep(v,1,v1,v2);
  
  KinVec  p  = mk * v + mc* (v1/2+v2/3); // Combined linear impulse of kite and cable
  KinVec v0i =vel(f1,f0,f3);v0i[1]=0;v0i[2]=0;
  
  KinVec eq1=p.der().expressedIn(f3)-(F+dmc*v0i);
  
  		// Now we are going to solve for ddphi and ddtheta
		// We use lines 1 and 2 of eq1 for this
		
		std::vector<int> di;
		di.push_back(2);di.push_back(3); // use ddphi and ddtheta
		std::vector<int> ri;
		
		ri.push_back(0);ri.push_back(1); // use lines 1 and 2 of eq1

		SXMatrix rpt=eq1.explicitize(di,ri); // find ddphi and ddtheta

		ddqe[2]=rpt[0];                         // fill it in ddq vector
		ddqe[3]=rpt[1];
	
		// Now we are ready to find Fc
		SXMatrix Fc=getRow(eq1.J,2)*ddqe + getRow(eq1.c,2);
		
  numeval(rpt[0],q,dq,ddq);
    numeval(rpt[1],q,dq,ddq);
      numeval(Fc,q,dq,ddq);

  T+=cross(KinVec(x_a,y_a,z_a,1,f4),-Fc*ex(f2,f4));


  KinVec w=rotVel(f4,f0,f4);
  KinVec a1=rotAcc(f4,f0,f4);a1.setDDQ(ddqe);
  KinVec a=rotAcc(f4,f0,f4);

  KinVec eq2=T - (cross(w,I*w) + I*a); // Euler equation
  eq2.setDDQ(ddqe);

	  // Now we are going to solve for ddR and ddP ddY
	  // We use all of eq2 for this
	  std::vector<int> dic;
	  dic.push_back(4);// use ddR and ddP and ddY
	  dic.push_back(5);
	  dic.push_back(6);

	  SXMatrix rRPY=eq2.explicitize(dic); // find ddR and ddP and ddY
	  ddqe[4]=rRPY[0];                      
	  ddqe[5]=rRPY[1];
	  ddqe[6]=rRPY[2];
  
  numeval(rRPY[0],q,dq,ddq);
    numeval(rRPY[1],q,dq,ddq);
      numeval(rRPY[2],q,dq,ddq);
      
  //cout << F[0]; better not print this


  return 0;

} catch (const char * str){
  cerr << str << endl;
  return 1;
}
}
