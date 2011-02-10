#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/fx/mx_function.hpp>
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "interfaces/sundials/cvodes_integrator.hpp"
#include "interfaces/sundials/idas_integrator.hpp"
#include "casadi/fx/jacobian.hpp"
#include "casadi/fx/parallelizer.hpp"
#include "casadi/fx/c_function.hpp"

using namespace CasADi;
using namespace CasADi::Sundials;
using namespace std;

//Final time (fixed)
double tf = 20.0;

//Control bounds
double u_min = -0.75;
double u_max = 1.0;
double u_init = 0.0;

//State bounds and initial guess
double inf = numeric_limits<double>::infinity();
double x_min = -inf, x_max =  inf, x_init = 0;
double y_min = -inf,  y_max =  inf,  y_init = 0;
double L_min = -inf, L_max =  inf, L_init = 0;

//State bounds at the initial time
double x0_min = 0,        x0_max =  0;
double y0_min = 1,        y0_max =  1;
double L0_min = 0,        L0_max =  0;

//State bounds at the final time
double xf_min = 0,        xf_max =  0;
double yf_min = 0,        yf_max =  0;
double Lf_min = -inf,     Lf_max =  inf;

//Numboer of shooting nodes
int NS = 50;

class Jac{
  public:
    
    Jac(){
        // Declare variables (use simple, efficient DAG)
      SX t("t"); //time
      SX x("x"), y("y"), u("u"), L("cost");

      // All states
      vector<SX> xx(3);  xx[0] = x;  xx[1] = y;  xx[2] = L;

      //ODE right hand side
      vector<SX> f(3);
      f[0] = (1 - y*y)*x - y + u;
      f[1] = x;
      f[2] = x*x + y*y + u*u;
      
      //ODE right hand side function
      vector<SXMatrix> rhs_in(ODE_NUM_IN);
      rhs_in[ODE_T] = t;
      rhs_in[ODE_Y] = xx;
      rhs_in[ODE_P] = u;
      SXFunction rhs(rhs_in,f);

      //Create an integrator (CVodes)
      I = CVodesIntegrator(rhs);
      I.setOption("abstol",1e-8); //abs. tolerance
      I.setOption("reltol",1e-8); //rel. tolerance
      I.setOption("steps_per_checkpoint",500);
      I.setOption("stop_at_end",true);
      I.init();

      FX IjacX = Jacobian(I,INTEGRATOR_X0,INTEGRATOR_XF);
      IjacX.setOption("ad_mode","forward");
      PX = Parallelizer(vector<FX>(NS,IjacX));
      
      IjacP = Jacobian(I,INTEGRATOR_P,INTEGRATOR_XF);
      IjacP.setOption("ad_mode","forward");
      PP = Parallelizer(vector<FX>(NS,IjacP));
      
      PX.init();
      PP.init();

      //Number of discretized controls
      int NU = NS;

      //Number of discretized states
      int NX = 3*(NS+1);

      //Declare variable vector
      int NV = NU+NX;
      MX V("V",NV);

      //Disretized control
      vector<MX> U;
      
      //Disretized state
      vector<MX> X;

      for(int i=0; i<NS; ++i){
        X.push_back(V(range(i*4,i*4+3),range(1)));
        UX_init.push_back(x_init);
        UX_init.push_back(y_init);
        UX_init.push_back(L_init);
        
        if(i==0){
          UX_min.push_back(x0_min);       UX_min.push_back(y0_min);       UX_min.push_back(L0_min);
          UX_max.push_back(x0_max);       UX_max.push_back(y0_max);       UX_max.push_back(L0_max);
        } else {
          UX_min.push_back(x_min);        UX_min.push_back(y_min);        UX_min.push_back(L_min);
          UX_max.push_back(x_max);        UX_max.push_back(y_max);        UX_max.push_back(L_max);
        }
        
        U.push_back(V[i*4+3]);
        UX_min.push_back(u_min);
        UX_max.push_back(u_max);
        UX_init.push_back(u_init);
      }
      
      X.push_back(V(range(NS*4,NS*4+3),range(1)));
      
      UX_init.push_back(x_init);
      UX_init.push_back(y_init);
      UX_init.push_back(L_init);

      UX_min.push_back(xf_min);
      UX_min.push_back(yf_min);
      UX_min.push_back(Lf_min);
      
      UX_max.push_back(xf_max);
      UX_max.push_back(yf_max);
      UX_max.push_back(Lf_max);

        
      //Beginning of each shooting interval (shift time horizon)
      vector<MX> T0(NS,MX(0));

      //End of each shooting interval (shift time horizon)
      vector<MX> TF(NS,MX(tf/NS));

      //The initial state (x=0, y=1, L=0)

      //State derivative (not used)
      MX XP;

      //Algebraic state (not used)
      MX Z;

      MX XF;

      
      //Constraint function with upper and lower bounds
      vector<MX> g;
      vector<SXMatrix> C;
      SXMatrix JJ(3*NS,3*NS+NS+3);

      //Build up a graph of integrator calls
      for(int k=0; k<NS; ++k){
        MX I_in[] = {T0[k],TF[k],X[k],U[k],XP,Z};
        
        //call the integrator
        vector<MX> Ires = I.call(vector<MX>(I_in,I_in+6));
        XF = Ires[0];
        XP = Ires[1];
        Z = Ires[2];
        
        //append continuity constraints
        g.push_back(XF-X[k+1]);
        for(int j=0; j<XF.numel(); ++j){
          g_min.push_back(0);
          g_max.push_back(0);
        }
        
        
        // Create a block
        stringstream ss;
        ss << k;
        
        SXMatrix CX = symbolic(string("X_")+ss.str(),3,3);
        SXMatrix CP = symbolic(string("P_")+ss.str(),3,1);
        C.push_back(CX);
        C.push_back(CP);
        for(int ii=0; ii<3; ++ii){
          for(int jj=0; jj<3; ++jj){
            JJ(3*k+ii, 4*k+jj) = SX(CX(ii,jj)); // syntax workaround
          }
          
          for(int jj=0; jj<1; ++jj){
            JJ(3*k+ii, 4*k+3+jj) = SX(CP(ii,jj)); // error!
          }
          
          JJ(3*k+ii, 4*k+4+ii) = -1;
        }
      }
      
/*      cout << JJ << endl;*/
      
      J3 = SXFunction(C,JJ);
      J3.init();

      //State at the final time
      XF = X[NS-1];

      //Objective function: L(T)
      F = MXFunction(V,XF[2]);

      //Terminal constraints: 0<=[x(T);y(T)]<=0
      G = MXFunction(V,vertcat(g));
      
      J = Jacobian(G);
      J.setOption("ad_mode","forward");
      J.init();
      
      J2 = CFunction(Jac::jacobian_wrapper);
      J2.setUserData(this);
      J2.setNumInputs(J3.getNumInputs());
      for(int i=0; i<J2.getNumInputs(); ++i){
        J2.input(i) = J3.input(i);
      }
      J2.setNumOutputs(J3.getNumOutputs());
      for(int i=0; i<J2.getNumOutputs(); ++i){
        J2.output(i) = J3.output(i);
      }
      
    }
    
    static void jacobian_wrapper(CFunction &f, int fsens_order, int asens_order, void* user_data){  
      casadi_assert(fsens_order==0 && asens_order==0);
      Jac* J = (Jac*)user_data;
      J->jacobian(f,fsens_order,asens_order);
    }

    void jacobian(CFunction &f, int fsens_order, int asens_order){
      for(int i=0; i<NS; ++i){
        PX.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
        PX.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf/NS; 
        for(int j=0; j<3; ++j){
          PX.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*4+j];
        }
        PX.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[0] = f.input()[i*4+3];
      }

      for(int i=0; i<NS; ++i){
        PP.input(INTEGRATOR_T0 + i*INTEGRATOR_NUM_IN)[0] = 0; 
        PP.input(INTEGRATOR_TF + i*INTEGRATOR_NUM_IN)[0] = tf/NS; 
        for(int j=0; j<3; ++j){
          PP.input(INTEGRATOR_X0 + i*INTEGRATOR_NUM_IN)[j] = f.input()[i*4+j];
        }
        PP.input(INTEGRATOR_P  + i*INTEGRATOR_NUM_IN)[0] = f.input()[i*4+3];
      }

      PX.evaluate();
      PP.evaluate();
      
      for(int i=0; i<NS; ++i){
        J3.setInput(PX.output(i),i*2);
        J3.setInput(PP.output(i),i*2+1);
      }
      
      J3.evaluate();
      
      for(int i=0; i<f.getNumOutputs(); ++i){
        J3.getOutput(f.output(i),i);
      }
    }


    MXFunction F, G;
    CFunction J2;
    Jacobian J;

    vector<double> UX_min;
    vector<double> UX_max;
    vector<double> UX_init;

      vector<double> g_min;
      vector<double> g_max;
      Parallelizer PX,PP;
      
      SXFunction J3;
CVodesIntegrator I;

FX IjacP;
};



int main(){

  Jac jac;
  
  IpoptSolver solver(jac.F,jac.G,FX(),jac.J2);
  solver.setOption("tol",1e-5);
  solver.setOption("hessian_approximation", "limited-memory");
  solver.setOption("max_iter",100);
  solver.setOption("linear_solver","ma57");
//  solver.setOption("derivative_test","first-order");

  //solver.setOption("verbose",true);
  solver.init();

  //Set bounds and initial guess
  solver.setInput(jac.UX_min,  NLP_LBX);
  solver.setInput(jac.UX_max,  NLP_UBX);
  solver.setInput(jac.UX_init,  NLP_X_INIT);
  
  solver.setInput(jac.g_min,NLP_LBG);
  solver.setInput(jac.g_max,NLP_UBG);

  //Solve the problem
  solver.solve();

  return 0;
}



