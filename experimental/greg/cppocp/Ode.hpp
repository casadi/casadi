// Ode.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include <string>
#include <map>

#include <symbolic/sx/sx_tools.hpp>
#include <symbolic/fx/sx_function.hpp>

#define IDX_ODE_STEP_STATE  0
#define IDX_ODE_STEP_ACTION 1
#define IDX_ODE_STEP_T0     2
#define IDX_ODE_STEP_DT     3
#define NUM_ODE_STEP_INPUTS 4

class Ode
{

public:
     ~Ode(void);
     Ode(std::string _name);

     std::string name;

     // states/actions/outputs
     std::map<std::string,int> states;
     std::map<std::string,int> actions;
     std::map<std::string,int> outputs;
	
     void addState(std::string _newState);
     void addAction(std::string _newAction);
     void addOutput(std::string _newOutput);

     void (*dxdt)(std::map<std::string,CasADi::SX> &xDot,
		  std::map<std::string,CasADi::SX> &outputs,
		  std::map<std::string,CasADi::SX> state,
		  std::map<std::string,CasADi::SX> action,
		  std::map<std::string,CasADi::SX> param,
		  CasADi::SX t);

     void assertUnlocked(void);
     void init(void);
  void setupIntegrators(void);
  void setupIntegrators(std::map<std::string,double> & params);

     std::map<std::string,CasADi::SX> getOutputFromDxdt( CasADi::SX x,
							 CasADi::SX u,
							 std::map<std::string,CasADi::SX> & p,
							 CasADi::SX t );

     CasADi::SX dxVectorDt( CasADi::SX x,
				  CasADi::SX u,
				  std::map<std::string,CasADi::SX> & p,
				  CasADi::SX t );

     CasADi::DMatrix rk4Step( CasADi::DMatrix & xk,
			      CasADi::DMatrix & uk,
			      double t0,
			      double dt);

     CasADi::DMatrix eulerStep( CasADi::DMatrix & xk,
				CasADi::DMatrix & uk,
				double t0,
				double dt);

     CasADi::SX rk4Step( CasADi::SX x0Vec,
			       CasADi::SX u0Vec,
			       CasADi::SX u1Vec,
			       std::map<std::string,CasADi::SX> & p,
			       CasADi::SX t0,
			       CasADi::SX dt);

     CasADi::SX eulerStep( CasADi::SX x0Vec,
				 CasADi::SX u0Vec,
				 std::map<std::string,CasADi::SX> & p,
				 CasADi::SX t0,
				 CasADi::SX dt);
	
     CasADi::SX simpsonsRuleError( CasADi::SX x0Vec,
					 CasADi::SX x1Vec,
					 CasADi::SX u0Vec,
					 CasADi::SX u1Vec,
					 std::map<std::string,CasADi::SX> & p,
					 CasADi::SX t0,
					 CasADi::SX dt);

     int nx(void);
     int nu(void);
     int nxu(void);
     int no(void);

     void assertUniqueName(std::string newName);

     std::map<std::string,CasADi::SX> getStateMap( CasADi::SX & x);
     std::map<std::string,CasADi::SX> getActionMap( CasADi::SX & u);

private:
     int isState(std::string stateName);
     int isAction(std::string actionName);
     int isOutput(std::string outputName);
     int locked;

     CasADi::SXFunction rk4StepFcn;
     CasADi::SXFunction eulerStepFcn;

protected:
};
