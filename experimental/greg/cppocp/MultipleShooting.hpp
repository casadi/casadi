// MultipleShooting.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <symbolic/sx/sx_tools.hpp>

class MultipleShooting
{
public:
     MultipleShooting(std::string name_,
		      Ode & _ode,
		      CasADi::SX t0_,
		      CasADi::SX tf_,
		      int N_,
		      std::vector<double>&lb_,
		      std::vector<double>&ub_,
		      std::vector<double>&guess_,
		      CasADi::SX & dv_,
		      int idx0_,
		      std::map<std::string,CasADi::SX>&params_);
     ~MultipleShooting();

     int getBigN(void);

     void boundStateAction(std::string xu, double lb_, double ub_, int timeStep);
     void boundStateAction(std::string xu, double lb_, double ub_);
     void setStateActionGuess(std::string xu, double guess_, int timeStep);
     void setStateActionGuess(std::string xu, double guess_);

     void writeOctaveOutput( std::ostream & f, std::vector<double> & xopt);

     CasADi::SX getState(std::string x, int timeStep);
     CasADi::SX getAction(std::string u, int timeStep);
     CasADi::SX getOutput(std::string o, int timeStep);
     CasADi::SX getOutput(std::string o);

     CasADi::SX getStateMat(int timeStep);
     CasADi::SX getActionMat(int timeStep);

     CasADi::DMatrix getState(int timeStep, std::vector<double> & xopt);
     CasADi::DMatrix getAction(int timeStep, std::vector<double> & xopt);


     int N;
	
     CasADi::SX getDynamicsConstraintError(int timeStep);

     Ode & ode;
     int getIdx(std::string xu, int timeStep);
  
     std::string name;

private:
     CasADi::SX & dv;
     std::vector<double>&lb;
     std::vector<double>&ub;
     std::vector<double>&guess;
     std::map<std::string,CasADi::SX> & params;
     int idx0;

     CasADi::SX t0;
     CasADi::SX tf;

};
