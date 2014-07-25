// MultipleShooting.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include "Ode.hpp"

#include <string>
#include <map>
#include <fstream>

#include <core/sx/sx_tools.hpp>

class MultipleShooting
{
public:
     MultipleShooting(std::string name_,
		      Ode & _ode,
		      casadi::SX t0_,
		      casadi::SX tf_,
		      int N_,
		      std::vector<double>&lb_,
		      std::vector<double>&ub_,
		      std::vector<double>&guess_,
		      casadi::SX & dv_,
		      int idx0_,
		      std::map<std::string,casadi::SX>&params_);
     ~MultipleShooting();

     int getBigN(void);

     void boundStateAction(std::string xu, double lb_, double ub_, int timeStep);
     void boundStateAction(std::string xu, double lb_, double ub_);
     void setStateActionGuess(std::string xu, double guess_, int timeStep);
     void setStateActionGuess(std::string xu, double guess_);

     void writeOctaveOutput( std::ostream & f, std::vector<double> & xopt);

     casadi::SX getState(std::string x, int timeStep);
     casadi::SX getAction(std::string u, int timeStep);
     casadi::SX getOutput(std::string o, int timeStep);
     casadi::SX getOutput(std::string o);

     casadi::SX getStateMat(int timeStep);
     casadi::SX getActionMat(int timeStep);

     casadi::DMatrix getState(int timeStep, std::vector<double> & xopt);
     casadi::DMatrix getAction(int timeStep, std::vector<double> & xopt);


     int N;
	
     casadi::SX getDynamicsConstraintError(int timeStep);

     Ode & ode;
     int getIdx(std::string xu, int timeStep);
  
     std::string name;

private:
     casadi::SX & dv;
     std::vector<double>&lb;
     std::vector<double>&ub;
     std::vector<double>&guess;
     std::map<std::string,casadi::SX> & params;
     int idx0;

     casadi::SX t0;
     casadi::SX tf;

};
