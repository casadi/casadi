// Ocp.hpp
// Greg Horn
// Casadi 2011

#pragma once

#include <iostream>
#include <string>
#include <map>

#include "Ode.hpp"
#include "MultipleShooting.hpp"

#include <symbolic/sx/sx_tools.hpp>
#include <symbolic/std_vector_tools.hpp>

class Ocp
{
public:
     Ocp(void);
     ~Ocp(void);

     void addNonlconIneq(casadi::SX gNew);
     void addNonlconIneq(casadi::SX gNew, std::string name);      
     void addNonlconEq(casadi::SX gNew);
     void addNonlconEq(casadi::SX gNew, std::string name);

     casadi::SX designVariables;
     std::vector<double>lb;
     std::vector<double>ub;
     std::vector<double>guess;

     casadi::SX objFun;
     casadi::SX g;
     std::vector<int> gSizes;
     std::vector<std::string> gLabels;
     std::vector<double> gMin;
     std::vector<double> gMax;

     MultipleShooting & addMultipleShooting(std::string name, Ode & ode, casadi::SX t0, casadi::SX tf, int N);

     casadi::SX & addParam(std::string _newParam);
     casadi::SX & getParam(std::string p);
     double getParamSolution(std::string p);
     void boundParam(std::string p, double lb_, double ub_);
     void setParamGuess(std::string p, double guess_);
        
     casadi::DMatrix getStateSolution(int timestep);
     casadi::DMatrix getActionSolution(int timestep);

     void setState(  casadi::DMatrix x, int timestep );
     void setAction( casadi::DMatrix u, int timestep );
     void setStates(  std::vector<casadi::DMatrix> & x);
     void setActions( std::vector<casadi::DMatrix> & u);

     void writeSolution( std::ostream & f );
     void writeSolution( const char * filename );
     void writeOctaveOutput( std::string name );
     void writeOctaveOutput( std::string name, std::string add );
     void writeOctaveOutput( std::string name, std::string add, std::string filename );
     void writeOctaveOutput( std::string name, std::string add, std::ostream & f );
     void loadGuess( std::istream & f );
     void loadGuess( const char * filename );

     // multiple shooting instances
     std::map<std::string,MultipleShooting*> ms;

     std::vector<double> xopt;

private:

     void assertUniqueName(std::string s);
     int isMultipleShooting(std::string msName);



     // params
     std::map<std::string,int> paramsIdx;
     std::map<std::string,casadi::SX> params;
     int np(void);
     int isParam(std::string paramName);


     casadi::SX getParamMat(void);

     int getParamIdx(std::string p);



protected:
};
