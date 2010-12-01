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

#include <toolboxes/kinematics/test/kiteode.hpp>

#include <casadi/expression_tools.hpp>
#include <toolboxes/kinematics/test/simulator.hpp>

using namespace KINEMATICS;
using namespace std;

ModelSimulator::ModelSimulator() {

    int nStates=N_x;
    SXMatrix q("q",N_x,1);
    SXMatrix u("u",4,1);
    
    SX time=0.0;

    SXMatrix p = 0.0;
    SXMatrix d = 0.0;
    SXMatrix h = zeros(N_h,1);
    SXMatrix F=ode(0,q,p,d,u,h);
    
    std::vector<SXMatrix> arg;
    arg.push_back(q);
    arg.push_back(u);
  
    f=SXFunction(arg,F);

    SXMatrix J=jacobian(F,q);

    SXMatrix arg2;arg2 << q << u;
    j=SXFunction(arg,J);

    g=SXFunction(q,h);

    state= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) state[i]=0;

    input= ublas::vector<double>(N_u);for (int i=0;i<N_u;i++) input[i]=0;

    k1= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) k1[i]=0;
    k2= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) k2[i]=0;
    k3= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) k3[i]=0;
    k4= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) k4[i]=0;
    state_= ublas::vector<double>(N_x);for (int i=0;i<N_x;i++) state_[i]=0;



    A=ublas::matrix<double>(N_x,N_x);
    b=ublas::vector<double>(N_x);
    t=0;
    deltat=5e-4;
    Deltat=1e-2;

}

void ModelSimulator::linearize(double t,const ublas::vector<double> &q,const ublas::vector<double> &u) {
	b=evaluate(t,q,u);
	ublas::vector<double> J((N_x+N_u)*N_x);
	for(int i=0;i<N_x;i++) j->input[0].data()[i]=q(i);
	for(int i=0;i<N_u;i++) j->input[0].data()[N_x+i]=u(i);
	j.evaluate();
	j.getOutputDense(&J[0]); // passing a ublas::matrix here &J(0,0) seems to not fully work. yes it does

	for(int i=0;i<N_x;i++) {
		for(int k=0;k<N_x;k++) {
			A(i,k)=J(i*N_x+k);
		}
	}
	q0=q;
	u0=u;

	ublas::vector<double> dq(N_x);
	for(int i=0;i<N_x;i++) dq[i]=1e-4;

}

void ModelSimulator::setdeltat(double t) {
	deltat=t;
}
ublas::vector<double> ModelSimulator::output() {
	ublas::vector<double> res(N_h);
	for(int i=0;i<17;i++) g->input[0].data()[i]=state(i);
	g.evaluate();
	g.getOutputDense(&res[0]);
	return res;
}

void ModelSimulator::setInput(ublas::vector<double> input_) {
	input=input_;
}

ublas::vector<double> ModelSimulator::evaluate(double t,const ublas::vector<double> &q,const ublas::vector<double> &u) {
	ublas::vector<double> res(N_x);
	for(int i=0;i<N_x;i++) f->input[0].data()[i]=q(i);
	for(int i=0;i<N_u;i++) f->input[1].data()[i]=u(i);
	f.evaluate();
	f.getOutputDense(&res[0]);
	return res;
}

ublas::vector<double> ModelSimulator::evaluatelin(double t,const ublas::vector<double> &dq) {
	return prod(A,dq)+b;
}


void ModelSimulator::integrateOnce(double dt) {
	state_=state-q0;
	k1=dt*evaluatelin(t,state_);
	k2=dt*evaluatelin(t+1/2.0*dt,state_+1/2.0*k1);
	k3=dt*evaluatelin(t+1/2.0*dt,state_+1/2.0*k2);
	k4=dt*evaluatelin(t+1/2.0*dt,state_+1/2.0*k3);
	state+=k1/6.0+k2/3.0+k3/3.0+k4/6.0;
	t+=dt;
}

void ModelSimulator::integratelin(double dt_) {
	linearize(t,state,input);
	int N=ceil(dt_/deltat);
        double dt=dt_/N;
	for (int i=0;i<N;i++) integrateOnce(dt);
}

void ModelSimulator::integrate(double dt_) {
	int N=ceil(dt_/Deltat);
        double dt=dt_/N;
	for (int i=0;i<N;i++) integratelin(dt);
}

void ModelSimulator::init(const ublas::vector<double> &q) {
	state=q;
}



int main( ){
	cout << "This is a simulator" << endl;
	double x0test[17]={0, 1, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0.0,0,0,0};

	ublas::vector<double> ini(N_x);
	for (int i=0;i<N_x;i++) ini[i]=x0test[i];

        cout << "ini:" << ini << endl;

	ModelSimulator M;
	M.init(ini);

	

	ublas::vector<double> out=M.output();
        cout << "output:" << out << endl;
	out=M.output();
        cout  << "output:" << out << endl;
	M.integrate(1);
	out=M.output();
        cout  << "output:" << out << endl;
	M.integrate(1);
	out=M.output();
        cout  << "output:" << out << endl;

}


void ModelSimulator::test() {
	cout << "Test is succesful" << endl;
}

ModelSimulator::~ModelSimulator() {
	cout << "It was nice to meet you" << endl;
	f->clear();
	g->clear();
	j->clear();
}
