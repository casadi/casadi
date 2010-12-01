/*
* kite0ode.h
*/
#ifndef KITE0ODE_H
#define KITE0ODE_H


#include <toolboxes/kinematics/kinetics.hpp>
#include <casadi/stl_vector_tools.hpp>

namespace KINEMATICS{
using namespace CasADi;

#define i_delta 0
#define i_r 1
#define i_phi 2
#define i_theta 3
#define i_R 4
#define i_P 5
#define i_Y 6

#define i_ddelta 7
#define i_dr 8
#define i_dphi 9
#define i_dtheta 10
#define i_dR 11
#define i_dP 12
#define i_dY 13

#define i_E 14
#define i_mu 15
#define i_nu 16

#define N_x 17

#define i_udddelta 0
#define i_uddr 1
#define i_udmuL 2
#define i_udmuR 3
#define i_udnu 4

#define N_u 5

#define N_h 9


/**
* \param t time
* \param q states of the system 
*        delta, r,  phi, theta, R, P, Y;
*   	 ddelta, dr,  dphi, dtheta, dR, dP, dY;
*        E;
*        mu,nu
* \param p parameters of the system
* \param d disturbances on the system
* \param u Input of the system
*        dddelta,ddr,dmu,dnu
* \param h internal expressions which the user may wish to inspect
*/

SXMatrix ode(const SX &t_,const SXMatrix &q_,const SXMatrix &p_,const SXMatrix &d_,const SXMatrix &u_,SXMatrix &h_);

} // namespace KINEMATICS

#endif
