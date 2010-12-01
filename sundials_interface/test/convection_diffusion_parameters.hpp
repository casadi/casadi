#ifndef CONVECTION_DIFFUSION_PARAMETERS_HPP
#define CONVECTION_DIFFUSION_PARAMETERS_HPP

#include <limits>

using namespace std;

// Spatial discretization in storage
static int NZ = 100;

// receiver temperature [C]
static double T_receiver = 680;

// return air temperature [C]
static double T_return = 120;

// diffusion coefficient
static double diff = 1e-3;

// height of the storage
static double h = 1;

// outputs discretization
static int nt_out = 100; 
static int nz_out = 100;

// be verbose
static bool verbose = false;

// infinity
static double infty = numeric_limits<double>::infinity();

// Number of times
static int nk_per_hour = 2;

// Time for a control
static double DT = 1.0/nk_per_hour;

// Discretize space
static const int x_dim = 2; // number of state dimensions
static const int u_dim = 1; // number of control dimensions
  
// Discretization for each dimension
static int n_x_disc[] = {100, 100};

// Lower and upper bound on each variable
static double x_min[] = {-1, 0.01};
static double x_max[] = { 1, 1   };

// Discretize the control
static const int n_u_pos = 20; // in the positive direction direction
 //  const int n_u_pos = 4; // in the positive direction direction
static const int n_u_disc = 2*n_u_pos+1; // n_u_pos downwind, zero and n_u_pos upwind

static int num_days = 19;

static int NK = 24*nk_per_hour*num_days;

// types of objective functions
enum OBJECTIVE{LEAST_SQUARE, ENERGY};
static OBJECTIVE obj = ENERGY;

#endif // CONVECTION_DIFFUSION_PARAMETERS_HPP
