/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#include <casadi/casadi.hpp>
/** Solve a QP using the low-level (conic) interface
  * The example below is QRECIPE from CUTE, borrowed from the qpOASES examples
  * Joel Andersson, 2016
  */

using namespace casadi;
using namespace std;

// Matrix H in sparse triplet format
const int H_nrow = 2;
const int H_ncol = 2;
const vector<casadi_int> H_colind = {
  0,  1,  2
};
const vector<casadi_int> H_row = {
  0, 1
};
const vector<double> H_nz = {
  2, 2
};

// Matrix A in sparse triplet format
const int A_nrow = 1;
const int A_ncol = 2;
const vector<casadi_int> A_colind = {
  0, 1, 2
};
const vector<casadi_int> A_row = {
  0, 0
};
const vector<double> A_nz = {
  1, 1
};

const vector<double> g = {
  0, 0
};
const vector<double> lbx = {
  1, 1
};
const vector<double> ubx = {
  inf, inf
};

const vector<double> lba = {
  2.0001
};
const vector<double> uba = {
  inf
};
const vector<double> x0 = {2, 2};

const vector<double> lam_x0 = {-1, -1};
const vector<double> lam_a0 = {-1};

int main(){
  // Create QP matrices
  DM H(Sparsity(H_nrow, H_ncol, H_colind, H_row), H_nz);
  DM A(Sparsity(A_nrow, A_ncol, A_colind, A_row), A_nz);

  // Create conic solver
  SpDict qp = {{"h", H.sparsity()}, {"a", A.sparsity()}};
  Dict opts;
  opts["verbose"] = true;
  opts["max_iter"] = 10;
  opts["print_iter"] = true;
  Function F = conic("F", "qrqp", qp, opts);

  // Get the optimal solution
  DMDict arg = {{"h", H}, {"a", A}, {"g", g},
                {"lbx", lbx}, {"ubx", ubx},
                {"lba", lba}, {"uba", uba},
                {"x0", x0}, {"lam_x0", lam_x0}, {"lam_a0", lam_a0}};
  DMDict res = F(arg);
  cout << "res = " << res << endl;

  return 0;
}
