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


#include "nlp_tools.hpp"

#include "casadi_misc.hpp"

namespace casadi {

template <class T>
void detect_simple_bounds_gen(const T& x, const T& p,
      const T& g, const T& lbg, const T& ubg,
      std::vector<casadi_int>& gi,
      T& lbx, T& ubx,
      Function& lam_forward,
      Function& lam_backward) {

    // Read dimensions
    casadi_int ng = g.size1();
    casadi_int nx = x.size1();

    casadi_assert(lbg.size1()==ng, "Dimension mismatch");
    casadi_assert(ubg.size1()==ng, "Dimension mismatch");
    casadi_assert(g.is_column(), "Dimension mismatch");
    casadi_assert(lbg.is_column(), "Dimension mismatch");
    casadi_assert(ubg.is_column(), "Dimension mismatch");
    casadi_assert(x.is_column(), "Dimension mismatch");

    // Get constraint Jacobian sparsity
    Function temp("temp", {x, p}, {g});
    Sparsity sp = temp.sparsity_jac(0, 0);
    Sparsity spT = sp.T();

    // Reset result vector
    std::vector<bool> is_simple(ng, true);

    // Check nonlinearity
    std::vector<bool> is_nonlin = which_depends(g, x, 2, true);

    const casadi_int* row = spT.colind();
    for (casadi_int i=0;i<ng;++i) {
      // Check if each row of jac_g_x only depends on one column
      bool single_dependency = row[i+1]-row[i]==1;
      is_simple[i] = single_dependency && !is_nonlin[i];
    }

    // Full-indices of all simple constraints
    std::vector<casadi_int> sgi = boolvec_to_index(is_simple);
    gi = boolvec_to_index(boolvec_not(is_simple));
    T g_bounds = g(sgi);
    std::vector<casadi_int> sgi_to_sgsi = lookupvector(sgi, ng);

    // Detect  f2(p)x+f1(p)==0
    Function gf = Function("gf", std::vector<T>{x, p},
                                  std::vector<T>{g_bounds, jtimes(g_bounds, x, T::ones(nx, 1))});

    std::vector<T> res;
    gf.call(std::vector<T>{0, p}, res, true, false);
    T f1 = res[0];
    T f2 = res[1];

    T lb = (lbg(sgi)-f1)/abs(f2);
    T ub = (ubg(sgi)-f1)/abs(f2);

    // Start without simple bounds
    lbx = -T::inf(nx);
    ubx = T::inf(nx);

    const casadi_int* xi = spT.row();


    /* We will process in groups
     *  (=collection of bounds that have been re-defined a certain number of times)
     *  instead of element-wise to keep the MX graph size minimal.
    */

    std::vector< std::vector<casadi_int> > sgsi_groups, sgi_groups, sxi_groups;
    casadi_int n_groups = 0;

    // How often has a certain variable been encountered?
    std::vector<casadi_int> group_count(nx);

    // Loop over all constraints
    for (casadi_int i=0;i<ng;++i) {
      // Only treat simple ones
      if (!is_simple[i]) continue;

      casadi_int j = xi[row[i]];
      casadi_int& k = group_count[j];
      // Grow as needed
      if (k==n_groups) {
        sgsi_groups.emplace_back();
        sgi_groups.emplace_back();
        sxi_groups.emplace_back();
        n_groups++;
      }
      sxi_groups[k].push_back(j);
      sgi_groups[k].push_back(i);
      sgsi_groups[k].push_back(sgi_to_sgsi[i]);
      k++;
    }

    casadi_assert(n_groups>=1, "mismatch");
    // Take min/max group-wise to determine simple bounds
    for (casadi_int i=0;i<n_groups;++i) {
      lbx(sxi_groups[i]) = fmax(lbx(sxi_groups[i]), lb(sgsi_groups[i]));
      ubx(sxi_groups[i]) = fmin(ubx(sxi_groups[i]), ub(sgsi_groups[i]));
    }


    // Choose the first group as reference
    const std::vector<casadi_int>& sxi_ref = sxi_groups[0];
    // Indices into the first group that reproduce the contents of other groups
    std::vector< std::vector<casadi_int> > xsub(n_groups);
    // Identity map for first entry of xsub
    xsub[0] = range(sxi_ref.size());

    // Computations for other xsub entries
    std::vector<casadi_int> lookup = lookupvector(sxi_ref, nx);
    for (casadi_int i=1;i<n_groups;++i) {
      const std::vector<casadi_int>& sxi = sxi_groups[i];
      for (casadi_int j=0;j<sxi.size();++j) {
        xsub[i].push_back(lookup[sxi[j]]);
      }
    }

    // Determine multiplier maps

    // Symbols
    T lam_g_forward = T::sym("lam_g", ng);
    T lam_sg_backward = T::sym("lam_g", gi.size());
    T lam_x_backward = T::sym("lam_x", nx);


    T lam_sg_forward = lam_g_forward(sgi); // NOLINT(cppcoreguidelines-slicing)
    T lam_x_forward = T::zeros(nx, 1);
    T lam_g_backward = T::zeros(ng, 1);
    lam_g_backward(gi) = lam_sg_backward;

    // Comparison expression per group
    std::vector<T> lcomp, ucomp;
    for (casadi_int i=0;i<n_groups;++i) {
      lcomp.push_back(lbx(sxi_groups[i])==lb(sgsi_groups[i]));
      ucomp.push_back(ubx(sxi_groups[i])==ub(sgsi_groups[i]));
    }

    // How many lb/ub are active?
    T count_lb = T::zeros(nx);
    T count_ub = T::zeros(nx);
    for (casadi_int i=0;i<n_groups;++i) {
      count_lb(sxi_groups[i]) += lcomp[i];
      count_ub(sxi_groups[i]) += ucomp[i];
    }

    // Compute lam_x from lam_g
    for (casadi_int i=0;i<n_groups;++i) {
      T l = lam_sg_forward(sgsi_groups[i]); // NOLINT(cppcoreguidelines-slicing)
      T lt = (l<0);
      T gt = !lt;
      lam_x_forward(sxi_groups[i]) += l*(lt*lcomp[i]+gt*ucomp[i]);
    }

    // Compute lam_g from lam_x
    T l = lam_x_backward(sxi_ref); // NOLINT(cppcoreguidelines-slicing)
    T lt = (l<0);
    T gt = !lt;
    T lam_xl = l*lt/count_lb(sxi_ref);
    T lam_xu = l*gt/count_ub(sxi_ref);

    for (casadi_int i=0;i<n_groups;++i) {
      lam_g_backward(sgi_groups[i]) = lam_xl(xsub[i])*lcomp[i]+lam_xu(xsub[i])*ucomp[i];
    }

    // Construct Functions for mappings
    lam_forward = Function("lam_forward",
      {lam_g_forward, p}, {lam_g_forward(gi), lam_x_forward},// NOLINT(cppcoreguidelines-slicing)
      {"lam_g", "p"}, {"lam_sg", "lam_x"});
    casadi_assert_dev(!lam_forward.has_free());
    lam_backward = Function("lam_backward",
      {lam_sg_backward, lam_x_backward, p}, {lam_g_backward},
      {"lam_sg", "lam_x", "p"}, {"lam_g"});
    casadi_assert_dev(!lam_backward.has_free());
  }

void detect_simple_bounds(const SX& x, const SX& p,
      const SX& g, const SX& lbg, const SX& ubg,
      std::vector<casadi_int>& gi,
      SX& lbx, SX& ubx,
      Function& lam_forward,
      Function& lam_backward) {
  detect_simple_bounds_gen(x, p, g, lbg, ubg,
    gi, lbx, ubx, lam_forward, lam_backward);
}

void detect_simple_bounds(const MX& x, const MX& p,
      const MX& g, const MX& lbg, const MX& ubg,
      std::vector<casadi_int>& gi,
      MX& lbx, MX& ubx,
      Function& lam_forward,
      Function& lam_backward) {
  detect_simple_bounds_gen(x, p, g, lbg, ubg,
    gi, lbx, ubx, lam_forward, lam_backward);
}

} // namespace casadi
