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


#include "bspline.hpp"
#include "interpolant_impl.hpp"
#include "casadi_low.hpp"

using namespace std;
namespace casadi {

  BSpline::BSpline(const MX& x, const MX& knots,
          const std::vector<casadi_int>& offset,
          const MX& coeffs,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const std::vector<casadi_int>& lookup_mode) :
          knots_(knots.eval()), coeffs_(coeffs.eval()), offset_(offset), degree_(degree),
          m_(m), lookup_mode_(lookup_mode) {

    if (knots_->type_ & MXNode::MX_PAR) knots_ = knots_->get_gate();
    if (coeffs_->type_ & MXNode::MX_PAR) coeffs_ = coeffs_->get_gate();

    prepare(m_, offset_, degree_, coeffs_size_, coeffs_dims_, strides_);
    casadi_assert_dev(x.numel()==degree.size());
    set_dep(x);
    set_sparsity(Sparsity::dense(m, 1));

    try {
      knots_ptr_ = &static_cast<const std::vector<double> &>(knots_);
    } catch (...) {
      knots_vec_ = static_cast<DM>(knots_).nonzeros();
      knots_ptr_ = &knots_vec_;
    }
    try {
      coeffs_ptr_ = &static_cast<const std::vector<double> &>(coeffs_);
    } catch (...) {
      coeffs_vec_ = static_cast<DM>(coeffs_vec_).nonzeros();
      coeffs_ptr_ = &coeffs_vec_;
    }

    casadi_assert_dev(knots_ptr_);
    casadi_assert_dev(coeffs_ptr_);
  }

  size_t BSpline::sz_iw() const {
    return n_iw(degree_);
  }

  size_t BSpline::sz_w() const {
    return n_w(degree_);
  }

  size_t BSpline::n_iw(const std::vector<casadi_int>& degree) {
    casadi_int n_dims = degree.size();
    casadi_int sz = 0;
    sz += n_dims+1; // boor_offset
    sz += n_dims; // starts
    sz += n_dims; // index
    sz += n_dims+1; // coeff_offset
    return sz;
  }

  size_t BSpline::n_w(const std::vector<casadi_int>& degree) {
    casadi_int n_dims = degree.size();
    casadi_int sz = 0;
    for (casadi_int k=0;k<n_dims-1;++k) {
      sz += degree[k]+1; // boor
    }
    sz += 2*degree[n_dims-1]+1;
    sz += n_dims+1;
    return sz;
  }

  casadi_int BSpline::get_coeff_size(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree) {
    casadi_int ret = m;
    for (casadi_int i=0;i<degree.size();++i) ret*= offset[i+1]-offset[i]-degree[i]-1;
    return ret;
  }

  void BSpline::prepare(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree, casadi_int &coeffs_size,
      std::vector<casadi_int>& coeffs_dims, std::vector<casadi_int>& strides) {

    casadi_int n_dims = degree.size();
    coeffs_size = get_coeff_size(m, offset, degree);
    coeffs_dims.resize(n_dims+1);
    coeffs_dims[0] = m;
    for (casadi_int i=0;i<n_dims;++i) coeffs_dims[i+1] = offset[i+1]-offset[i]-degree[i]-1;

    // Prepare strides
    strides.resize(n_dims);
    strides[0] = m;
    for (casadi_int i=0;i<n_dims-1;++i) {
      strides[i+1] = strides[i]*coeffs_dims[i+1];
    }
  }

  void BSpline::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("BSpline::knots", knots_);
    s.pack("BSpline::coeffs", coeffs_);
    s.pack("BSpline::offset", offset_);
    s.pack("BSpline::degree", degree_);
    s.pack("BSpline::m", m_);
    s.pack("BSpline::lookup_mode", lookup_mode_);
    s.pack("BSpline::strides", strides_);
    s.pack("BSpline::coeffs_dims", coeffs_dims_);
    s.pack("BSpline::coeffs_size", coeffs_size_);
    s.pack("BSpline::jac_cache_", jac_cache_);
  }

  BSpline::BSpline(DeserializingStream& s) : MXNode(s) {
    s.unpack("BSpline::knots", knots_);
    s.unpack("BSpline::coeffs", coeffs_);
    s.unpack("BSpline::offset", offset_);
    s.unpack("BSpline::degree", degree_);
    s.unpack("BSpline::m", m_);
    s.unpack("BSpline::lookup_mode", lookup_mode_);
    s.unpack("BSpline::strides", strides_);
    s.unpack("BSpline::coeffs_dims", coeffs_dims_);
    s.unpack("BSpline::coeffs_size", coeffs_size_);
    s.unpack("BSpline::jac_cache_", jac_cache_);

    try {
      knots_ptr_ = &static_cast<const std::vector<double> &>(knots_);
    } catch (...) {
      knots_vec_ = static_cast<DM>(knots_).nonzeros();
      knots_ptr_ = &knots_vec_;
    }
    try {
      coeffs_ptr_ = &static_cast<const std::vector<double> &>(coeffs_);
    } catch (...) {
      coeffs_vec_ = static_cast<DM>(coeffs_vec_).nonzeros();
      coeffs_ptr_ = &coeffs_vec_;
    }
  }

  void get_boor(const MX& x, const MX& knots, casadi_int degree, casadi_int lookup_mode,
      MX& start, MX& boor) {
    MX knots_clipped = knots(range(degree, knots.size1()-degree));

    Dict low_opts;
    low_opts["lookup_mode"] = Low::lookup_mode_from_enum(lookup_mode);
    MX L = low(knots_clipped, x, low_opts);
    start = fmin(L, knots.size1()-2*degree-2);

    DM boor_init = DM::zeros(x.size2(), 2*degree+1);
    boor_init(Slice(), degree) = 1;
    std::vector<MX> boor_full = horzsplit(MX(boor_init));
    casadi_int n_knots = 2*degree+2;

    MX kn;
    MX(knots).get_nz(kn, false, start, MX(range(n_knots)));

    std::vector<MX> knv = horzsplit(kn);

    MX xt = x.T();

    for (casadi_int d=1;d<degree+1;++d) {
      for (casadi_int i=0;i<n_knots-d-1;++i) {
        MX bottom = knv[i+d]-knv[i];
        MX b = if_else_zero(bottom, (xt-knv[i])*boor_full[i]/(bottom+1e-100));
        bottom = knv[i+d+1]-knv[i + 1];
        b += if_else_zero(bottom, (knv[i+d+1]-xt)*boor_full[i+1]/(bottom+1e-100));
        boor_full[i] = b;
      }
    }

    boor = horzcat(std::vector<MX>(boor_full.begin(), boor_full.begin()+degree+1));
  }

  MX do_inline(const MX& x,
                const std::vector< MX >& knots,
                const MX& coeffs,
                casadi_int m,
                const std::vector<casadi_int>& degree,
                const std::vector<casadi_int>& lookup_mode) {

    casadi_int batch_x = x.size2();

    // Number of grid points
    casadi_int N = knots.size();
    std::vector<MX> xs = vertsplit(x);

    // Compute De Boor vector in each direction
    std::vector<MX> starts(N);
    std::vector< std::vector<MX> > boors(N);
    for (casadi_int i=0;i<N;++i) {
      MX boor;
      get_boor(xs[i], knots[i], degree[i], lookup_mode[i], starts[i], boor);
      boors[i] = horzsplit(boor.T());
    }

    // Compute strides
    std::vector<casadi_int> strides = {m};
    for (casadi_int i=0;i<N-1;++i) {
      strides.push_back(strides.back()*(knots[i].numel()-degree[i]-1));
    }

    // Start index of subtensor: row vector
    MX start = mtimes(DM(strides).T(), vertcat(starts));

    // Elements of subtensor
    DM core = DM(range(m));
    for (casadi_int i=0;i<N;++i) {
      casadi_int n = degree[i]+1;
      core = vec(repmat(core, 1, n)+repmat(strides[i]*DM(range(n)).T(), core.size1(), 1));
    }

    std::vector<MX> res;

    for (casadi_int k=0;k<batch_x;++k) {

      // Flattened subtensor of coefficients
      MX c = reshape(coeffs(start(k)+core), m, -1);

      // Compute outer product of De Boor vectors
      MX boor = 1;
      for (casadi_int i=0;i<N;++i) {
        boor = vec(mtimes(boor, boors[i][k].T()));
      }

      res.push_back(mtimes(c, boor));
    }

    return horzcat(res);
  }

  MX BSpline::create(const MX& x, const std::vector< MX >& knots,
          const MX& coeffs,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const Dict& opts) {


    bool do_inline_flag = false;
    std::vector<std::string> lookup_mode;

    for (auto&& op : opts) {
      if (op.first=="inline") {
        do_inline_flag = op.second;
      } else if (op.first=="lookup_mode") {
        lookup_mode = op.second;
      }
    }

    std::vector<casadi_int> offset;
    MX stacked;
    Interpolant::stack_grid(knots, offset, stacked);

    std::vector<casadi_int> mode =
      Interpolant::interpret_lookup_mode(lookup_mode, stacked, offset, degree, degree);

    if (do_inline_flag) {
      return do_inline(x, knots, coeffs, m, degree, mode);
    } else {
      return x->get_bspline(coeffs, stacked, offset, degree, m, mode);
    }
  }

  std::string BSpline::disp(const std::vector<std::string>& arg) const {
    return "BSpline(" + arg.at(0) + ")";
  }

  void BSpline::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_bspline(knots_, coeffs_, offset_, degree_, m_, lookup_mode_);
  }

  MX BSpline::jac_cached() const {
    if (jac_cache_.is_empty()) {
      jac_cache_ = jac(dep(0), coeffs_, knots_);
    }
    return jac_cache_;
  }

  void BSpline::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    MX J = jac_cached();

    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = mtimes(J, fseed[d][0]);
    }
  }

  void BSpline::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    MX JT = jac_cached().T();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += mtimes(JT, aseed[d][0]);
    }
  }

  int BSpline::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (!res[0]) return 0;

    casadi_clear(res[0], m_);
    casadi_nd_boor_eval(res[0], degree_.size(), get_ptr(*knots_ptr_), get_ptr(offset_),
      get_ptr(degree_), get_ptr(strides_), get_ptr(*coeffs_ptr_), m_, arg[0], get_ptr(lookup_mode_),
      iw, w);
    return 0;
  }

  void BSpline::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    casadi_int n_dims = offset_.size()-1;

    g.add_auxiliary(CodeGenerator::AUX_ND_BOOR_EVAL);
    g.add_auxiliary(CodeGenerator::AUX_FILL);
    g << g.clear(g.work(res[0], m_), m_) << "\n";

    // Input and output buffers
    g << "CASADI_PREFIX(nd_boor_eval)(" << g.work(res[0], m_) << "," << n_dims << ","
      << generate_grid(g) << "," << g.constant(offset_) << "," <<  g.constant(degree_)
      << "," << g.constant(strides_) << "," << generate_coeff(g) << "," << m_  << ","
      << g.work(arg[0], n_dims) << "," <<  g.constant(lookup_mode_) << ", iw, w);\n";
  }

  std::string BSpline::generate_coeff(CodeGenerator& g) const {
    return (coeffs_->type_ & MXNode::MX_PAR) ? g.rw_double(coeffs_.get()) : g.constant(*coeffs_ptr_);
  }

  std::string BSpline::generate_grid(CodeGenerator& g) const {
    return (knots_->type_ & MXNode::MX_PAR) ? g.rw_double(knots_.get()) : g.constant(*knots_ptr_);
  }

  MX BSpline::dual(const MX& x,
          const std::vector< MX >& knots,
          const std::vector<casadi_int>& degree,
          const Dict& opts) {

    std::vector<casadi_int> offset;
    MX stacked;
    Interpolant::stack_grid(knots, offset, stacked);

    std::vector<std::string> lookup_mode;
    auto it = opts.find("lookup_mode");
    if (it!=opts.end()) lookup_mode = it->second;
    std::vector<casadi_int> lookup_mode_int =
      Interpolant::interpret_lookup_mode(lookup_mode, stacked, offset, degree, degree);

    std::vector<double> stacked_numeric = static_cast<DM>(stacked).nonzeros();
    std::vector<double> x_numeric = static_cast<DM>(x).nonzeros();

    casadi_int n_dims = degree.size();
    casadi_int N = x.numel()/n_dims;
    casadi_assert_dev(N*n_dims==x.numel());

    casadi_int coeffs_size;
    std::vector<casadi_int> coeffs_dims, strides;
    prepare(1, offset, degree, coeffs_size, coeffs_dims, strides);

    // Size of coefficients
    std::vector<double> contribution(coeffs_size);
    std::vector<casadi_int> nz(coeffs_size);

    std::vector<double> data;
    std::vector<casadi_int> row, col;

    std::vector<double> w(n_w(degree));
    std::vector<casadi_int> iw(n_iw(degree));

    for (casadi_int i=0;i<N;++i) {
      std::fill(contribution.begin(), contribution.end(), 0.0);
      casadi_int nnz = casadi_nd_boor_dual_eval(get_ptr(contribution), get_ptr(nz),
        degree.size(), get_ptr(stacked_numeric), get_ptr(offset),
        get_ptr(degree), get_ptr(strides), get_ptr(x_numeric)+i*n_dims, get_ptr(lookup_mode_int),
        get_ptr(iw), get_ptr(w));
      data.insert(data.end(), contribution.begin(), contribution.begin()+nnz);
      col.insert(col.end(), nz.begin(), nz.begin()+nnz);
      row.insert(row.end(), nnz, i);
    }

    return DM(Sparsity::triplet(coeffs_size, N, col, row), data).T();
  }

    MX BSpline::derivative_coeff(casadi_int i, const MX& coeffs, const MX& knots) const {
      casadi_int n_dims = degree_.size();

      casadi_int n_knots = offset_[i+1]-offset_[i];
      casadi_int n = n_knots-degree_[i]-1;
      MX delta_knots = knots(range(offset_[i]+1+degree_[i], offset_[i]+n_knots-1))
            - knots(range(offset_[i]+1, offset_[i]+n_knots-degree_[i]-1));
      Sparsity sp_diag = vertsplit(Sparsity::diag(n), {0, n-1, n})[0];
      Sparsity sp_band = vertsplit(Sparsity::band(n, -1), {0, n-1, n})[0];
      MX delta_knots_inv = 1/delta_knots;
      MX T = MX(sp_diag, -delta_knots_inv) + MX(sp_band, delta_knots_inv);
      T*= degree_[i];

      std::vector<casadi_int> coeffs_dims_new = coeffs_dims_;
      coeffs_dims_new[i+1] = T.size1();

      // Apply transformation T on axis i

      // Bring axis i to the back
      std::vector<casadi_int> order = range(n_dims+1);
      std::swap(order.back(), order[i+1]);
      std::vector<casadi_int> mapping = tensor_permute_mapping(coeffs_dims_, order);
      MX coeff_matrix = coeffs.nz(mapping); // NOLINT(cppcoreguidelines-slicing)

      // Cast as matrix
      coeff_matrix = reshape(coeff_matrix, -1, T.size2());

      // Apply the transformation matrix from the right
      coeff_matrix = mtimes(coeff_matrix, T.T());

      // Bring axis i back to the original place
      mapping = tensor_permute_mapping(permute(coeffs_dims_new, order), order);
      coeff_matrix = coeff_matrix.nz(mapping); // NOLINT(cppcoreguidelines-slicing)

      // Return the flat vector
      return coeff_matrix;
    }

    MX BSpline::jac(const MX& x, const MX& coeffs, const MX& kn) const {
    casadi_int n_dims = degree_.size();
    std::vector<MX> parts;

    Dict opts;
    std::vector<std::string> lookup_mode;
    for (auto e : lookup_mode_) lookup_mode.push_back(Low::lookup_mode_from_enum(e));
    opts["lookup_mode"] = lookup_mode;
    
    // Loop over dimensions
    for (casadi_int k=0;k<n_dims;++k) {
      std::vector< MX > knots;
      std::vector< casadi_int> degree;
      for (casadi_int i=0;i<degree_.size();++i) {
        if (i==k) {
          knots.push_back(kn(range(offset_[i]+1, offset_[i+1]-1)));
          degree.push_back(degree_[i]-1);
        } else {
          knots.push_back(kn(range(offset_[i], offset_[i+1])));
          degree.push_back(degree_[i]);
        }
      }
      MX d = MX::bspline(x, derivative_coeff(k, coeffs, kn), Interpolant::parse_grid(knots), degree, m_, opts);
      parts.push_back(d);
    }

    return horzcat(parts);
  }

} // namespace casadi
