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


#include "bspline_impl.hpp"

using namespace std;
namespace casadi {

  MXNode* BSplineCommon::deserialize(DeserializingStream& s) {
    char t;
    s.unpack("BSpline::type", t);
    switch (t) {
      case 'n':
        return new BSpline(s);
      case 'p':
        return new BSplineParametric(s);
      default:
        casadi_error("Unknown BSpline type");
    }
  }

  /// Constructor
  BSplineCommon::BSplineCommon(const std::vector<double>& knots,
            const std::vector<casadi_int>& offset,
            const std::vector<casadi_int>& degree,
            casadi_int m,
            const std::vector<casadi_int>& lookup_mode) :
          knots_(knots), offset_(offset), degree_(degree),
          m_(m), lookup_mode_(lookup_mode) {
    prepare(m_, offset_, degree_, coeffs_size_, coeffs_dims_, strides_);
  }

  size_t BSplineCommon::sz_iw() const {
    return n_iw(degree_);
  }

  size_t BSplineCommon::sz_w() const {
    return n_w(degree_);
  }

  size_t BSplineCommon::n_iw(const std::vector<casadi_int>& degree) {
    casadi_int n_dims = degree.size();
    casadi_int sz = 0;
    sz += n_dims+1; // boor_offset
    sz += n_dims; // starts
    sz += n_dims; // index
    sz += n_dims+1; // coeff_offset
    return sz;
  }

  size_t BSplineCommon::n_w(const std::vector<casadi_int>& degree) {
    casadi_int n_dims = degree.size();
    casadi_int sz = 0;
    for (casadi_int k=0;k<n_dims-1;++k) {
      sz += degree[k]+1; // boor
    }
    sz += 2*degree[n_dims-1]+1;
    sz += n_dims+1;
    return sz;
  }

  casadi_int BSplineCommon::get_coeff_size(casadi_int m, const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& degree) {
    casadi_int ret = m;
    for (casadi_int i=0;i<degree.size();++i) ret*= offset[i+1]-offset[i]-degree[i]-1;
    return ret;
  }

  void BSplineCommon::prepare(casadi_int m, const std::vector<casadi_int>& offset,
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

  void BSpline::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("BSpline::type", 'n');
  }

  void BSplineParametric::serialize_type(SerializingStream& s) const {
    MXNode::serialize_type(s);
    s.pack("BSpline::type", 'p');
  }

  void BSplineCommon::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    s.pack("BSplineCommon::knots", knots_);
    s.pack("BSplineCommon::offset", offset_);
    s.pack("BSplineCommon::degree", degree_);
    s.pack("BSplineCommon::m", m_);
    s.pack("BSplineCommon::lookup_mode", lookup_mode_);
    s.pack("BSplineCommon::strides", strides_);
    s.pack("BSplineCommon::coeffs_dims", coeffs_dims_);
    s.pack("BSplineCommon::coeffs_size", coeffs_size_);
    s.pack("BSplineCommon::jac_cache_", jac_cache_);
  }

  BSplineCommon::BSplineCommon(DeserializingStream& s) : MXNode(s) {
    s.unpack("BSplineCommon::knots", knots_);
    s.unpack("BSplineCommon::offset", offset_);
    s.unpack("BSplineCommon::degree", degree_);
    s.unpack("BSplineCommon::m", m_);
    s.unpack("BSplineCommon::lookup_mode", lookup_mode_);
    s.unpack("BSplineCommon::strides", strides_);
    s.unpack("BSplineCommon::coeffs_dims", coeffs_dims_);
    s.unpack("BSplineCommon::coeffs_size", coeffs_size_);
    s.unpack("BSplineCommon::jac_cache_", jac_cache_);
  }

  void BSpline::serialize_body(SerializingStream& s) const {
    BSplineCommon::serialize_body(s);
    s.pack("BSpline::coeffs", coeffs_);
  }

  BSpline::BSpline(DeserializingStream& s) : BSplineCommon(s) {
    s.unpack("BSpline::coeffs", coeffs_);
  }

  BSpline::BSpline(const MX& x, const std::vector<double>& knots,
          const std::vector<casadi_int>& offset,
          const std::vector<double>& coeffs,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const std::vector<casadi_int>& lookup_mode) :
          BSplineCommon(knots, offset, degree, m, lookup_mode), coeffs_(coeffs) {
    casadi_assert_dev(x.numel()==degree.size());
    set_dep(x);
    set_sparsity(Sparsity::dense(m, 1));
  }

  BSplineParametric::BSplineParametric(const MX& x,
          const MX& coeffs,
          const std::vector<double>& knots,
          const std::vector<casadi_int>& offset,
          const std::vector<casadi_int>& degree,
          casadi_int m,
          const std::vector<casadi_int>& lookup_mode) :
          BSplineCommon(knots, offset, degree, m, lookup_mode) {
    casadi_assert_dev(x.size1()==degree.size());
    set_dep(x, coeffs);
    set_sparsity(Sparsity::dense(m, 1));
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
                const std::vector< std::vector<double> >& knots,
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
      strides.push_back(strides.back()*(knots[i].size()-degree[i]-1));
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

  MX BSpline::create(const MX& x, const std::vector< std::vector<double> >& knots,
          const std::vector<double>& coeffs,
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
    std::vector<double> stacked;
    Interpolant::stack_grid(knots, offset, stacked);

    std::vector<casadi_int> mode =
      Interpolant::interpret_lookup_mode(lookup_mode, stacked, offset, degree, degree);

    if (do_inline_flag) {
      return do_inline(x, knots, coeffs, m, degree, mode);
    } else {
      return x->get_bspline(stacked, offset, coeffs, degree, m, mode);
    }
  }

  MX BSplineParametric::create(const MX& x,
          const MX& coeffs,
          const std::vector< std::vector<double> >& knots,
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
    std::vector<double> stacked;
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

  std::string BSplineParametric::disp(const std::vector<std::string>& arg) const {
    return "BSplineParametric(" + arg.at(0) + ")";
  }

  void BSpline::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_bspline(knots_, offset_, coeffs_, degree_, m_, lookup_mode_);
  }

  void BSplineParametric::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    res[0] = arg[0]->get_bspline(arg[0], knots_, offset_, degree_, m_, lookup_mode_);
  }

  MX BSpline::jac_cached() const {
    if (jac_cache_.is_empty()) {
      jac_cache_ = jac(dep(0), DM(coeffs_));
    }
    return jac_cache_;
  }

  MX BSplineParametric::jac_cached() const {
    if (jac_cache_.is_empty()) {
      jac_cache_ = jac(dep(0), dep(1));
    }
    return jac_cache_;
  }

  void BSplineCommon::ad_forward(const std::vector<std::vector<MX> >& fseed,
                          std::vector<std::vector<MX> >& fsens) const {
    MX J = jac_cached();

    for (casadi_int d=0; d<fsens.size(); ++d) {
      fsens[d][0] = mtimes(J, fseed[d][0]);
    }
  }

  void BSplineCommon::ad_reverse(const std::vector<std::vector<MX> >& aseed,
                          std::vector<std::vector<MX> >& asens) const {
    MX JT = jac_cached().T();
    for (casadi_int d=0; d<aseed.size(); ++d) {
      asens[d][0] += mtimes(JT, aseed[d][0]);
    }
  }

  int BSpline::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (!res[0]) return 0;

    casadi_clear(res[0], m_);
    casadi_nd_boor_eval(res[0], degree_.size(), get_ptr(knots_), get_ptr(offset_),
      get_ptr(degree_), get_ptr(strides_), get_ptr(coeffs_), m_, arg[0], get_ptr(lookup_mode_),
      iw, w);
    return 0;
  }

  int BSplineParametric::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    if (!res[0]) return 0;

    casadi_clear(res[0], m_);
    casadi_nd_boor_eval(res[0], degree_.size(), get_ptr(knots_), get_ptr(offset_),
      get_ptr(degree_), get_ptr(strides_), arg[1], m_, arg[0], get_ptr(lookup_mode_),
      iw, w);
    return 0;
  }

  void BSplineCommon::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg,
                      const std::vector<casadi_int>& res) const {
    casadi_int n_dims = offset_.size()-1;

    g.add_auxiliary(CodeGenerator::AUX_ND_BOOR_EVAL);
    g.add_auxiliary(CodeGenerator::AUX_FILL);
    g << g.clear(g.work(res[0], m_), m_) << "\n";

    // Input and output buffers
    g << "CASADI_PREFIX(nd_boor_eval)(" << g.work(res[0], m_) << "," << n_dims << ","
      << g.constant(knots_) << "," << g.constant(offset_) << "," <<  g.constant(degree_)
      << "," << g.constant(strides_) << "," << generate(g, arg) << "," << m_  << ","
      << g.work(arg[0], n_dims) << "," <<  g.constant(lookup_mode_) << ", iw, w);\n";
  }

  std::string BSpline::generate(CodeGenerator& g, const std::vector<casadi_int>& arg) const {
    return g.constant(coeffs_);
  }

  std::string BSplineParametric::generate(CodeGenerator& g,
                      const std::vector<casadi_int>& arg) const {
    return g.work(arg[1], dep(1).nnz());
  }

  DM BSpline::dual(const std::vector<double>& x,
          const std::vector< std::vector<double> >& knots,
          const std::vector<casadi_int>& degree,
          const Dict& opts) {

    std::vector<casadi_int> offset;
    std::vector<double> stacked;
    Interpolant::stack_grid(knots, offset, stacked);

    std::vector<std::string> lookup_mode;
    auto it = opts.find("lookup_mode");
    if (it!=opts.end()) lookup_mode = it->second;
    std::vector<casadi_int> lookup_mode_int =
      Interpolant::interpret_lookup_mode(lookup_mode, stacked, offset, degree, degree);

    casadi_int n_dims = degree.size();
    casadi_int N = x.size()/n_dims;
    casadi_assert_dev(N*n_dims==x.size());

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
        degree.size(), get_ptr(stacked), get_ptr(offset),
        get_ptr(degree), get_ptr(strides), get_ptr(x)+i*n_dims, get_ptr(lookup_mode_int),
        get_ptr(iw), get_ptr(w));
      data.insert(data.end(), contribution.begin(), contribution.begin()+nnz);
      col.insert(col.end(), nz.begin(), nz.begin()+nnz);
      row.insert(row.end(), nnz, i);
    }

    return DM(Sparsity::triplet(coeffs_size, N, col, row), data).T();
  }

} // namespace casadi
