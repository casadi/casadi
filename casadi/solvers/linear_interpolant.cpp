/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "linear_interpolant.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTERPOLANT_LINEAR_EXPORT
  casadi_register_interpolant_linear(Interpolant::Plugin* plugin) {
    plugin->creator = LinearInterpolant::creator;
    plugin->name = "linear";
    plugin->doc = LinearInterpolant::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &LinearInterpolant::options_;
    plugin->deserialize = &LinearInterpolant::deserialize;
    plugin->exposed.do_inline = &LinearInterpolant::do_inline;
    return 0;
  }

  extern "C"
  void CASADI_INTERPOLANT_LINEAR_EXPORT casadi_load_interpolant_linear() {
    Interpolant::registerPlugin(casadi_register_interpolant_linear);
  }

  const Options LinearInterpolant::options_
  = {{&Interpolant::options_},
     {{"lookup_mode",
       {OT_STRINGVECTOR,
        "Sets, for each grid dimenion, the lookup algorithm used to find the correct index. "
        "'linear' uses a for-loop + break; "
        "'exact' uses floored division (only for uniform grids)."}}
     }
  };

  LinearInterpolant::
  LinearInterpolant(const std::string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const std::vector<double>& values,
                    casadi_int m)
                    : Interpolant(name, grid, offset, values, m) {
  }

  LinearInterpolant::~LinearInterpolant() {
    clear_mem();
  }

  LinearInterpolantJac::~LinearInterpolantJac() {
    clear_mem();
  }

  void LinearInterpolant::init(const Dict& opts) {
    // Call the base class initializer
    Interpolant::init(opts);

    lookup_mode_ = Interpolant::interpret_lookup_mode(lookup_modes_, grid_, offset_);

    // Needed by casadi_interpn
    alloc_w(ndim_, true);
    alloc_iw(2*ndim_, true);
  }

  int LinearInterpolant::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    setup(mem, arg, res, iw, w);
    if (res[0]) {
      const double* values = has_parametric_values() ? arg[arg_values()] : get_ptr(values_);
      const double* grid = has_parametric_grid() ? arg[arg_grid()] : get_ptr(grid_);
      for (casadi_int i=0;i<ndim_;++i) {
        if (extrapolation_mode_[i]==INTERP_EXTRAPOLATION_ERROR) {
          casadi_assert(arg[0][i]>=grid[offset_[i]] && arg[0][i]<=grid[offset_[i+1]-1],
                        "Extrapolation explicitly forbidden (option extrapolation_mode 'error). "
                        "Linear interpolation on axis i=" + str(i) + ": "
                        "Grid runs from " + str(grid[offset_[i]]) + " to " + str(grid[offset_[i+1]-1]) + ", "
                        "but got " + str(arg[0][i]) + " as value instead.");
        }
      }

      casadi_interpn(res[0], ndim_, grid, get_ptr(offset_),
                    values, arg[0], get_ptr(lookup_mode_), get_ptr(extrapolation_mode_), m_, iw, w);
    }
    return 0;
  }

class LinearInterpolantIntervalPropagator : public FunctionInternal {
public:
  Function orig_;

  /** \brief  Constructor */
  LinearInterpolantIntervalPropagator(const std::string& name, const Function& orig) :
    FunctionInternal(name), orig_(orig) {
  }

  /** \brief  Destructor */
  ~LinearInterpolantIntervalPropagator() override {
    clear_mem();
  }

  void init(const Dict& opts) override {
    // Call the initialization method of the base class
    FunctionInternal::init(opts);

    alloc_w(orig_.sz_w());
    alloc_iw(orig_.sz_iw());
  }

  /** \brief Get type name */
  std::string class_name() const override {return "LinearInterpolantIntervalPropagator";}

  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return orig_.n_in()*2;}
  size_t get_n_out() override { return orig_.n_out()*2;}

  /** \brief Sparsities of function inputs and outputs */
  Sparsity get_sparsity_in(casadi_int i) override {
    return orig_.sparsity_in(i % orig_.n_in());
  }
  Sparsity get_sparsity_out(casadi_int i) override {
    return orig_.sparsity_out(i % orig_.n_out());
  }

  ///@{
  /** \brief Names of function input and outputs */
  std::string get_name_in(casadi_int i) override {
    return orig_.name_in(i % orig_.n_in()) + "_L";
  }
  std::string get_name_out(casadi_int i) override {
    return orig_.name_out(i % orig_.n_out()) + "_R";
  }
  /// @}

  /** \brief  Evaluate numerically */
  virtual int eval(const double** arg, double** res,
    casadi_int* iw, double* w, void* mem) const override {
    LinearInterpolant* orig = orig_.get<LinearInterpolant>();

    casadi_int ndim = orig->ndim_;
    casadi_int m = orig->m_;

    int orig_mem = orig->checkout();

    double* orig_w = w;
    casadi_int* orig_iw = iw;

    const std::vector<double>& grid = orig->grid_;
    const std::vector<casadi_int>& offset = orig->offset_;
    const std::vector<double>& values = orig->values_;

    std::vector<casadi_int> strides(ndim, m);
    for (casadi_int i=1;i<ndim;++i) {
      strides[i] = strides[i-1]*(offset[i]-offset[i-1]);
    }

    // Pairs of these define a subset of the grid that lies in the interior, per coordinate
    casadi_int interior_numel = 1;
    std::vector<casadi_int> interior_start(ndim), interior_size(ndim);
    for (casadi_int i=0;i<ndim;++i) {
      for (casadi_int j=0;j<offset[i+1]-offset[i];++j) {
      if (grid[offset[i]+j]>arg[0][i] && grid[offset[i]+j]<arg[1][i]) {
          if (interior_size[i]==0) {
            interior_start[i] = j;
          }
          interior_size[i]++;
        }
      }
      interior_numel*=interior_size[i];
    }

    // General approach:
    //   The lower and upper interval bounds X_L[i] and X_R[i] define a hypercube.
    //   That hypercube together with the knots defines a set of critical points.
    //   Get the min/max of the interpolating surface at all critical points.
    //
    // Critical points are given by the cartesian product of
    //      [X_L[i], X_R[i], knots[i] that lie in (X_L[i], X_R[i])]
    //
    // At knots, the interpolant values are given directly by the values array.
    // For efficiency, we treat these 'interior' critical points separately,
    // and consider only the remaining 'surfacial' critical points for interpolation


    // current multi-index, initialized to start[]
    std::vector<casadi_int> idx = interior_start;

    // Linear index into values
    casadi_int lin = 0;
    for (casadi_int i = 0; i < ndim; ++i) {
      lin += interior_start[i] * strides[i];
    }

    // Initialize the L and R outputs
    for (casadi_int k=0;k<m;++k) {
      res[0][k] = inf;
      res[1][k] = -inf;
    }

    if (interior_numel) {
      for (casadi_int k=0;k<m;++k) {
        res[0][k] = fmin(res[0][k], values[lin+k]);
        res[1][k] = fmax(res[1][k], values[lin+k]);
      }
    }

    // Visit cartesian product of all interior grids
    while (true) {
      bool advanced = false;
      for (casadi_int i=0; i < ndim; ++i) {
        // Can we advance axis i by one and stay inside?
        if (idx[i] + 1 < interior_start[i]+interior_size[i]) {
          ++idx[i];
          lin += strides[i];        // take one step along axis i

          // Visit values[lin]
          for (casadi_int k=0;k<m;++k) {
            res[0][k] = fmin(res[0][k], values[lin+k]);
            res[1][k] = fmax(res[1][k], values[lin+k]);
          }
          advanced = true;
          break;                    // done with this iteration
        } else {
          // wrap axis i back to start: adjust from CURRENT idx[i] (== stop_i-1) to start
          lin += (interior_start[i] - idx[i]) * strides[i];
          idx[i] = interior_start[i];
          // continue to try slower axis (i+1)
        }
      }
      if (!advanced) break;
    }

    std::vector< std::vector<double> > critical_points(ndim);

    for (casadi_int i=0;i<ndim;++i) {
      critical_points[i].push_back(arg[0][i]);
      for (casadi_int j=interior_start[i];j<interior_start[i]+interior_size[i];++j) {
        critical_points[i].push_back(grid[offset[i]+j]);
      }
      critical_points[i].push_back(arg[1][i]);
    }

    // Odometer state
    std::vector<casadi_int> pos(ndim, 0);
    std::vector<double> point(ndim);
    std::vector<casadi_int> len(ndim);
    for (casadi_int i = 0; i < ndim; ++i) len[i] = critical_points[i].size();

    // Initialize first point (all X_L)
    for (casadi_int i = 0; i < ndim; ++i) point[i] = critical_points[i][pos[i]];

    std::vector<const double *> orig_arg(1);
    std::vector<double *> orig_res(1);

    std::vector<double> orig_output(m);

    orig_arg[0] = get_ptr(point);
    orig_res[0] = get_ptr(orig_output);

    // Iterate product; skip the single “all-interior” tuple
    while (true) {
      bool all_interior = true;
      for (casadi_int d = 0; d < ndim; ++d) {
        if (pos[d] == 0 || pos[d] + 1 == len[d]) { all_interior = false; break; }
      }
      if (!all_interior) {
        orig_(get_ptr(orig_arg), get_ptr(orig_res), orig_iw, orig_w, orig_mem);
        // Process the interpolant outputs
        for (casadi_int k=0;k<m;++k) {
          res[0][k] = fmin(res[0][k], orig_output[k]);
          res[1][k] = fmax(res[1][k], orig_output[k]);
        }
      }

      // advance odometer (axis 0 fastest; flip if you prefer)
      bool advanced = false;
      for (casadi_int i = 0; i < ndim; ++i) {
        if (pos[i] + 1 < len[i]) {
          ++pos[i];
          point[i] = critical_points[i][pos[i]];
          // reset faster axes 0..i-1
          for (casadi_int j = 0; j < i; ++j) {
            pos[j] = 0;
            point[j] = critical_points[j][pos[j]];
          }
          advanced = true;
          break;
        }
      }
      if (!advanced) break; // done
    }

    orig->release(orig_mem);

    return 0;
  }
};

Function LinearInterpolant::get_interval_propagator(const Dict& opts) const {
    Function ret;
    ret.own(
      new LinearInterpolantIntervalPropagator(name_ + "_interval_propagator", self()));
    ret->construct(opts);

    return ret;
  }

  void LinearInterpolant::codegen_body(CodeGenerator& g) const {
    std::string values = has_parametric_values() ? g.arg(arg_values()) : g.constant(values_);
    std::string grid = has_parametric_grid() ? g.arg(arg_grid()) : g.constant(grid_);
    g << "  if (res[0]) {\n"
      << "    " << g.interpn("res[0]", ndim_, grid, g.constant(offset_),
      values, "arg[0]", g.constant(lookup_mode_), g.constant(extrapolation_mode_), m_,  "iw", "w") << "\n"
      << "  }\n";
  }

  Function LinearInterpolant::
  get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    Function ret;
    ret.own(new LinearInterpolantJac(name));
    ret->construct(opts);
    return ret;
  }

  Function LinearInterpolantJac::
  get_jacobian(const std::string& name,
                  const std::vector<std::string>& inames,
                  const std::vector<std::string>& onames,
                  const Dict& opts) const {
    std::vector<MX> args = mx_in();
    std::vector<MX> res(n_out_);
    for (casadi_int i=0;i<n_out_;++i)
      res[i] = DM(size1_out(i), size2_out(i));
    Function f("f", args, res);

    return f->get_jacobian(name, inames, onames, Dict());
  }

  void LinearInterpolantJac::init(const Dict& opts) {
    // Call the base class initializer
    FunctionInternal::init(opts);

    // Needed by casadi_interpn
    auto m = derivative_of_.get<LinearInterpolant>();
    alloc_w(2*m->ndim_ + m->m_, true);
    alloc_iw(2*m->ndim_, true);
  }

  int LinearInterpolantJac::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = derivative_of_.get<LinearInterpolant>();

    const double* values = has_parametric_values() ? arg[m->arg_values()] : get_ptr(m->values_);
    const double* grid = has_parametric_grid() ? arg[m->arg_grid()] : get_ptr(m->grid_);

    casadi_interpn_grad(res[0], m->ndim_, grid, get_ptr(m->offset_),
                      values, arg[0], get_ptr(m->lookup_mode_), get_ptr(m->extrapolation_mode_),
                      m->m_, iw, w);
    return 0;
  }

  bool LinearInterpolantJac::has_parametric_values() const {
    auto m = derivative_of_.get<LinearInterpolant>();
    return m->has_parametric_values();
  }

  bool LinearInterpolantJac::has_parametric_grid() const {
    auto m = derivative_of_.get<LinearInterpolant>();
    return m->has_parametric_grid();
  }

  void LinearInterpolantJac::codegen_body(CodeGenerator& g) const {

    auto m = derivative_of_.get<LinearInterpolant>();
    std::string values = has_parametric_values() ? g.arg(m->arg_values()) : g.constant(m->values_);
    std::string grid = has_parametric_grid() ? g.arg(m->arg_grid()) : g.constant(m->grid_);

    g << "  " << g.interpn_grad("res[0]", m->ndim_,
      grid, g.constant(m->offset_), values,
      "arg[0]", g.constant(m->lookup_mode_), g.constant(m->extrapolation_mode_),
      m->m_, "iw", "w") << "\n";
  }


  LinearInterpolant::LinearInterpolant(DeserializingStream& s) : Interpolant(s) {
    s.unpack("LinearInterpolant::lookup_mode", lookup_mode_);
  }

  ProtoFunction* LinearInterpolant::deserialize(DeserializingStream& s) {
    s.version("LinearInterpolant", 1);
    char type;
    s.unpack("LinearInterpolant::type", type);
    switch (type) {
      case 'f': return new LinearInterpolant(s);
      case 'j': return new LinearInterpolantJac(s);
      default:
        casadi_error("LinearInterpolant::deserialize error");
    }
  }

  void LinearInterpolant::serialize_body(SerializingStream &s) const {
    Interpolant::serialize_body(s);
    s.pack("LinearInterpolant::lookup_mode", lookup_mode_);
  }

  void LinearInterpolant::serialize_type(SerializingStream &s) const {
    Interpolant::serialize_type(s);
    s.version("LinearInterpolant", 1);
    s.pack("LinearInterpolant::type", 'f');
  }

  void LinearInterpolantJac::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    auto m = derivative_of_.get<LinearInterpolant>();
    m->PluginInterface<Interpolant>::serialize_type(s);
    s.version("LinearInterpolant", 1);
    s.pack("LinearInterpolant::type", 'j');
  }


  Function LinearInterpolant::do_inline(const std::string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const std::vector<double>& values,
                    casadi_int m,
                    const Dict& opts) {

    // Number of grid points
    casadi_int ndim = offset.size()-1;

    MX x = MX::sym("x", ndim);
    MX g;
    if (grid.empty()) {
      g = MX::sym("g", offset.back());
    } else {
      g = MX(DM(grid));
    }
    MX c;
    if (values.empty()) {
      c = MX::sym("c", Interpolant::coeff_size(offset, m));
    } else {
      c = MX(DM(values));
    }

    MX f = MX::interpn_linear(vertsplit(g, offset), c, vertsplit(x), opts);

    std::vector<MX> args = {x};
    std::vector<std::string> arg_names = {"x"};
    if (grid.empty()) {
      args.push_back(g);
      arg_names.push_back("g");
    }
    if (values.empty()) {
      args.push_back(c);
      arg_names.push_back("c");
    }

    return Function(name, args, {f.T()}, arg_names, {"f"});
  }

} // namespace casadi
