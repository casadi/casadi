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


#include "linear_interpolant.hpp"

using namespace std;
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
  LinearInterpolant(const string& name,
                    const std::vector<double>& grid,
                    const std::vector<casadi_int>& offset,
                    const vector<double>& values,
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
    if (res[0]) {
      const double* values = has_parametric_values() ? arg[arg_values()] : get_ptr(values_);
      const double* grid = has_parametric_grid() ? arg[arg_grid()] : get_ptr(grid_);
      casadi_interpn(res[0], ndim_, grid, get_ptr(offset_),
                    values, arg[0], get_ptr(lookup_mode_), m_, iw, w);
    }
    return 0;
  }

  void LinearInterpolant::codegen_body(CodeGenerator& g) const {
    std::string values = has_parametric_values() ? g.arg(arg_values()) : g.constant(values_);
    std::string grid = has_parametric_grid() ? g.arg(arg_grid()) : g.constant(grid_);
    g << "  if (res[0]) {\n"
      << "    " << g.interpn("res[0]", ndim_, grid, g.constant(offset_),
      values, "arg[0]", g.constant(lookup_mode_), m_,  "iw", "w") << "\n"
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
                      values, arg[0], get_ptr(m->lookup_mode_), m->m_, iw, w);
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
      "arg[0]", g.constant(m->lookup_mode_), m->m_, "iw", "w") << "\n";
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

} // namespace casadi
