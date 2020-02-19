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


#include "interpolant_impl.hpp"
#include "casadi_misc.hpp"
#include "mx_node.hpp"
#include "casadi_low.hpp"
#include <typeinfo>

using namespace std;
namespace casadi {

  bool has_interpolant(const string& name) {
    return Interpolant::has_plugin(name);
  }

  void load_interpolant(const string& name) {
    Interpolant::load_plugin(name);
  }

  string doc_interpolant(const string& name) {
    return Interpolant::getPlugin(name).doc;
  }

  void Interpolant::stack_grid(const std::vector< std::vector<double> >& grid,
    std::vector<casadi_int>& offset, std::vector<double>& stacked) {

    // Get offset for each input dimension
    offset.clear();
    offset.reserve(grid.size()+1);
    offset.push_back(0);
    for (auto&& g : grid) offset.push_back(offset.back()+g.size());

    // Stack input grids
    stacked.clear();
    stacked.reserve(offset.back());
    for (auto&& g : grid) stacked.insert(stacked.end(), g.begin(), g.end());
  }

  void Interpolant::stack_grid(const std::vector< MX> & grid,
    std::vector<casadi_int>& offset, MX& stacked) {

    // Get offset for each input dimension
    offset.clear();
    offset.reserve(grid.size()+1);
    offset.push_back(0);
    for (auto&& g : grid) offset.push_back(offset.back()+g.numel());

    // Stack input grids
    stacked = vertcat(grid);
  }

  void Interpolant::check_grid(const std::vector<MX> & grid) {
    // Dimension at least 1
    casadi_assert(!grid.empty(), "At least one dimension required");

    // Grid must be strictly increasing
    for (const auto& e : grid) {
      casadi_assert(e.is_dense() && e.is_vector(), "Grid component must be dense vector.");
      casadi_assert(e.numel()>=2, "Need at least two grid points for every input");
      //casadi_assert(is_increasing(g), "Gridpoints must be strictly increasing");
      //casadi_assert(is_regular(g), "Gridpoints must be regular");
      //casadi_assert(g.size()>=2, "Need at least two grid points for every input");
    }
  }

  std::vector< std::vector<double> > Interpolant::parse_grid(const std::vector< DM >& grid) {
    std::vector< std::vector<double> > ret;
    for (const auto& e : grid) {
      ret.push_back(parse_grid(e));
    }
    return ret;
  }

  std::vector<double> Interpolant::parse_grid(const DM& grid) {
    casadi_assert_dev(grid.is_vector());
    casadi_assert_dev(grid.is_dense());
    return grid.nonzeros();
  }

  std::vector<double> Interpolant::meshgrid(const std::vector< std::vector<double> >& grid) {
    std::vector<casadi_int> cnts(grid.size()+1, 0);
    std::vector<casadi_int> sizes(grid.size(), 0);
    for (casadi_int k=0;k<grid.size();++k) sizes[k]= grid[k].size();

    casadi_int total_iter = 1;
    for (casadi_int k=0;k<grid.size();++k) total_iter*= sizes[k];

    casadi_int n_dims = grid.size();

    std::vector<double> ret(total_iter*n_dims);
    for (casadi_int i=0;i<total_iter;++i) {

      for (casadi_int j=0;j<grid.size();++j) {
        ret[i*n_dims+j] = grid[j][cnts[j]];
      }

      cnts[0]++;
      casadi_int j = 0;
      while (j<n_dims && cnts[j]==sizes[j]) {
        cnts[j] = 0;
        j++;
        cnts[j]++;
      }

    }

    return ret;
  }

  casadi_int Interpolant::coeff_size() const  {
    return coeff_size(offset_, m_);
  }

  casadi_int Interpolant::coeff_size(const std::vector<casadi_int>& offset, casadi_int m) {
    casadi_int ret = 1;
    for (casadi_int k=0;k<offset.size()-1;++k) {
      ret *= offset[k+1]-offset[k];
    }
    return m*ret;
  }

  Function interpolant(const std::string& name,
                       const std::string& solver,
                       const std::vector<MX>& grid,
                       const MX& values,
                       const Dict& opts) {
    casadi_assert(values.is_dense() && values.is_vector(), "Values must be a dense vector.");

    MX values_p = values;
    if (values_p->type_==MXNode::MX_NONE) values_p = evalf(values_p);

    std::vector<MX> grid_p = grid;
    for (auto& e : grid_p) {
      if (e->type_==MXNode::MX_NONE) e = evalf(e);
    }

    Interpolant::check_grid(grid);
    // Get offset for each input dimension
    vector<casadi_int> offset;
    // Stack input grids
    MX stacked;

      // Consistency check, number of elements
    casadi_uint nel=1;
      for (auto&& g : grid) nel *= g.numel();
      casadi_assert(values.numel() % nel== 0,
        "Inconsistent number of elements. Must be a multiple of " +
        str(nel) + ", but got " + str(values.numel()) + " instead.");

    Interpolant::stack_grid(grid, offset, stacked);

    if (stacked->type_==MXNode::MX_NONE) stacked = evalf(stacked);

    casadi_int m = values.numel()/nel;
    return Interpolant::construct(solver, name, stacked, offset, values, m, opts);
  }

  Function Interpolant::construct(const std::string& solver,
                    const std::string& name,
                    const MX& grid,
                    const std::vector<casadi_int>& offset,
                    const MX& values,
                    casadi_int m,
                    const Dict& opts) {
    bool do_inline = false;
    Dict options = extract_from_dict(opts, "inline", do_inline);
    if (do_inline && !Interpolant::getPlugin(solver).exposed.do_inline) {
      options["inline"] = true;
      do_inline = false;
    }
    if (do_inline && Interpolant::getPlugin(solver).exposed.do_inline) {
      return Interpolant::getPlugin(solver).exposed.
        do_inline(name, grid, offset, values, m, options);
    } else {
      return Function::create(Interpolant::getPlugin(solver)
              .creator(name, grid, offset, values, m), options);
    }
  }

  Interpolant::
  Interpolant(const std::string& name,
              const MX& grid,
              const std::vector<casadi_int>& offset,
              const MX& values,
              casadi_int m)
              : FunctionInternal(name), m_(m), offset_(offset), grid_(grid), values_(values) {


    if (values_->type_ & MXNode::MX_PAR) values_ = values_->get_gate();
    if (grid_->type_ & MXNode::MX_PAR) grid_ = grid_->get_gate();

    // Number of grid points
    ndim_ = offset_.size()-1;
  }

  Interpolant::~Interpolant() {
  }

  char Interpolant::type() const {
    return values_->type_ | grid_->type_;
  }

  Sparsity Interpolant::get_sparsity_in(casadi_int i) {
    casadi_assert_dev(i==0);
    return Sparsity::dense(ndim_, batch_x_);
  }

  Sparsity Interpolant::get_sparsity_out(casadi_int i) {
    casadi_assert_dev(i==0);
    return Sparsity::dense(m_, batch_x_);
  }

  std::string Interpolant::get_name_in(casadi_int i) {
    if (i==0) return "x";
    casadi_assert_dev(false);
  }

  std::string Interpolant::get_name_out(casadi_int i) {
    casadi_assert_dev(i==0);
    return "f";
  }

  std::map<std::string, Interpolant::Plugin> Interpolant::solvers_;

  const std::string Interpolant::infix_ = "interpolant";

  const Options Interpolant::options_
  = {{&FunctionInternal::options_},
     {{"lookup_mode",
       {OT_STRINGVECTOR,
        "Specifies, for each grid dimenion, the lookup algorithm used to find the correct index. "
        "'linear' uses a for-loop + break; (default when #knots<=100), "
        "'exact' uses floored division (only for uniform grids), "
        "'binary' uses a binary search. (default when #knots>100)."}},
      {"inline",
       {OT_BOOL,
        "Implement the lookup table in MX primitives. "
        "Useful when you need derivatives with respect to grid and/or coefficients. "
        "Such derivatives are fundamentally dense, so use with caution."}},
      {"batch_x",
       {OT_INT,
        "Evaluate a batch of different inputs at once (default 1)."}}
     }
  };

  void Interpolant::init(const Dict& opts) {

    batch_x_ = 1;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="lookup_mode") {
        lookup_modes_ = op.second;
      } else if (op.first=="batch_x") {
        batch_x_ = op.second;
      }
    }

    // Call the base class initializer
    FunctionInternal::init(opts);

    // Needed by casadi_interpn
    alloc_w(ndim_, true);
    alloc_iw(2*ndim_, true);
  }


  void Interpolant::finalize() {

    try {
      grid_ptr_ = &static_cast<const std::vector<double> &>(grid_);
    } catch (...) {
      grid_vec_ = static_cast<DM>(grid_).nonzeros();
      grid_ptr_ = &grid_vec_;
    }
    try {
      values_ptr_ = &static_cast<const std::vector<double> &>(values_);
    } catch (...) {
      values_vec_ = static_cast<DM>(values_vec_).nonzeros();
      values_ptr_ = &values_vec_;
    }

    // Recursive call
    FunctionInternal::finalize();
  }

/**
  std::vector<casadi_int> Interpolant::interpret_lookup_mode(
      const std::vector<std::string>& modes, const std::vector<double>& knots,
      const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& margin_left, const std::vector<casadi_int>& margin_right) {
    casadi_assert_dev(modes.empty() || modes.size()==offset.size()-1);

    std::vector<casadi_int> ret;
    for (casadi_int i=0;i<offset.size()-1;++i) {
      casadi_int n = offset[i+1]-offset[i];
      ret.push_back(Low::interpret_lookup_mode(modes.empty() ? "auto": modes[i], n));
    }

    for (casadi_int i=0;i<offset.size()-1;++i) {
      if (ret[i]==LOOKUP_EXACT) {
        if (!knots.empty()) {
          casadi_int m_left  = margin_left.empty() ? 0 : margin_left[i];
          casadi_int m_right = margin_right.empty() ? 0 : margin_right[i];

          std::vector<double> grid(
              knots.begin()+offset[i]+m_left,
              knots.begin()+offset[i+1]-m_right);
          casadi_assert_dev(is_increasing(grid) && is_equally_spaced(grid));
        }
      }
    }
    return ret;
  }
*/
  std::vector<casadi_int> Interpolant::interpret_lookup_mode(
      const std::vector<std::string>& modes, const MX& knots,
      const std::vector<casadi_int>& offset,
      const std::vector<casadi_int>& margin_left, const std::vector<casadi_int>& margin_right) {
    casadi_assert_dev(modes.empty() || modes.size()==offset.size()-1);

    std::vector<casadi_int> ret;
    for (casadi_int i=0;i<offset.size()-1;++i) {
      casadi_int n = offset[i+1]-offset[i];
      ret.push_back(Low::interpret_lookup_mode(modes.empty() ? "auto": modes[i], n));
    }

    return ret;
  }

  void Interpolant::serialize_body(SerializingStream &s) const {
    FunctionInternal::serialize_body(s);
    s.version("Interpolant", 2);
    s.pack("Interpolant::ndim", ndim_);
    s.pack("Interpolant::m", m_);
    s.pack("Interpolant::grid", grid_);
    s.pack("Interpolant::offset", offset_);
    s.pack("Interpolant::values", values_);
    s.pack("Interpolant::lookup_modes", lookup_modes_);
    s.pack("Interpolant::batch_x", batch_x_);
  }

  void Interpolant::serialize_type(SerializingStream &s) const {
    FunctionInternal::serialize_type(s);
    PluginInterface<Interpolant>::serialize_type(s);
  }

  ProtoFunction* Interpolant::deserialize(DeserializingStream& s) {
    return PluginInterface<Interpolant>::deserialize(s);
  }

  Interpolant::Interpolant(DeserializingStream & s) : FunctionInternal(s) {
    int version = s.version("Interpolant", 1, 2);
    s.unpack("Interpolant::ndim", ndim_);
    s.unpack("Interpolant::m", m_);
    s.unpack("Interpolant::grid", grid_);
    s.unpack("Interpolant::offset", offset_);
    s.unpack("Interpolant::values", values_);
    s.unpack("Interpolant::lookup_modes", lookup_modes_);
    if (version==1) {
      batch_x_ = 1;
    } else {
      s.unpack("Interpolant::batch_x", batch_x_);
    }
  }

  void Interpolant::codegen_declarations(CodeGenerator& g) const {
    if (values_->type_ & MXNode::MX_PAR) {
      g.define_rw_double(values_.get(), values_.nnz());
    }
    if (grid_->type_ & MXNode::MX_PAR) {
      g.define_rw_double(grid_.get(), grid_.nnz());
    }
  }

  std::string Interpolant::codegen_values(CodeGenerator& g) const {
    return (values_->type_ & MXNode::MX_PAR) ? g.rw_double(values_.get()) : g.constant(*values_ptr_);
  }

  std::string Interpolant::codegen_grid(CodeGenerator& g) const {
    return (grid_->type_ & MXNode::MX_PAR) ? g.rw_double(grid_.get()) : g.constant(*grid_ptr_);
  }

} // namespace casadi
