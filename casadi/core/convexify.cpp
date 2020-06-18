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


#include "convexify.hpp"

using namespace std;
namespace casadi {

  std::string strategy_to_string(casadi_convexify_strategy_t s) {
    switch (s) {
      case CVX_REGULARIZE: return "regularize";
      case CVX_EIGEN_REFLECT: return "eigen-reflect";
      case CVX_EIGEN_CLIP: return "eigen-clip";
    }
    return "unknown";
  }

  Convexify::Convexify(const MX& H, const Dict& opts) {
    set_dep(H);
    set_sparsity(setup(convexify_data_, H.sparsity(), opts, false));
  }

  size_t Convexify::sz_iw() const {
    return convexify_data_.sz_iw;
  }

  size_t Convexify::sz_w() const {
    return convexify_data_.sz_w;
  }

  std::string Convexify::disp(const std::vector<std::string>& arg) const {
    return "convexify(" + arg.at(0) + ")";
  }

  void Convexify::eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const {
    Dict options;
    options["strategy"] = strategy_to_string(convexify_data_.config.strategy);
    options["margin"] = convexify_data_.config.margin;
    options["max_iter_eig"] = convexify_data_.config.max_iter_eig;
    res[0] = convexify(arg[0], options);
  }

  int Convexify::eval(const double** arg, double** res, casadi_int* iw, double* w) const {
    int ret = convexify_eval(&convexify_data_.config, arg[0], res[0], iw, w);
    casadi_assert(!ret, "Failure in convexification.");
    return 0;
  }

  std::string Convexify::generate(CodeGenerator& g,
    const ConvexifyData &d,
    const std::string& Hin, const std::string& Hout,
    const std::string& iw, const std::string& w) {
    g.local("cvx_config", "struct casadi_convexify_config");
    if (d.config.strategy==CVX_REGULARIZE) {
      g << "cvx_config.strategy = CVX_REGULARIZE;\n";
    } else if (d.config.strategy==CVX_EIGEN_CLIP) {
      g << "cvx_config.strategy = CVX_EIGEN_CLIP;\n";
    } else if (d.config.strategy==CVX_EIGEN_REFLECT) {
      g << "cvx_config.strategy = CVX_EIGEN_REFLECT;\n";
    }
    if (d.config.type_in==CVX_SYMM) {
      g << "cvx_config.type_in = CVX_SYMM;\n";
    } else if (d.config.type_in==CVX_TRIL) {
      g << "cvx_config.type_in = CVX_TRIL;\n";
    } else if (d.config.type_in==CVX_TRIU) {
      g << "cvx_config.type_in = CVX_TRIU;\n";
    }
    g << "cvx_config.Hsp = " << g.sparsity(d.Hsp) << ";\n";
    g << "cvx_config.Hrsp = " << g.sparsity(d.Hrsp) << ";\n";
    g << "cvx_config.margin = " << d.config.margin << ";\n";
    g << "cvx_config.Hsp_project = " << d.config.Hsp_project << ";\n";
    g << "cvx_config.scc_transform = " << d.config.scc_transform << ";\n";
    g << "cvx_config.scc_offset = " << g.constant(d.scc_offset) << ";\n";
    g << "cvx_config.scc_mapping = " << g.constant(d.scc_mapping) << ";\n";
    g << "cvx_config.scc_offset_size = " << d.scc_offset.size() << ";\n";
    g << "cvx_config.max_iter_eig = " << d.config.max_iter_eig << ";\n";
    g << "cvx_config.verbose = " << d.config.verbose << ";\n";
    return "convexify_eval(&cvx_config, " + Hin + "," + Hout + "," + iw + "," + "w)";
  }

  void Convexify::serialize_body(SerializingStream& s) const {
    MXNode::serialize_body(s);
    serialize(s, "", convexify_data_);
  }

  Convexify::Convexify(DeserializingStream& s) : MXNode(s) {
    deserialize(s, "", convexify_data_);
  }

  void Convexify::serialize(SerializingStream& s, const std::string& prefix,
      const ConvexifyData& d) {
    s.version(prefix + "Convexify", 1);
    s.pack(prefix + "Convexify::type_in", static_cast<int>(d.config.type_in));
    s.pack(prefix + "Convexify::strategy", static_cast<int>(d.config.strategy));
    s.pack(prefix + "Convexify::margin", d.config.margin);
    s.pack(prefix + "Convexify::max_iter_eig", d.config.max_iter_eig);
    s.pack(prefix + "Convexify::scc_offset", d.scc_offset);
    s.pack(prefix + "Convexify::scc_mapping", d.scc_mapping);
    s.pack(prefix + "Convexify::Hsp_project", d.config.Hsp_project);
    s.pack(prefix + "Convexify::scc_transform", d.config.scc_transform);
    s.pack(prefix + "Convexify::verbose", d.config.verbose);
    s.pack(prefix + "Convexify::Hsp", d.Hsp);
    s.pack(prefix + "Convexify::Hrsp", d.Hrsp);
  }

  void Convexify::deserialize(DeserializingStream& s, const std::string& prefix,
      ConvexifyData& d) {
    s.version(prefix + "Convexify", 1);
    int type_in;
    s.unpack(prefix + "Convexify::type_in", type_in);
    d.config.type_in = static_cast<casadi_convexify_type_in_t>(type_in);
    int strategy;
    s.unpack(prefix + "Convexify::strategy", strategy);
    d.config.strategy = static_cast<casadi_convexify_strategy_t>(strategy);
    s.unpack(prefix + "Convexify::margin", d.config.margin);
    s.unpack(prefix + "Convexify::max_iter_eig", d.config.max_iter_eig);
    s.unpack(prefix + "Convexify::scc_offset", d.scc_offset);
    s.unpack(prefix + "Convexify::scc_mapping", d.scc_mapping);
    s.unpack(prefix + "Convexify::Hsp_project", d.config.Hsp_project);
    s.unpack(prefix + "Convexify::scc_transform", d.config.scc_transform);
    s.unpack(prefix + "Convexify::verbose", d.config.verbose);
    s.unpack(prefix + "Convexify::Hsp", d.Hsp);
    s.unpack(prefix + "Convexify::Hrsp", d.Hrsp);


    d.config.scc_offset_size = d.scc_offset.size();

    // Set pointers
    d.config.Hsp = d.Hsp;
    d.config.Hrsp = d.Hrsp;;
    d.config.scc_offset = get_ptr(d.scc_offset);
    d.config.scc_mapping = get_ptr(d.scc_mapping);
  }

  void Convexify::generate(CodeGenerator& g,
                       const std::vector<casadi_int>& arg,
                       const std::vector<casadi_int>& res) const {
    std::string ret = g.convexify_eval(convexify_data_,
      g.work(arg[0], dep(0).nnz()), g.work(res[0], nnz()), "iw", "w");
    g << "if (" << ret << ") return 1;\n";
  }

  Sparsity Convexify::setup(ConvexifyData& d, const Sparsity& H, const Dict& opts, bool inplace) {
    // Validate and categorize matrix input sparsity
    casadi_assert(H.is_square(), "Convexify ");
    if (H.is_symmetric()) {
      d.config.type_in = CVX_SYMM;
    } else if (H.is_tril()) {
      d.config.type_in = CVX_TRIL;
    } else if (H.is_triu()) {
      d.config.type_in = CVX_TRIU;
    } else {
      casadi_error("Convexify operation requires symmetric or triangular input");
    }

    // Read options
    d.config.margin = 1e-7;
    d.config.max_iter_eig = 200;
    std::string strategy = "eigen-clip";
    d.config.verbose = false;

    for (auto&& op : opts) {
      if (op.first=="strategy") {
        strategy = op.second.to_string();
      } else if (op.first=="margin") {
        d.config.margin = op.second;
        casadi_assert(d.config.margin>=0, "Margin must be >=0");
      } else if (op.first=="max_iter_eig") {
        d.config.max_iter_eig = op.second;
      } else if (op.first=="verbose") {
        d.config.verbose = op.second;
      } else {
        casadi_error("Unknown option '" + op.first + "'.");
      }
    }

    // Interpret strategy
    if (strategy=="regularize") {
      d.config.strategy = CVX_REGULARIZE;
      casadi_assert(d.config.type_in==CVX_SYMM, "Only truly symmetric matrices supported");
    } else if (strategy=="eigen-reflect") {
      d.config.strategy = CVX_EIGEN_REFLECT;
    } else if (strategy=="eigen-clip") {
      d.config.strategy = CVX_EIGEN_CLIP;
    } else {
      casadi_error("Invalid convexify strategy. "
        "Choose from regularize|eigen-reflect|eigen-clip. Got '" + strategy + "'.");
    }

    d.Hrsp = H;

    d.config.scc_transform = 0;

    Sparsity Hrsp = H+H.T();

    casadi_int block_size = 0;
    Sparsity& Hsp = d.Hsp;
    if (d.config.strategy==CVX_EIGEN_REFLECT || d.config.strategy==CVX_EIGEN_CLIP) {
      // Uncover strongly connected components
      std::vector<casadi_int> scc_index;
      casadi_int scc_nb = Hrsp.scc(scc_index, d.scc_offset);

      // Represent Hessian as block-dense in permuted space
      std::vector<Sparsity> sp;
      for (casadi_int i=0;i<scc_nb;++i) {
        casadi_int block = d.scc_offset.at(i+1)-d.scc_offset.at(i);
        Sparsity stencil;
        if (d.config.type_in==CVX_SYMM) {
          stencil = Sparsity::dense(block, block);
        } else if (d.config.type_in==CVX_TRIL) {
          stencil = Sparsity::lower(block);
        } else {
          stencil = Sparsity::upper(block);
        }
        sp.push_back(stencil);
      }

      std::vector<casadi_int> ssc_perm = lookupvector(scc_index);
      std::vector<casadi_int> mapping_dummy;
      Hsp = diagcat(sp).sub(ssc_perm, ssc_perm, mapping_dummy);
      Hsp.sub(scc_index, scc_index, d.scc_mapping);

      // Find out size of maximum block
      for (casadi_int i=0;i<scc_nb;++i) {
        casadi_int block = d.scc_offset.at(i+1)-d.scc_offset.at(i);
        if (block>block_size) block_size = block;
      }

      d.config.scc_transform = d.scc_offset!=range(H.size1());

      if (d.config.verbose) casadi_message("Identified " + str(scc_nb) + " blocks "
        "with maximum size " + str(block_size) + ".");
    } else if (d.config.strategy==CVX_REGULARIZE) {
      Hsp = Hrsp + Sparsity::diag(H.size1());
    } else {
      Hsp = Hrsp;
    }
    d.config.Hsp_project = Hsp!=Hrsp;

    if (d.config.type_in==CVX_TRIL) {
      Hsp = Sparsity::tril(Hsp);
    } else if (d.config.type_in==CVX_TRIU) {
      Hsp = Sparsity::triu(Hsp);
    }

    d.sz_iw = 0;
    if (d.config.strategy==CVX_EIGEN_REFLECT || d.config.strategy==CVX_EIGEN_CLIP) {
       d.sz_iw = 1+3* d.config.max_iter_eig;
    }

    d.sz_w = 0;
    if (d.config.strategy==CVX_EIGEN_REFLECT || d.config.strategy==CVX_EIGEN_CLIP) {
      d.sz_w = max(block_size, 2*(block_size-1)*d.config.max_iter_eig);
      if (d.config.Hsp_project) d.sz_w = max(d.sz_w, Hsp.size1());
      if (d.config.scc_transform) d.sz_w += block_size*block_size;
      if (inplace) d.sz_w = max(d.sz_w, Hsp.size1()+d.Hrsp.nnz());
    }
    d.sz_w = Hsp.size1()+d.sz_w;

    d.config.scc_offset_size = d.scc_offset.size();

    // Set pointers
    d.config.Hsp = Hsp;
    d.config.Hrsp = d.Hrsp;
    d.config.scc_offset = get_ptr(d.scc_offset);
    d.config.scc_mapping = get_ptr(d.scc_mapping);

    return Hsp;
  }

} // namespace casadi
